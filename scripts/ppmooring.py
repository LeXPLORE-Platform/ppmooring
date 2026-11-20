# -*- coding: utf-8 -*-
from __future__ import annotations

#### Import Packages
import logging
import sys
import os

import netCDF4
import pandas as pd
import numpy as np
import fnmatch
import datetime as dt
from datetime import datetime, timezone
from copy import deepcopy
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from numpy.typing import ArrayLike

# Import configuration management
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config import get_config

logger = logging.getLogger(__name__)


#### Time grid helpers

def datetimes_to_microseconds(values) -> np.ndarray:
    """Convert naive datetime objects to integer microseconds since 1970-01-01."""
    return np.array(values, dtype='datetime64[us]').astype(np.int64)


def match_to_grid(data_us: np.ndarray, grid_us: np.ndarray, tolerance_us: int = 150_000_000):
    """Select the measurements that sit on a time grid.

    A measurement is kept when it falls exactly on a grid point, or when it is less
    than ``tolerance_us`` (2.5 minutes) away from its nearest grid point and closer to
    that grid point than the measurements just before and just after it.

    Parameters
    ----------
    data_us : int array, sorted, without duplicates
        Measurement times in microseconds.
    grid_us : int array, sorted
        Grid times in microseconds.

    Returns
    -------
    grid_mask : bool array, True for the grid points that receive a measurement
    data_mask : bool array, True for the measurements that are kept
    """
    data_us = np.asarray(data_us, dtype=np.int64)
    grid_us = np.asarray(grid_us, dtype=np.int64)
    n_data, n_grid = len(data_us), len(grid_us)
    grid_mask = np.full(n_grid, False)
    data_mask = np.full(n_data, False)
    if n_data == 0 or n_grid == 0:
        return grid_mask, data_mask

    # Nearest grid point of each measurement; the earlier one when two are equally close
    pos = np.searchsorted(grid_us, data_us, side='left')
    before = np.clip(pos - 1, 0, n_grid - 1)
    after = np.clip(pos, 0, n_grid - 1)
    dist_before = np.abs(data_us - grid_us[before])
    dist_after = np.abs(grid_us[after] - data_us)
    nearest = np.where(dist_before <= dist_after, before, after)
    dist = np.minimum(dist_before, dist_after)
    candidate = grid_us[nearest]

    exact = dist == 0
    near = (dist > 0) & (dist < tolerance_us)

    # The measurement must be closer to its grid point than its two neighbours are.
    # The first measurement has no previous neighbour, the last one no next neighbour.
    closer_than_previous = np.full(n_data, True)
    closer_than_next = np.full(n_data, True)
    if n_data > 1:
        closer_than_previous[1:] = dist[1:] < np.abs(candidate[1:] - data_us[:-1])
        closer_than_next[:-1] = dist[:-1] < np.abs(candidate[:-1] - data_us[1:])

    data_mask = exact | (near & closer_than_previous & closer_than_next)
    grid_mask[nearest[data_mask]] = True
    return grid_mask, data_mask



#### Quality control tests (QC v2)
#
# Independent test functions. Each one is pure: it takes data in and returns a
# boolean mask. No side effects, no re-initialization, no cascading.
#
# Flag bits:
#     bit 0 (1):  Range test
#     bit 1 (2):  Maintenance mask
#     bit 2 (4):  Climatological range
#     bit 3 (8):  Spike test
#     bit 4 (16): Gradient test
#     bit 5 (32): Flat-line test

# Result of the apply_* functions: bit flags, 0/1 quality, mask of each test
QCResult = Tuple[np.ndarray, np.ndarray, Dict[str, np.ndarray]]



# ── Flag bit constants ────────────────────────────────────────────────────────
FLAG_RANGE = 0          # bit 0
FLAG_MAINTENANCE = 1    # bit 1
FLAG_CLIMATOLOGICAL = 2 # bit 2
FLAG_SPIKE = 3          # bit 3
FLAG_GRADIENT = 4       # bit 4
FLAG_FLATLINE = 5       # bit 5

FLAG_NAMES = {
    FLAG_RANGE: "range",
    FLAG_MAINTENANCE: "maintenance",
    FLAG_CLIMATOLOGICAL: "climatological",
    FLAG_SPIKE: "spike",
    FLAG_GRADIENT: "gradient",
    FLAG_FLATLINE: "flatline",
}


def test_range(values: ArrayLike, vmin: float, vmax: float) -> np.ndarray:
    """
    Test 1: Flag values outside hard physical limits [vmin, vmax].
    
    Parameters
    ----------
    values : np.ndarray
        DO values in mgO2/L.
    vmin : float
        Minimum valid value.
    vmax : float
        Maximum valid value.
    
    Returns
    -------
    np.ndarray (bool)
        True where the test fails (value is out of range).
    """
    values = np.asarray(values, dtype=float)
    return (values < vmin) | (values > vmax) | np.isnan(values) | np.isinf(values)


def test_maintenance(timestamps: ArrayLike, maintenance_df: Optional[pd.DataFrame]) -> np.ndarray:
    """
    Test 2: Flag data points falling within maintenance windows.
    
    Parameters
    ----------
    timestamps : np.ndarray
        Unix timestamps of measurements.
    maintenance_df : pd.DataFrame
        DataFrame with 'debut' and 'end' columns (unix timestamps).
    
    Returns
    -------
    np.ndarray (bool)
        True where the test fails (point is within a maintenance window).
    """
    timestamps = np.asarray(timestamps, dtype=float)
    mask = np.zeros(len(timestamps), dtype=bool)
    
    if maintenance_df is None or len(maintenance_df) == 0:
        return mask
    
    for _, row in maintenance_df.iterrows():
        window_mask = (timestamps >= row["debut"]) & (timestamps <= row["end"])
        mask |= window_mask
    
    return mask


def test_climatological_range(values: ArrayLike, timestamps: ArrayLike,
                              border_val_df: Optional[pd.DataFrame], depth_cm: int) -> np.ndarray:
    """
    Test 3: Flag values outside month- and depth-specific climatological range.
    
    Uses the border_value_lower.csv file to determine expected range.
    Applies min_margin (0.95) and max_margin (1.20) as in the current code.
    
    Parameters
    ----------
    values : np.ndarray
        DO values in mgO2/L.
    timestamps : np.ndarray
        Unix timestamps.
    border_val_df : pd.DataFrame
        DataFrame from border_value_lower.csv with columns:
        month, min_depth (in meters), min, max.
    depth_cm : int
        Sensor depth in cm.
    
    Returns
    -------
    np.ndarray (bool)
        True where the test fails.
    """
    values = np.asarray(values, dtype=float)
    timestamps = np.asarray(timestamps, dtype=float)
    mask = np.zeros(len(values), dtype=bool)
    
    if border_val_df is None or len(border_val_df) == 0:
        return mask
    
    # Ensure numeric columns
    bv = border_val_df.copy()
    bv[["min_depth", "min", "max"]] = bv[["min_depth", "min", "max"]].apply(pd.to_numeric)
    # Convert min_depth from meters to cm for matching
    bv["min_depth_cm"] = bv["min_depth"] * 100
    
    # Get unique months in the data
    datetimes = np.array([datetime.utcfromtimestamp(ts) for ts in timestamps])
    months = np.array([dt.month for dt in datetimes])
    
    for month in np.unique(months):
        month_mask = months == month
        
        # Find matching row: month + depth
        rows = bv[(bv["month"] == month) & (bv["min_depth_cm"] == depth_cm)]
        if len(rows) == 0:
            # Fallback to depth=0 (same logic as current code)
            rows = bv[(bv["month"] == month) & (bv["min_depth_cm"] == 0)]
        if len(rows) == 0:
            continue
        
        row = rows.iloc[0]
        vo2_min = np.floor(row["min"] * 0.95)
        vo2_max = np.round(row["max"] * 1.20)
        
        # Flag values outside climatological range for this month
        month_fail = (values < vo2_min) | (values > vo2_max)
        mask |= (month_mask & month_fail)
    
    return mask


def test_spike(values: ArrayLike, threshold: float, window: int = 1) -> np.ndarray:
    """
    Test 4: Flag spikes — short deviations that return to baseline.
    
    A spike at index i is detected when:
    - |v[i] - mean(v[i-w], v[i+w])| > threshold  (point deviates from neighbors)
    - |v[i-w] - v[i+w]| < threshold                (neighbors agree — not a step change)
    
    This correctly distinguishes spikes from legitimate step changes.
    
    Parameters
    ----------
    values : np.ndarray
        DO values.
    threshold : float
        Maximum allowed deviation from neighbors (mgO2/L).
    window : int
        Number of neighbors on each side to consider.
    
    Returns
    -------
    np.ndarray (bool)
        True where spike is detected.
    """
    values = np.asarray(values, dtype=float)
    n = len(values)
    mask = np.zeros(n, dtype=bool)
    
    if n < 2 * window + 1:
        return mask
    
    for i in range(window, n - window):
        # Get neighbor values (skip NaN neighbors)
        before = values[i - window:i]
        after = values[i + 1:i + window + 1]
        
        valid_before = before[~np.isnan(before)]
        valid_after = after[~np.isnan(after)]
        
        if len(valid_before) == 0 or len(valid_after) == 0:
            continue
        
        neighbor_mean = 0.5 * (np.mean(valid_before) + np.mean(valid_after))
        neighbor_spread = abs(np.mean(valid_before) - np.mean(valid_after))
        
        deviation = abs(values[i] - neighbor_mean)
        
        # Spike: large deviation from neighbors, but neighbors agree with each other
        if deviation > threshold and neighbor_spread < threshold:
            mask[i] = True
    
    return mask


def test_gradient(values: ArrayLike, timestamps: ArrayLike, max_rate: float) -> np.ndarray:
    """
    Test 5: Flag points where rate of change exceeds physical maximum.
    
    Parameters
    ----------
    values : np.ndarray
        DO values in mgO2/L.
    timestamps : np.ndarray
        Unix timestamps.
    max_rate : float
        Maximum allowed |dO2/dt| in mgO2/L per hour.
    
    Returns
    -------
    np.ndarray (bool)
        True where gradient test fails. Both endpoints of an excessive
        gradient are flagged.
    """
    values = np.asarray(values, dtype=float)
    timestamps = np.asarray(timestamps, dtype=float)
    n = len(values)
    mask = np.zeros(n, dtype=bool)
    
    if n < 2:
        return mask
    
    # Compute rate of change in mgO2/L per hour
    dt_hours = np.diff(timestamps) / 3600.0
    dv = np.diff(values)
    
    # Avoid division by zero
    dt_hours[dt_hours == 0] = np.nan
    rate = np.abs(dv / dt_hours)
    
    # Flag both endpoints of excessive gradients
    excessive = rate > max_rate
    # Suppress NaN comparisons
    excessive[np.isnan(rate)] = False
    
    mask[:-1] |= excessive  # flag the point before
    mask[1:] |= excessive   # flag the point after
    
    return mask


def test_flatline(values: ArrayLike, tolerance: float, min_count: int) -> np.ndarray:
    """
    Test 6: Flag runs of consecutive near-identical values (sensor failure).
    
    Parameters
    ----------
    values : np.ndarray
        DO values.
    tolerance : float
        Values within this tolerance are considered "identical".
    min_count : int
        Minimum consecutive identical values to trigger the flag.
    
    Returns
    -------
    np.ndarray (bool)
        True where flat-line is detected.
    """
    values = np.asarray(values, dtype=float)
    n = len(values)
    mask = np.zeros(n, dtype=bool)
    
    if n < min_count:
        return mask
    
    run_start = 0
    for i in range(1, n):
        if abs(values[i] - values[run_start]) > tolerance or np.isnan(values[i]):
            # Check if the run was long enough
            run_length = i - run_start
            if run_length >= min_count:
                mask[run_start:i] = True
            run_start = i
    
    # Check final run
    run_length = n - run_start
    if run_length >= min_count:
        mask[run_start:n] = True
    
    return mask


# ── Helper functions ──────────────────────────────────────────────────────────

def apply_all_tests(values: ArrayLike, timestamps: ArrayLike, depth_cm: int,
                    maintenance_df: Optional[pd.DataFrame], border_val_df: Optional[pd.DataFrame],
                    config: dict) -> QCResult:
    """
    Run all 6 QC tests and return a bitfield integer array.
    
    Parameters
    ----------
    values : np.ndarray
        DO values in mgO2/L.
    timestamps : np.ndarray
        Unix timestamps.
    depth_cm : int
        Sensor depth in cm.
    maintenance_df : pd.DataFrame or None
        Maintenance windows (debut, end columns).
    border_val_df : pd.DataFrame or None
        Climatological bounds.
    config : dict
        QC configuration from qc_config.yaml['DO'].
    
    Returns
    -------
    qc_flags : np.ndarray (int)
        Bitfield integer array. Each bit represents a test.
    qual : np.ndarray (int)
        Simple 0/1 array for backward compatibility.
        1 if any test failed, 0 otherwise.
    flag_details : dict
        Per-test boolean masks for diagnostics.
    """
    n = len(values)
    qc_flags = np.zeros(n, dtype=np.int32)
    flag_details = {}
    
    # Test 1: Range
    cfg = config.get("range_test", {})
    mask = test_range(values, cfg.get("min", 0.0), cfg.get("max", 25.0))
    qc_flags[mask] |= (1 << FLAG_RANGE)
    flag_details["range"] = mask
    
    # Test 2: Maintenance
    mask = test_maintenance(timestamps, maintenance_df)
    qc_flags[mask] |= (1 << FLAG_MAINTENANCE)
    flag_details["maintenance"] = mask
    
    # Test 3: Climatological range
    mask = test_climatological_range(values, timestamps, border_val_df, depth_cm)
    qc_flags[mask] |= (1 << FLAG_CLIMATOLOGICAL)
    flag_details["climatological"] = mask
    
    # Test 4: Spike
    cfg = config.get("spike_test", {})
    mask = test_spike(values, cfg.get("threshold", 3.0), cfg.get("window", 1))
    qc_flags[mask] |= (1 << FLAG_SPIKE)
    flag_details["spike"] = mask
    
    # Test 5: Gradient
    cfg = config.get("gradient_test", {})
    mask = test_gradient(values, timestamps, cfg.get("max_rate", 2.0))
    qc_flags[mask] |= (1 << FLAG_GRADIENT)
    flag_details["gradient"] = mask
    
    # Test 6: Flat-line
    cfg = config.get("flatline_test", {})
    mask = test_flatline(values, cfg.get("tolerance", 0.001), cfg.get("min_count", 20))
    qc_flags[mask] |= (1 << FLAG_FLATLINE)
    flag_details["flatline"] = mask
    
    # Backward-compatible qual: 1 if any flag is set
    qual = (qc_flags > 0).astype(np.int32)
    
    return qc_flags, qual, flag_details


def decode_flags(qc_value: int) -> List[str]:
    """
    Decode a bitfield QC value into a list of test names that failed.
    
    Parameters
    ----------
    qc_value : int
        Bitfield QC value.
    
    Returns
    -------
    list of str
        Names of tests that failed.
    """
    failed = []
    for bit, name in FLAG_NAMES.items():
        if qc_value & (1 << bit):
            failed.append(name)
    return failed


def summarize_flags(qc_flags: np.ndarray) -> None:
    """
    Print a summary of how many points were flagged by each test.
    
    Parameters
    ----------
    qc_flags : np.ndarray (int)
        Bitfield QC array.
    """
    total = len(qc_flags)
    flagged = np.sum(qc_flags > 0)
    logger.info("  Total points: %d, flagged: %d (%.1f%%)", total, flagged, 100*flagged/total)
    
    for bit, name in FLAG_NAMES.items():
        count = np.sum((qc_flags & (1 << bit)) > 0)
        if count > 0:
            logger.info("    %-20s: %6d (%.1f%%)", name, count, 100*count/total)


# ── PAR-specific QC ──────────────────────────────────────────────────────────

def apply_par_tests(values: ArrayLike, timestamps: ArrayLike, depth_cm: int,
                    maintenance_df: Optional[pd.DataFrame], config: dict) -> QCResult:
    """
    Run PAR QC tests: range + maintenance only.
    
    Parameters
    ----------
    values : np.ndarray
        PAR values in µmol/m²/s.
    timestamps : np.ndarray
        Unix timestamps.
    depth_cm : int
        Sensor depth in cm.
    maintenance_df : pd.DataFrame or None
        Maintenance windows (debut, end columns).
    config : dict
        QC configuration from qc_config.yaml['PAR'].
    
    Returns
    -------
    qc_flags : np.ndarray (int)
        Bitfield integer array.
    qual : np.ndarray (int)
        Simple 0/1 array for backward compatibility.
    flag_details : dict
        Per-test boolean masks.
    """
    n = len(values)
    qc_flags = np.zeros(n, dtype=np.int32)
    flag_details = {}
    
    # Test 1: Range (with depth-specific max)
    cfg = config.get("range_test", {})
    vmin = cfg.get("min", 0.0)
    vmax = cfg.get("max", 5000.0)
    
    # Use depth-specific max if available
    depth_max = cfg.get("depth_max", {})
    if depth_cm in depth_max:
        vmax = depth_max[depth_cm]
    
    mask = test_range(values, vmin, vmax)
    qc_flags[mask] |= (1 << FLAG_RANGE)
    flag_details["range"] = mask
    
    # Test 2: Maintenance
    mask = test_maintenance(timestamps, maintenance_df)
    qc_flags[mask] |= (1 << FLAG_MAINTENANCE)
    flag_details["maintenance"] = mask
    
    # Backward-compatible qual
    qual = (qc_flags > 0).astype(np.int32)
    
    return qc_flags, qual, flag_details


# ── Temperature-specific QC ──────────────────────────────────────────────────

def apply_temperature_tests(values: ArrayLike, timestamps: ArrayLike, depth_cm: int,
                            maintenance_df: Optional[pd.DataFrame], config: dict) -> QCResult:
    """
    Run Temperature QC tests: range, maintenance, spike, gradient, flat-line.
    
    Parameters
    ----------
    values : np.ndarray
        Temperature values in °C.
    timestamps : np.ndarray
        Unix timestamps.
    depth_cm : int
        Sensor depth in cm.
    maintenance_df : pd.DataFrame or None
        Maintenance windows (debut, end columns).
    config : dict
        QC configuration from qc_config.yaml['Temperature'].
    
    Returns
    -------
    qc_flags : np.ndarray (int)
        Bitfield integer array.
    qual : np.ndarray (int)
        Simple 0/1 array for backward compatibility.
    flag_details : dict
        Per-test boolean masks.
    """
    n = len(values)
    qc_flags = np.zeros(n, dtype=np.int32)
    flag_details = {}
    
    # Test 1: Range (with depth-specific limits)
    cfg = config.get("range_test", {})
    vmin = cfg.get("min", 3.5)
    vmax = cfg.get("max", 30.0)
    
    # Use depth-specific range if available
    depth_ranges = cfg.get("depth_ranges", {})
    if depth_cm in depth_ranges:
        vmin, vmax = depth_ranges[depth_cm]
    
    mask = test_range(values, vmin, vmax)
    qc_flags[mask] |= (1 << FLAG_RANGE)
    flag_details["range"] = mask
    
    # Test 2: Maintenance
    mask = test_maintenance(timestamps, maintenance_df)
    qc_flags[mask] |= (1 << FLAG_MAINTENANCE)
    flag_details["maintenance"] = mask
    
    # Test 3: Spike
    cfg = config.get("spike_test", {})
    mask = test_spike(values, cfg.get("threshold", 2.0), cfg.get("window", 1))
    qc_flags[mask] |= (1 << FLAG_SPIKE)
    flag_details["spike"] = mask
    
    # Test 4: Gradient
    cfg = config.get("gradient_test", {})
    mask = test_gradient(values, timestamps, cfg.get("max_rate", 2.0))
    qc_flags[mask] |= (1 << FLAG_GRADIENT)
    flag_details["gradient"] = mask
    
    # Test 5: Flat-line
    cfg = config.get("flatline_test", {})
    mask = test_flatline(values, cfg.get("tolerance", 0.001), cfg.get("min_count", 30))
    qc_flags[mask] |= (1 << FLAG_FLATLINE)
    flag_details["flatline"] = mask
    
    # Backward-compatible qual
    qual = (qc_flags > 0).astype(np.int32)
    
    return qc_flags, qual, flag_details


#### Calibration correction of the DO data
#
# Applies time- and depth-specific corrections to O2 and O2Sat measurements.

class CalibrationCorrector:
    """
    Manages calibration data and applies corrections to DO measurements

    This class loads sensor calibration data from CSV and applies time-
    and depth-specific corrections to dissolved oxygen measurements.

    Attributes:
        calibration_csv_path: Path to calibration CSV file
        calibration_data: DataFrame containing processed calibration records
    """

    def __init__(self, calibration_file_path):
        """
        Initialize calibration corrector

        Args:
            calibration_file_path: Path to calibration file (.xlsx or .csv)
        """
        self.calibration_file_path = calibration_file_path
        self.calibration_data = self.load_calibration_data()

    def load_calibration_data(self):
        """
        Load and preprocess calibration data from CSV

        Returns:
            DataFrame with processed calibration data

        Raises:
            FileNotFoundError: If calibration CSV file doesn't exist
            ValueError: If CSV is missing required columns
        """
        # Check if file exists
        if not Path(self.calibration_file_path).exists():
            raise FileNotFoundError(f"Calibration file not found: {self.calibration_file_path}")

        # Load file depending on extension
        file_path_str = str(self.calibration_file_path)
        if file_path_str.endswith('.xlsx') or file_path_str.endswith('.xls'):
            df_raw = pd.read_excel(self.calibration_file_path, header=None)
            # Find the row containing 'Index' in the first column
            header_idx = df_raw[df_raw.iloc[:, 0] == 'Index'].index[0]
            df = df_raw.iloc[header_idx+1:].reset_index(drop=True)
            df.columns = df_raw.iloc[header_idx]
        else:
            df = pd.read_csv(self.calibration_file_path, sep=',', encoding='utf-8-sig')

        # Verify required columns exist
        required_columns = ['Index', 'Beginning deployment', 'End Deployment',
                          'Correction_percent', 'Depth_m']
        missing_columns = [col for col in required_columns if col not in df.columns]
        if missing_columns:
            raise ValueError(f"CSV missing required columns: {missing_columns}")

        # Remove rows with missing Index
        df = df.dropna(subset=['Index'])

        # Convert depth from meters to centimeters
        df["Depth_cm"] = df["Depth_m"].astype(float) * 100

        # Convert date columns to datetime
        df['Beginning_deployment'] = pd.to_datetime(df['Beginning deployment'], errors='coerce')
        df['End_Deployment'] = pd.to_datetime(df['End Deployment'], errors='coerce')

        # Remove rows with missing critical data
        df = df.dropna(subset=['Beginning_deployment', 'End_Deployment', 'Correction_percent', 'Depth_cm'])

        return df

    @staticmethod
    def calculate_o2sat_garcia_gordon_1992(temperature_C, salinity=0):
        """
        Calculate oxygen saturation concentration using Garcia & Gordon (1992)

        This is the most precise method for calculating dissolved oxygen saturation
        concentrations, accounting for both temperature and salinity effects.

        Args:
            temperature_C: Water temperature in degrees Celsius (scalar or array)
            salinity: Water salinity in psu (default=0 for freshwater)

        Returns:
            O2 saturation concentration in mg/L

        Reference:
            Garcia, H. E., & Gordon, L. I. (1992). Oxygen solubility in seawater:
            Better fitting equations. Limnology and Oceanography, 37(6), 1307-1312.

        Notes:
            - For Lake Geneva freshwater, use salinity=0
            - Formula calculates µmol/kg then converts to mg/L
            - Assumes water density: 1.0 kg/L (freshwater) or 1.025 kg/L (seawater)
        """
        # Coefficients for µmol/kg (from Garcia & Gordon 1992, Table 1)
        A0 = 5.80871
        A1 = 3.20291
        A2 = 4.17887
        A3 = 5.10006
        A4 = -9.86643e-2
        A5 = 3.80369

        B0 = -7.01577e-3
        B1 = -7.70028e-3
        B2 = -1.13864e-2
        B3 = -9.51519e-3

        C0 = -2.75915e-7

        # Calculate scaled temperature
        # Ts = ln[(298.15 - T) / (273.15 + T)]
        Ts = np.log((298.15 - temperature_C) / (273.15 + temperature_C))

        # Calculate ln(C_O2*) - Garcia & Gordon equation
        ln_C_O2_star = (A0 + A1*Ts + A2*Ts**2 + A3*Ts**3 + A4*Ts**4 + A5*Ts**5 +
                        salinity * (B0 + B1*Ts + B2*Ts**2 + B3*Ts**3) +
                        C0 * salinity**2)

        # Convert to µmol/kg
        C_O2_star_umol_kg = np.exp(ln_C_O2_star)

        # Convert from µmol/kg to mg/L
        # Molecular weight of O2 = 31.999 g/mol
        # For freshwater: density ≈ 1.0 kg/L
        # For seawater: density ≈ 1.025 kg/L
        density_kg_L = 1.0 if salinity == 0 else 1.025

        O2sat_mg_L = C_O2_star_umol_kg * (31.999 / 1000) * density_kg_L

        return O2sat_mg_L

    def get_correction_records(self, depth_cm, time_start, time_end, sensor_id=None):
        """
        Get applicable calibration correction records

        Filters calibration data by depth and time range to find applicable
        correction records for a given sensor deployment.

        Args:
            depth_cm: Sensor depth in centimeters
            time_start: Start of data time range (datetime or timestamp)
            time_end: End of data time range (datetime or timestamp)
            sensor_id: Optional sensor ID for filtering

        Returns:
            DataFrame with matching calibration records (may be empty)

        Notes:
            - A calibration period overlaps if:
              calibration_end >= data_start AND calibration_start <= data_end
            - Depth must match exactly (no tolerance)
        """
        # Convert timestamps to datetime if needed
        if isinstance(time_start, (int, float)):
            time_start = pd.to_datetime(time_start, unit='s')
        if isinstance(time_end, (int, float)):
            time_end = pd.to_datetime(time_end, unit='s')

        # Filter by depth (exact match)
        df_filtered = self.calibration_data[self.calibration_data['Depth_cm'] == depth_cm].copy()

        # Filter by overlapping time range
        # A calibration period overlaps if:
        # - Calibration end date >= data start date AND
        # - Calibration beginning date <= data end date
        df_filtered = df_filtered[
            (df_filtered['End_Deployment'] >= time_start) &
            (df_filtered['Beginning_deployment'] <= time_end)
        ]

        # Optionally filter by sensor ID
        if sensor_id is not None:
            df_filtered = df_filtered[df_filtered['Sensor'] == sensor_id]

        return df_filtered

    def apply_correction(self, do_obj):
        """
        Apply calibration correction to DO_data object

        This method applies time- and depth-specific calibration corrections
        to dissolved oxygen measurements. Corrections are applied in-place
        to the DO_data object.

        Modifies in-place:
            - do_obj.O2: Corrected oxygen concentration (mg/L)
            - do_obj.O2Sat: Corrected oxygen saturation (%)

        Args:
            do_obj: DO_data instance with O2, O2Sat, temperature, depth, and time data

        Returns:
            dict: Summary of applied corrections with keys:
                - depth_cm: Sensor depth
                - n_corrections: Number of calibration periods applied
                - corrections: List of correction details
            Returns None if no corrections were applied

        Raises:
            ValueError: If required data (temperature, O2, O2Sat) is missing
        """
        # Get depth from DO_data object
        depth_cm = do_obj.depth if hasattr(do_obj, 'depth') else None

        if depth_cm is None:
            return None

        # Verify required data exists
        if not hasattr(do_obj, 'O2') or do_obj.O2 is None:
            raise ValueError("DO_data object missing O2 data")
        if not hasattr(do_obj, 'O2Sat') or do_obj.O2Sat is None:
            raise ValueError("DO_data object missing O2Sat data")
        if not hasattr(do_obj, 'temperature') or do_obj.temperature is None:
            raise ValueError("DO_data object missing temperature data (required for Garcia & Gordon 1992)")
        if not hasattr(do_obj, 'time') or do_obj.time is None:
            raise ValueError("DO_data object missing time data")

        # Get time range
        time_start = pd.to_datetime(do_obj.time.min(), unit='s')
        time_end = pd.to_datetime(do_obj.time.max(), unit='s')

        # Find matching calibration records
        calibration_records = self.get_correction_records(depth_cm, time_start, time_end)

        if calibration_records.empty:
            # No calibration data available for this depth/time
            return None

        # Apply corrections for each calibration period
        corrections_applied = []

        for idx, cal_record in calibration_records.iterrows():
            # Find data points within this calibration period
            # Create time series with matching index for proper boolean masking
            time_series = pd.Series(
                pd.to_datetime(do_obj.time.values, unit='s'),
                index=do_obj.time.index
            )
            time_mask = (
                (time_series >= cal_record['Beginning_deployment']) &
                (time_series <= cal_record['End_Deployment'])
            )

            n_points = time_mask.sum()
            if n_points == 0:
                continue

            # Get correction factor (in percentage points)
            correction_percent = cal_record['Correction_percent']

            # ========================================
            # CORRECTION FORMULA (MODIFIABLE)
            # ========================================
            # Formula: Subtract sensor error from O2Sat, then recalculate O2 using Garcia & Gordon (1992)
            #
            # The Correction_percent represents the sensor ERROR (how much the sensor is off by).
            # We SUBTRACT the error to get the true value.
            #
            # Step 1: Correct O2Sat by subtracting sensor error (in percentage points)
            # Example: O2Sat = 85%, correction_percent = -2.70% (sensor reads 2.70% too LOW)
            #          O2Sat_corrected = 85 - (-2.70) = 87.70% (corrected UP)
            # Example: O2Sat = 85%, correction_percent = +1.54% (sensor reads 1.54% too HIGH)
            #          O2Sat_corrected = 85 - (1.54) = 83.46% (corrected DOWN)
            #
            # Step 2: Calculate O2sat (saturation concentration) using Garcia & Gordon (1992)
            #
            # Step 3: Recalculate O2 concentration from corrected saturation
            #          O2_corrected = (O2Sat_corrected / 100) * O2sat
            #
            # ⚠️ MODIFY THIS SECTION IF FORMULA CHANGES
            # ========================================

            # Get integer positions where time_mask is True to avoid index alignment issues
            mask_positions = np.where(time_mask)[0]

            # Extract data as numpy arrays using integer positions
            O2Sat_values = do_obj.O2Sat.iloc[mask_positions].values
            temperature_values = do_obj.temperature.iloc[mask_positions].values

            # Step 1: Correct O2Sat by subtracting sensor error (percentage points)
            O2Sat_corrected = O2Sat_values - correction_percent

            # Step 2: Calculate O2 saturation concentration using Garcia & Gordon (1992)
            # For Lake Geneva freshwater (salinity = 0)
            O2sat_mg_L = self.calculate_o2sat_garcia_gordon_1992(
                temperature_C=temperature_values,
                salinity=0  # ⚠️ MODIFIABLE: Freshwater assumption for Lake Geneva
            )

            # Step 3: Recalculate O2 concentration from corrected saturation percentage
            O2_corrected = (O2Sat_corrected / 100.0) * O2sat_mg_L

            # Apply corrections to data object using integer positions (avoids index alignment issues)
            do_obj.O2Sat.iloc[mask_positions] = O2Sat_corrected
            do_obj.O2.iloc[mask_positions] = O2_corrected

            # ========================================
            # END CORRECTION FORMULA
            # ========================================

            # Record applied correction
            corrections_applied.append({
                'sensor': cal_record['Sensor'],
                'correction_percent': correction_percent,
                'n_points': n_points,
                'period_start': cal_record['Beginning_deployment'].date(),
                'period_end': cal_record['End_Deployment'].date()
            })

        if corrections_applied:
            return {
                'depth_cm': depth_cm,
                'n_corrections': len(corrections_applied),
                'corrections': corrections_applied
            }

        return None


#### Define Classes
class ppmooring(object):
    ## Parent class for all the different sensors of the ppmooring module
    _config = None
    
    @classmethod
    def get_config(cls):
        if cls._config is None:
            cls._config = get_config()
        return cls._config

    def set_output(self, output_dir, level, instrument):
        ## Set output name and create output folder it not already existing
       
        # Set output folder and output file name
        if level == "L2":    
            outfolder = output_dir + self.Level2_outfolder
            date_time = datetime.utcfromtimestamp(self.time_grid.iloc[0]).strftime('%Y%m%d_%H%M%S')
            self.outfile = outfolder + "/" + level + "_" + instrument + "_" + date_time + ".nc"
        elif level == "L1B":
            outfolder = output_dir + self.outfolders[self.folder]

            preoutfile = self.name1B.replace("temp/","")
            preoutfile = preoutfile.replace("_temp", "")
            preoutfile = preoutfile.replace("Level1A", "Level1B")
            preoutfile = preoutfile.replace("L1A", "L1B")
            preoutfile = preoutfile.replace("1A", "1B")
            # Strip extension (.csv or .nc)
            preoutfile = os.path.splitext(preoutfile)[0]

            self.outfile = preoutfile + '.nc'
            self.outfilecsv = preoutfile + '.csv'
        else:
            outfolder = output_dir + self.outfolders[self.folder]
            date_time = datetime.utcfromtimestamp(self.time.iloc[0]).strftime('%Y%m%d_%H%M%S')
            self.outfile = outfolder + "/" + level + "_" + instrument + "_" + \
                str(self.depth) + "cm" + "_" + date_time + ".nc"
            self.outfilecsv = outfolder + "/" + level + "_" + instrument + "_" + \
                              str(self.depth) + "cm" + "_" + date_time + ".csv"
        
        # Make sure the output folder exists
        if not os.path.exists(outfolder):
            os.makedirs(outfolder)
            
    def to_NetCDF(self, *args):
        ## Save data to NetCDF file
        
        # Check which variables to store in NetCDF file
        if args:
            self.var_dict = {}
        else:
            logger.warning("No variables will be stored to NetCDF file")
            
        # Store variables and attributes
        for level in args:
            self.var_dict.update(getattr(self, level+"_dict"))
                
        # Sort dict alphabetically
        sorted_dict = {}
        for key, values in sorted(self.var_dict.items()): 
            sorted_dict[key] = values 
        self.var_dict = sorted_dict
            
        # Make sure NetCDF files are closed
        try: ncfile.close() 
        except: pass
                    
        # Check if file already exists
        if os.path.isfile(self.outfile):
            # Open file in 'append' mode
            ncfile = netCDF4.Dataset(self.outfile,mode='a',format='NETCDF4_CLASSIC')
            
            # Add Level 2 data
            if "L2" in args:
                # Add data to variables
                for key, values in self.var_dict.items():
                    self.addDataNetCDF(getattr(self,key),0,ncfile,"L2",**values)
                    
                # Close NetCDF file
                ncfile.close()
            
            # Add Level 1A/1B data
            else:
                # Check if last timestamp has already been added, if not data will be added
                if self.time.iloc[-1] in ncfile.variables['time'][:]:
                    logger.info("Duplicated run, no data added")
                    ncfile.close() 
                    return
                
                else:
                    # Get index of last saved timestamp
                    time = ncfile.variables['time']
                    idx = len(time[:])
                    
                    # Add data to coordinate variable 'time'
                    time[:] = np.append(time[:], self.time)
                    
                    # Add data to variables
                    for key, values in self.var_dict.items():
                        self.addDataNetCDF(getattr(self,key),idx,ncfile,**values)
                        
                    # Close NetCDF file
                    ncfile.close()
                
        else:           
            # Create new NetCDF file
            ncfile = netCDF4.Dataset(self.outfile, mode='w',format='NETCDF4_CLASSIC') 
            
            # Add global attributes
            for key in self.att:
                setattr(ncfile, key, self.att[key])
    
            # Create Dimensions
            for key, values in self.dim_dict.items():
                ncfile.createDimension(values['dim_name'], values['dim_size'])
                
            # Add data to coordinate variables
            for key, values in self.dim_var_dict.items():
                self.addVarNetCDF(getattr(self,key),ncfile,**values)
                
            # Add data to variables
            if "L2" in args:
                # Create variable for each depth
                for var in [*self.L2_dict]:
                    temp_dict = deepcopy(self.L2_dict[var])
                    for entry in list(self.depthdict.values()):
                        # Include variable name in NetCDF variable name to avoid duplicates
                        # Use depth in cm (integer) to avoid decimal points in variable names
                        depth_cm = int(entry)
                        temp_dict["var_name"] = var + "_" + str(depth_cm) + "cm"

                        var_temp = ncfile.createVariable(temp_dict["var_name"], np.float64, temp_dict["dim"])
                        var_temp.units = temp_dict["unit"]
                        var_temp.long_name = temp_dict["longname"]
                
                # Add data of current depth
                for key, values in self.var_dict.items():
                    self.addDataNetCDF(getattr(self,key),0,ncfile,"L2",**values)    
            else:
                for key, values in self.var_dict.items():
                    self.addVarNetCDF(getattr(self,key),ncfile,**values)
                
            # Close NetCDF file
            ncfile.close()
    
    def addDataNetCDF(self,data,idx,ncfile,*args,**kwargs):
        ## Add data to existing NetCDF file (only along time dimension)
         
        # Check if variable exists, if not create new one
        try: 
            var = ncfile.variables[kwargs["var_name"]]
        except:
            var = ncfile.createVariable(kwargs["var_name"], np.float64, kwargs["dim"])
            var.units = kwargs["unit"]
            var.long_name = kwargs["longname"]
        
        # Add Level 2 data if args specified, otherwise Level 1A/1B
        if args:   
            # Add data to variable
            var[self.idx_grid] = data
            
        else:
            # Check dimension size
            if "dim" in kwargs.keys():
                ndim = len(kwargs["dim"])
            else:
                ndim = 0        
                
            # Add data to variable
            if ndim == 2:
                var[:,idx:] = data
            elif ndim == 1:
                var[idx:] = data
            else:
                var[:] = data
            
    def addVarNetCDF(self,data,ncfile,**kwargs):
        ## Create variable in new NetCDF file
        
        # Create variable and check for dimension variable
        if "dim" in kwargs.keys():
            var = ncfile.createVariable(kwargs["var_name"], np.float64, kwargs["dim"])
            opt = 1
        else:
            var = ncfile.createVariable(kwargs["var_name"], np.float64)
            opt = 0
           
        # Add unit and longname
        var.units = kwargs["unit"]
        var.long_name = kwargs["longname"]
               
        # Add data to variable
        if opt ==1:
            var[:] = data # vector with dimension >1
        else:
            var[:] = np.ravel(data)[0]
        
    def get_files_to_read_Level1A_nc(self, input_dir):
        """Get list of Level 1A NetCDF files of this sensor, sorted by folder then name."""
        folders = os.listdir(input_dir)
        folders = fnmatch.filter(folders, self.level1A_folder_pattern)

        self.filestoread = []
        for folder in sorted(folders):
            folder_path = os.path.join(input_dir, folder)
            nc_files = fnmatch.filter(os.listdir(folder_path), '*.nc')
            nc_files = [os.path.join(folder_path, f).replace(os.sep, '/') for f in sorted(nc_files)]
            self.filestoread.extend(nc_files)

    def read_data_Level1A_nc(self, infile):
        """Read a Level 1A NetCDF file: time, depth and the variables of this sensor."""
        ds = netCDF4.Dataset(infile)

        # Extract variables
        time_arr = ds['time'][:].data
        values = {name: ds[name][:].data for name in self.level1A_variables}
        depth = float(ds['depth'][:].data) if ds['depth'].shape else float(ds['depth'][...])
        ds.close()

        # Sort by time
        sort_idx = np.argsort(time_arr)

        # Store as pandas Series (backward compatible with existing code)
        self.time = pd.Series(time_arr[sort_idx])
        for name, data in values.items():
            setattr(self, name, pd.Series(data[sort_idx]))
        self.depth = depth

        # Set folder name for output path resolution (key of self.outfolders)
        self.folder = f"{self.level1A_folder_prefix}{int(depth)}cm"
        self.name1B = infile.replace(os.sep, '/')

    def mask_data(self):
        ## Replace bad values from Level 1B by NaN
        
        for var in [*self.L2_dict]:
            idx = getattr(self, var+"_qual") > 0
            
            getattr(self, var)[idx] = np.nan
            
    def grid_data(self, start_period, end_period):
        ## Store data in grid before outputting it as Level 2
        
        # Define time grid: one point every 5 minutes from start_period, end_period excluded
        step = dt.timedelta(minutes = 5)
        grid = [start_period]
        while grid[-1] + step < end_period:
            grid.append(grid[-1] + step)
        time_grid = np.array(grid)

        # Timestamps of the measurements, as integer microseconds (exact arithmetic)
        date_us = datetimes_to_microseconds([datetime.utcfromtimestamp(time) for time in self.time])

        # Consider duplicate time values because of change from winter time to summer time
        date_us, idx_unique = np.unique(date_us, return_index=True)
        for var in [*self.L2_dict]:
            setattr(self, var, getattr(self, var).iloc[idx_unique])
            setattr(self, var, getattr(self, var).reset_index(drop = True))

        # Keep the measurements that are the closest ones to a time_grid value
        self.idx_grid, idx_data = match_to_grid(date_us, datetimes_to_microseconds(grid))

        # Update time and measurement data
        self.time_grid = pd.Series([value.replace(tzinfo=timezone.utc).timestamp() for value in time_grid])
        for var in [*self.L2_dict]:
            setattr(self, var, getattr(self, var)[idx_data])
            setattr(self, var, getattr(self, var).reset_index(drop = True))
        
    def Level2_to_NetCDF(self, directory, instrument):
        ## Save Level2 data to NetCDF files where the data is seperated
        ## s.t. each file spans over a time period of 10 days
        
        # Define time interval in days 
        interval = 10
        
        # Get the datetime objects corresponding to the timestamps of the measurements
        date_time = np.array([datetime.utcfromtimestamp(time) for time in self.time])
        
        # Set beginning of time grid
        grid_begin = datetime(2018, 1, 1)
                    
        # Set end of time grid
        grid_end = date_time[-1]

        # Periods of 10 days from grid_begin. Period p covers
        # [grid_begin + p*step, grid_begin + (p+1)*step); the last period is the one
        # that contains grid_end.
        step = dt.timedelta(days = interval)
        n_periods = max(0, (grid_end - grid_begin) // step + 1)

        # Period of each measurement (integer microseconds, exact arithmetic)
        date_us = datetimes_to_microseconds(date_time)
        begin_us = datetimes_to_microseconds([grid_begin])[0]
        step_us = step // dt.timedelta(microseconds = 1)
        period_of = (date_us - begin_us) // step_us

        # Loop over the periods that contain measurements, in chronological order
        for period in np.unique(period_of):
            period = int(period)

            # Slice values to store
            idx = np.where(period_of == period)[0]
            period_start = grid_begin + period * step
            period_end = period_start + step

            # Measurements outside the grid are not stored
            if period < 0 or period >= n_periods:
                continue
            else:
                # Create copy of current class instance to save current period to NetCDF
                period_obj = deepcopy(self)
                
                # Slice data
                period_obj.time = period_obj.time[idx]
                for var in [*period_obj.L2_dict]:
                    setattr(period_obj, var, getattr(period_obj,var)[idx])
                    # Include variable name in NetCDF variable name to avoid duplicates
                    # Use depth in cm (integer) to avoid decimal points in variable names
                    depth_cm = int(np.ravel(period_obj.depth)[0])
                    period_obj.L2_dict[var]["var_name"] = var + "_" + str(depth_cm) + "cm"
                    
                # Grid data
                period_obj.grid_data(period_start, period_end)
                    
                # Change dictonary key to guarantee compability with the .to_NetCDF method
                period_obj.dim_dict["time_grid"] = period_obj.dim_dict.pop("time")
                period_obj.dim_var_dict["time_grid"] = period_obj.dim_var_dict.pop("time")
                
                # Set output name
                period_obj.set_output(directory, "L2", instrument)            
                
                # Save to NetCDF
                period_obj.to_NetCDF("L2")

    def to_csv(self,*args):
        ## Save data to csv files
        # Check which variables to store in csv file
        if args:
            self.var_dict = {}
        else:
            logger.warning("No variables will be stored to CSV file")

        # Store variables and attributes
        for level in args:
            self.var_dict.update(getattr(self, level + "_dict"))

        # Sort dict alphabetically
        sorted_dict = {}
        for key, values in sorted(self.var_dict.items()):
            sorted_dict[key] = values
        self.var_dict = sorted_dict

        # Check if file already exists
        if os.path.isfile(self.outfilecsv):
            logger.info("Duplicated run, no csv created")
        else:
            data = {'time' : getattr(self,'time')}
            for key,values in self.var_dict.items():
                data[key] = getattr(self,key)

            data = pd.DataFrame.from_dict(data)
            data.to_csv(self.outfilecsv,index=False)

class DO_data(ppmooring):
    # Level 1A files of this sensor: folder name pattern, variables stored, and the
    # prefix of the key used in self.outfolders
    level1A_folder_pattern = 'LexplorePPMooringDO*'
    level1A_variables = ("O2", "O2Sat", "temperature", "O2_qual", "O2Sat_qual", "temperature_qual")
    level1A_folder_prefix = "DO_"

    ## Class for the DO data

    def __init__(self):
        ## Set attributes for this class

        #root
        self.root = str(ppmooring.get_config().ppmooring_root)
        # Output folders
        self.outfolders = {
            "DO_50cm": "/LexplorePPMooringDO_50cm",
            "DO_250cm": "/LexplorePPMooringDO_250cm",
            "DO_500cm": "/LexplorePPMooringDO_500cm",
            "DO_1000cm": "/LexplorePPMooringDO_1000cm",
            "DO_1500cm": "/LexplorePPMooringDO_1500cm",
            "DO_2000cm": "/LexplorePPMooringDO_2000cm",
            "DO_3000cm": "/LexplorePPMooringDO_3000cm",
            "DO_5000cm": "/LexplorePPMooringDO_5000cm",
            "DO_10000cm": "/LexplorePPMooringDO_10000cm"
        }

        self.Level2_outfolder = "/LexplorePPMooringDO"
        
        # Matrix depth - folder
        self.depthdict = {
            "DO_50cm": 50,
            "DO_250cm": 250,
            "DO_500cm": 500,
            "DO_1000cm": 1000,
            "DO_1500cm": 1500,
            "DO_2000cm": 2000,
            "DO_3000cm": 3000,
            "DO_5000cm": 5000,
            "DO_10000cm": 10000
        }
        
        # General attributes
        self.att = {
            "institution": "EPFL",
            "source": "Mooring DO",
            "references": "LéXPLORE commun instruments Viet Tran Khac <viet.tran-khac@inra.fr>",
            "history": "See history on Renku",
            "conventions": "CF 1.7",
            "comment": "Dissolved Oxygen from Mooring M2 on Lexplore Platform in Lake Geneva",
            "title": "Lexplore Mooring DO"
        }
        
        # Dimension attributes
        self.dim_dict = {
            'time': {'dim_name':'time', 'dim_size':None}
        }
        
        # Coordinate variable attributes
        self.dim_var_dict = {
            'time': {'var_name':'time', 'dim':('time',), 'unit':'seconds since 1970-01-01 00:00:00', 'longname':'time'}
        }
        
        # Level 1A - Variable attributes
        self.L1A_dict = {
            'depth': {'var_name':'depth','unit':'cm', 'longname':'depth'},
            'O2Sat': {'var_name':'O2Sat', 'dim':('time',), 'unit':'%', 'longname':'oxygen saturation'},
            'O2Sat_qual': {'var_name': 'O2Sat_qual', 'dim': ('time',),
                           'unit': '0 = nothing to report, 1 = more investigation',
                           'longname': 'oxygen saturation quality flag'},
            'O2': {'var_name':'O2', 'dim':('time',), 'unit':'mgO2/L', 'longname':'dissolved_oxygen'},
            'O2_qual': {'var_name': 'O2_qual', 'dim': ('time',),
                        'unit': '0 = nothing to report, 1 = more investigation',
                        'longname': 'dissolved_oxygen quality flag'},
            'temperature': {'var_name':'temperature', 'dim':('time',), 'unit':'degree_Celsius', 'longname':'water temperature'},
            'temperature_qual': {'var_name': 'temperature_qual', 'dim': ('time',),
                                'unit': '0 = nothing to report, 1 = more investigation',
                                'longname': 'water temperature quality flag'}
        }

        # Level 1B - Quality flag attributes
        self.L1B_dict = {
            'O2Sat_qual': {'var_name':'O2Sat_qual', 'dim':('time',), 'unit':'0 = nothing to report, 1 = more investigation', 'longname':'oxygen saturation quality flag'},
            'O2_qual': {'var_name':'O2_qual', 'dim':('time',), 'unit':'0 = nothing to report, 1 = more investigation', 'longname':'dissolved_oxygen quality flag'},
            'temperature_qual': {'var_name':'temperature_qual', 'dim':('time',), 'unit':'0 = nothing to report, 1 = more investigation', 'longname':'water temperature quality flag'}
        }

        # Level 2 - Variable attributes
        self.L2_dict = {
            'O2': {'var_name':'O2', 'dim':('time',), 'unit':'mgO2/L', 'longname':'dissolved_oxygen'},
            'O2Sat': {'var_name':'O2Sat', 'dim':('time',), 'unit':'%', 'longname':'oxygen saturation'},
            'temperature': {'var_name':'temperature', 'dim':('time',), 'unit':'degree_Celsius', 'longname':'water temperature'}
        }

    def read_data_rbr(self, infile):
        ## Read RBR data
        
        # Store infile variable
        self.infile = infile
        
        # Get input folder
        self.folder = self.infile.rsplit("/")[-3]
                
        # Import metadata and get indicies of metadata
        with open(self.infile, encoding="utf8", errors='ignore') as f:
            for idx, line in enumerate(f):
                if 'Serial' in line:
                    idx_serial = idx
                if 'Channel[1].calibration' in line:
                    idx_calibr = idx
                if 'NumberOfSamples' in line:
                    idx_meta = idx + 1
                    break
            f.seek(0)
            lines = f.readlines()
            line_name = lines[idx_meta + 1]
            meta_data = lines[0:idx_meta]
        
        # Extract meta-data
        numsensor = meta_data[idx_serial].rsplit("=")[1].rstrip("\n")
        numcalibration = meta_data[idx_calibr].rsplit("=")[1].rstrip("\n")
    
        # Read data
        #20221201-insert encoding = 'unicode_escape', because utf-8 encoding bugs

        df = pd.read_csv(self.infile, sep=r'\s+', header=None, skiprows=idx_meta + 2, encoding='unicode_escape')
        
        # Take into account that files do not always have the same columns
        if "doxy07" in line_name or "doxy01" in line_name:
            df = df.drop([3, 5, 6], axis=1)
            df.columns = ["Date", "Time", "Temp", "O2Sat", "O2"]
            
            cols = ['O2Sat', 'O2']
            df[cols] = df[cols].apply(pd.to_numeric, errors='coerce', axis=1)
        else:
            df = df.drop([4], axis = 1)
            df.columns = ["Date", "Time", "O2", "Temp", "O2Sat"]
            
            cols = ['O2Sat', 'O2']
            df[cols] = df[cols].apply(pd.to_numeric, errors='coerce', axis=1)
            
            # Convert O1 µmol/L to mg/L&4
            df["O2"] = df["O2"] * 32 / 1000
            
        df = df.reset_index(drop=True)
    
        # Replace date and time by timestamp
        datehour = df["Date"] + " " + df["Time"]
        datehour = datehour.apply(lambda x: datetime.strptime(x, '%Y-%m-%d %H:%M:%S.%f'))
        df["Timestamp"] = [datetime.timestamp(s) for s in datehour]
    
        df = df[["Timestamp", "O2Sat", "O2", "Temp"]]

        # Sort by date and reset index
        df = df.sort_values(by=['Timestamp'])
        self.df = df.reset_index(inplace=False, drop=True)

        # Sensor data, depth and point
        self.meta_data = {
            "sensor": numsensor,
            "calibration": numcalibration,
            "instrument": "rbr",
            "depth": self.depthdict[self.folder],
            "point": self.folder
        }

        # Store variables
        self.time = df["Timestamp"] # need to be checked
        self.depth = self.meta_data["depth"]
        self.O2Sat = df["O2Sat"]
        self.O2 = df["O2"]
        self.temperature = df["Temp"]

        # Create vector quality for data quality flags at Level 1A
        self.O2_qual = np.zeros(len(self.O2))
        self.O2Sat_qual = np.zeros(len(self.O2_qual))
        self.temperature_qual = np.zeros(len(self.temperature))

    def read_data_rbr_excel(self, infile):
        ## Read RBR data
        # Store infile variable
        self.infile = infile

        # Get input folder
        self.folder = self.infile.rsplit("/")[-3]

        # Import metadata and get indicies of metadata
        metadata = pd.read_excel(open(infile, 'rb'), sheet_name='Metadata')

        idx = np.where(metadata.iloc[:, 0].str.contains("Serial") == True)[0][0] + 1
        numsensor = metadata.iloc[idx, 0]

        idx = np.where(metadata.iloc[:, 0].str.contains("Atmospheric pressure") == True)[0][0] + 1
        atmospheric_pressure = metadata.iloc[idx, 0]
        density = metadata.iloc[idx, 1]
        meta_Temperature = metadata.iloc[idx, 3]
        pressure = metadata.iloc[idx, 4]
        meta_conductivity = metadata.iloc[idx, 5]
        salinity = metadata.iloc[idx, 6]
        conductivity_specific = metadata.iloc[idx, 7]

        # Read data
        df = pd.read_excel(open(infile, 'rb'), sheet_name='Data', skiprows=1)

        ncol = df.shape[1]
        # Take into account that files do not always have the same columns
        if ncol ==7:
            #Data rbr with Temperature	Pressure	O₂ saturation	Sea pressure	Depth	 O₂ concentration (umol/L)
            df.columns = ["Date","Temp","Pressure","O2Sat","Sea_pressure","Depth","O2"]

            ##df["O2"] =(df["O2"] +30.06) * 0.032 #Convert to mg/L
            termA = np.exp(7.7117 - 1.31403 * np.log((df["Temp"] + 45.93)))
            termB = np.exp(11.8571-(3840.7 / (df["Temp"]+273.15))-(216961/np.power((df["Temp"] + 273.15),2)))
            termC = (0.000975 - (0.00001426 * df["Temp"]) + (0.00000006436 * np.power(df["Temp"],2)))
            termABC = (termA * (df["Pressure"] / 10) * (1 - termB / df["Pressure"] / 10) * (1 - termC * df["Pressure"]))
            cpO2 = (termABC / (1 - np.exp(11.8571 - (3840.7 / (df["Temp"] + 273.15)) - (216961 / np.power((df["Temp"] + 273.15),2)))) / (1 - (0.000975 - (0.00001426 * df["Temp"]) + (0.00000006436 * np.power(df["Temp"],2)))))

            df["O2"]= df["O2Sat"]/100 * cpO2

            #remove all rows < or > than depth
            vdepth = self.depthdict[self.folder]
            idx = np.where((df["Depth"]<vdepth/100*0.8) | (df["Depth"]>vdepth/100*1.2))[0]
            df.loc[idx, "O2"] = 0
            df.loc[idx, "O2Sat"] = 0


            df = df.reset_index(drop=True)
            df["Timestamp"] = [datetime.timestamp(s) for s in df['Date']]

        if ncol == 5:
            df.columns = ["Date","O2","Temp","Corrected_phase","O2Sat"]
            df["O2"] = (df["O2"]+df["Corrected_phase"])*0.032
            ##df["O2Sat"] = df["O2Sat"]*1.0502 #1/p correction 400m

            df = df.reset_index(drop=True)

            # Replace date and time by timestamp
            df["Timestamp"] = [datetime.timestamp(s) for s in df['Date']]

        df = df[["Timestamp", "O2Sat", "O2", "Temp"]]

        # Sort by date and reset index
        df = df.sort_values(by=['Timestamp'])
        self.df = df.reset_index(inplace=False, drop=True)

        # Sensor data, depth and point
        self.meta_data = {
            "sensor": numsensor,
            "atmospheric_pressure": atmospheric_pressure,
            "density": density,
            "instrument": "rbr",
            "meta_conductivity": meta_conductivity,
            "depth": self.depthdict[self.folder],
            "point": self.folder
        }

        # Store variables
        self.time = df["Timestamp"]  # need to be checked
        self.depth = self.meta_data["depth"]
        self.O2Sat = df["O2Sat"]
        self.O2 = df["O2"]
        self.temperature = df["Temp"]

        # Create vector quality for data quality flags at Level 1A
        self.O2_qual = np.zeros(len(self.O2))
        self.O2Sat_qual = np.zeros(len(self.O2_qual))
        self.temperature_qual = np.zeros(len(self.temperature))

    def read_data_exo(self, infile):
        ## Read EXO data    
            
        # Store infile variable
        self.infile = infile
        
        # Get input folder
        self.folder = self.infile.rsplit("/")[-3]
            
        # Import metadata and get indicies of metadata
        with open(self.infile, encoding="utf8", errors='ignore') as f:
            for idx, line in enumerate(f):
                if 'Sonde ID' in line:
                    idx_sonde = idx
                    break
            f.seek(0)
            lines = f.readlines()
            meta_data = lines[0:26]

        # Extract meta-data
        numsensor = meta_data[idx_sonde].rsplit(";")[1].rsplit(" ")[1]    
        numcalibration = -99
        
        # Read data
        df = pd.read_csv(self.infile, sep=";", skiprows=26, header=None, encoding='unicode_escape')
        
        ncol = len(df.columns)
        # If data has 26 columns = no TSS, if 27 columns = TSS column included
        if ncol == 27:
            df = df.drop([2, 3, 4, 6], axis=1)
            df.columns = ["Date", "Time", "Batt", "pH", "pHmV", "NNO3mv", "NNO3",
                          "ChlaRFU", "Chla", "PCRFU", "PC", "Temp", "Cond", "SpCond", "Sal", "nLFCond",
                          "TDS", "Turbidity", "TSS", "O2Sat", "O2", "Pressure", "Depth"]
            
            # Add timestamp
            datehour = df["Date"] + " " + df["Time"]
            datehour = datehour.apply(lambda x: datetime.strptime(x, '%m/%d/%Y %H:%M:%S'))
            df["Timestamp"] = [datetime.timestamp(s) for s in datehour]
    
            df = df[["Timestamp", "Batt", "pH", "pHmV", "NNO3mv", "NNO3",
                     "ChlaRFU", "Chla", "PCRFU", "PC", "Temp", "Cond", "SpCond", "Sal", "nLFCond",
                     "TDS", "Turbidity", "TSS", "O2Sat", "O2", "Pressure", "Depth"]]
        else:
            df = df.drop([2, 3, 4, 6], axis=1)
            df.columns = ["Date", "Time", "Batt", "pH", "pHmV", "NNO3mv", "NNO3",
                          "ChlaRFU", "Chla", "PCRFU", "PC", "Temp", "Cond", "SpCond", "Sal", "nLFCond",
                          "TDS", "Turbidity", "O2Sat", "O2", "Pressure", "Depth"]
            df["TSS"] = -99
            
            # Add timestamp
            datehour = df["Date"] + " " + df["Time"]
            datehour = datehour.apply(lambda x: datetime.strptime(x, '%m/%d/%Y %H:%M:%S'))
            df["Timestamp"] = [datetime.timestamp(s) for s in datehour]
    
            df = df[["Timestamp", "Batt", "pH", "pHmV", "NNO3mv", "NNO3",
                     "ChlaRFU", "Chla", "PCRFU", "PC", "Temp", "Cond", "SpCond", "Sal", "nLFCond",
                     "TDS", "Turbidity", "TSS", "O2Sat", "O2", "Pressure", "Depth"]]
        
        df =df[["Timestamp", "O2Sat", "O2", "Temp"]]

        # Sort by date and reset index
        df = df.sort_values(by=['Timestamp'])
        self.df = df.reset_index(inplace=False, drop=True)

        # Sensor data, depth and point
        self.meta_data = {
            "sensor": numsensor,
            "calibration": numcalibration,
            "instrument": "exo",
            "depth": self.depthdict[self.folder],
            "point": self.folder
        }

        # Store variables
        self.time = df["Timestamp"] # need to be checked
        self.depth = self.meta_data["depth"]
        self.O2Sat = df["O2Sat"]
        self.O2 = df["O2"]
        self.temperature = df["Temp"]

        # Create vector quality for data quality flags at Level 1A
        self.O2_qual = np.zeros(len(self.O2))
        self.O2Sat_qual = np.zeros(len(self.O2_qual))
        self.temperature_qual = np.zeros(len(self.temperature))
        
    def read_data_minidot(self, infile):
        ## Read Minidot data       
        
        # Store infile variable
        self.infile = infile
        
        # Get input folder
        self.folder = self.infile.rsplit("/")[-3]
               
        # Import metadata and get indicies of metadata
        with open(self.infile, encoding="utf8", errors='ignore') as f:
            lines = f.readlines()
            meta_data = lines[0:3]
            
        # Extract meta-data
        numsensor = meta_data[0].rsplit("\n")[0]
        numcalibration = meta_data[1].rsplit(":")[2].rsplit('\n')[0]

        # Read data
        df = pd.read_csv(self.infile, sep=",", skiprows=3, header=None)
        df.columns= ["Timestamp", "Batt", "Temp", "O2", "Ratio"]
        
        # O2Sat column
        #Calculer O2Sat data
        DO_cp = 14.64 - (0.4227 * df["O2"]) + (0.009937 * np.power(df["O2"],2)) - (0.0001575 * np.power(df["O2"],3)) + (0.000001125 * np.power(df["O2"], 4))
        df["O2Sat"] = (df["O2"]/DO_cp)*100

        df =df[["Timestamp", "O2Sat", "O2", "Temp"]]

        # Sort by date and reset index
        df = df.sort_values(by=['Timestamp'])
        self.df = df.reset_index(inplace=False, drop=True)

        # Sensor data, depth and point
        self.meta_data = {
            "sensor": numsensor,
            "calibration": numcalibration,
            "instrument": "minidot",
            "depth": self.depthdict[self.folder],
            "point": self.folder
        }

        # Store variables
        self.time = df["Timestamp"] # need to be checked
        self.depth = self.meta_data["depth"]
        self.O2Sat = df["O2Sat"]
        self.O2 = df["O2"]
        self.temperature = df["Temp"]

        # Create vector quality for data quality flags at Level 1A
        self.O2_qual = np.zeros(len(self.O2))
        self.O2Sat_qual = np.zeros(len(self.O2_qual))
        self.temperature_qual = np.zeros(len(self.temperature))

    def quality_flags_v2(self):
        """
        Simplified QC with 6 independent, physically-motivated tests.

        Replaces both quality_flags_level0() and quality_flags_level1A()
        with a single call. Each test is independent — no cascading,
        no re-initialization between filters.

        Flag bits (stored in O2_qc integer array):
            bit 0 (1):  Range test         — hard physical limits
            bit 1 (2):  Maintenance mask   — known deployment windows
            bit 2 (4):  Climatological     — month+depth expected range
            bit 3 (8):  Spike test         — short deviations returning to baseline
            bit 4 (16): Gradient test      — rate of change exceeds physical max
            bit 5 (32): Flat-line test     — consecutive identical values

        O2_qual remains 0/1 for backward compatibility (1 if any bit is set).
        O2_qc is the full bitfield for diagnostics.
        """
        import yaml

        # Load QC configuration
        config_obj = get_config()
        qc_config_path = os.path.join(str(config_obj.scripts_path), 'qc_config.yaml')
        with open(qc_config_path, 'r') as f:
            qc_config = yaml.safe_load(f)
        do_config = qc_config.get('DO', {})

        # Prepare inputs
        values = np.copy(self.O2).astype(float)
        timestamps = np.asarray(self.time, dtype=float)
        depth_cm = int(np.ravel(self.depth)[0])

        # Load maintenance data
        try:
            maintenance_df = pd.read_csv(
                str(config_obj.maintenance_files['par_maintenance']), sep=";"
            )
        except Exception:
            maintenance_df = None
            
        # Try self.maintenance if available (set by main_ppmooring.py)
        if maintenance_df is None and hasattr(self, 'maintenance'):
            maintenance_df = self.maintenance

        # Load climatological border values
        try:
            border_val_df = pd.read_csv(
                str(config_obj.maintenance_files['border_values']), sep=";"
            )
        except Exception:
            border_val_df = None

        # Run all 6 tests
        qc_flags, qual, flag_details = apply_all_tests(
            values, timestamps, depth_cm,
            maintenance_df, border_val_df, do_config
        )

        # Store results
        self.O2_qc = qc_flags                        # bitfield (new)
        self.O2_qual = qual.astype(float)             # 0/1 backward compatible
        self.O2_flag_details = flag_details           # per-test masks (for diagnostics)

        # Propagate to O2Sat and temperature (same logic as original)
        self.O2Sat_qual = np.copy(self.O2_qual)
        self.temperature_qual = np.copy(self.O2_qual)

        # Print summary
        depth_str = f"{depth_cm}cm"
        logger.info("QC v2 - DO at %s:", depth_str)
        summarize_flags(qc_flags)

class PAR_data(ppmooring):
    # Level 1A files of this sensor: folder name pattern, variables stored, and the
    # prefix of the key used in self.outfolders
    level1A_folder_pattern = 'LexplorePPMooringPAR*'
    level1A_variables = ("PAR", "PAR_qual")
    level1A_folder_prefix = ""

    ## Class for the PAR data

    def __init__(self):
        ## Set attributes for this class
        self.root = str(ppmooring.get_config().ppmooring_root)
        # Output folders
        self.outfolders = {
            "RBR_PAR_50cm": "/LexplorePPMooringPAR_50cm",
            "RBR_PAR_250cm": "/LexplorePPMooringPAR_250cm",
            "RBR_PAR_500cm": "/LexplorePPMooringPAR_500cm",
            "RBR_PAR_1000cm": "/LexplorePPMooringPAR_1000cm",
            "RBR_PAR_2000cm": "/LexplorePPMooringPAR_2000cm",
            "RBR_PAR_3000cm": "/LexplorePPMooringPAR_3000cm",
            "50cm": "/LexplorePPMooringPAR_50cm",
            "250cm": "/LexplorePPMooringPAR_250cm",
            "500cm": "/LexplorePPMooringPAR_500cm",
            "1000cm": "/LexplorePPMooringPAR_1000cm",
            "2000cm": "/LexplorePPMooringPAR_2000cm",
            "3000cm": "/LexplorePPMooringPAR_3000cm"
        }
        self.Level2_outfolder = "/LexplorePPMooringPAR"

        # Matrix depth - folder
        self.depthdict = {
            "RBR_PAR_50cm":50,
            "RBR_PAR_250cm":250,
            "RBR_PAR_500cm":500,
            "RBR_PAR_1000cm":1000,
            "RBR_PAR_2000cm":2000,
            "RBR_PAR_3000cm":3000
        }

        # General attributes
        self.att = {
            "institution": "EPFL",
            "source": "Mooring DO",
            "references": "LéXPLORE commun instruments Viet Tran Khac <viet.tran-khac@inrae.fr>",
            "history": "See history on Renku",
            "conventions": "CF 1.7",
            "comment": "PAR (photosynthetically active radiation) from Mooring M1 on Lexplore Platform in Lake Geneva",
            "title": "Lexplore PPMooring PAR"
        }

        # Dimension attributes
        self.dim_dict = {
            'time': {'dim_name':'time', 'dim_size':None}
        }

        # Coordinate variable attributes
        self.dim_var_dict = {
            'time': {'var_name':'time', 'dim':('time',), 'unit':'seconds since 1970-01-01 00:00:00',
                     'longname':'time'}
        }

        # Level 1A - Variable attributes
        self.L1A_dict = {
            'depth': {'var_name':'depth', 'unit':'cm', 'longname':'depth'},
            'PAR': {'var_name':'PAR', 'dim':('time',), 'unit':'μmol/m²/s', 'longname':'PAR'},
            'PAR_qual': {'var_name': 'PAR_qual', 'dim': ('time',),
                         'unit': '0 = nothing to report, 1 = more investigation',
                         'longname': 'PAR quality flag'}

        }

        # Level 1B - Quality flag attributes
        self.L1B_dict = {
            'PAR_qual': {'var_name':'PAR_qual', 'dim':('time',), 'unit':'0 = nothing to report, 1 = more investigation',
                         'longname':'PAR quality flag'}
        }

        # Level 2 - Variable attributes
        self.L2_dict = {
            'PAR': {'var_name':'PAR', 'dim':('time',), 'unit':'μmol/m²/s', 'longname':'PAR'}
        }

    def read_data(self, infile):
        #Since update of software RBR version 6.13, data is converted into .xlsx, no longer .txt

        #Read txt files
        if (infile.endswith('txt')):
            ## Read PAR data

            # Store infile variable
            self.infile = infile

            # Get input folder
            self.folder = self.infile.rsplit("/")[-3]

            # Import metadata and get indicies of metadata
            with open(self.infile, encoding="utf8", errors='ignore') as f:
                for idx, line in enumerate(f):
                    if 'Serial' in line:
                        idx_serial = idx
                    if 'Channel[1].calibration' in line:
                        idx_calibr = idx
                    if 'NumberOfSamples' in line:
                        idx_meta = idx + 1
                        break
                f.seek(0)
                lines = f.readlines()
                meta_data = lines[0:idx_meta]

            # Extract meta-data
            numsensor = meta_data[idx_serial].rsplit("=")[1].rstrip("\n")
            numcalibration = meta_data[idx_calibr].rsplit("=")[1].rstrip("\n")

            # Read Data
            df = pd.read_csv(self.infile, sep=r'\s+', header=None, skiprows=idx_meta + 2,encoding="unicode_escape")
            df.columns = ["Date", "Time", "PAR"]

            # Replace date and time by timestamp
            datehour = df["Date"] + " " + df["Time"]
            datehour = datehour.apply(lambda x: datetime.strptime(x, '%Y-%m-%d %H:%M:%S.%f'))
            df["Timestamp"] = [datetime.timestamp(s) for s in datehour]

            df = df[["Timestamp", "PAR"]]

            # Sort by date and reset index
            df = df.sort_values(by=['Timestamp'])
            self.df = df.reset_index(inplace=False, drop=True)

            # Sensor data, depth and point
            self.meta_data = {
                "sensor": numsensor,
                "calibration": numcalibration,
                "instrument": "rbr",
                "depth": self.depthdict[self.folder],
                "point": self.folder
            }

            # Store variables
            self.time = df["Timestamp"]  # need to be checked
            self.depth = self.meta_data["depth"]
            self.PAR = df["PAR"]

            # Create vector quality for data quality flags at Level 1A
            self.PAR_qual = np.zeros(len(self.PAR))

        #Read xlsx files
        if (infile.endswith('xlsx')):
            ## Read PAR data

            # Store infile variable
            self.infile = infile

            # Get input folder
            self.folder = self.infile.rsplit("/")[-3]

            # Import metadata
            metadata = pd.read_excel(open(infile, 'rb'),sheet_name = 'Metadata')

            #Get metadata
            metadata_nrow = metadata.shape[0]
            metadata_ncol = metadata.shape[1]

            idx = np.where(metadata.iloc[:,0].str.contains("Serial")==True)[0][0] +1
            numsensor = metadata.iloc[idx,0]

            idx = np.where(metadata.iloc[:,0].str.contains("Atmospheric pressure")==True)[0][0] +1
            atmospheric_pressure = metadata.iloc[idx,0]
            density = metadata.iloc[idx,1]
            Altitude = metadata.iloc[idx,2]
            meta_Temperature = metadata.iloc[idx,3]
            pressure = metadata.iloc[idx,4]
            meta_conductivity = metadata.iloc[idx,5]
            salinity = metadata.iloc[idx,6]
            conductivity_specific = metadata.iloc[idx,7]


            # Read Data
            df = pd.read_excel(open(infile, 'rb'),sheet_name = 'Data',skiprows=1)
            df.columns = ["Time", "PAR"]

            # Replace date and time by timestamp
            df["Timestamp"] = [s.timestamp() for s in df["Time"]]

            df = df[["Timestamp", "PAR"]]

            # Sort by date and reset index
            df = df.sort_values(by=['Timestamp'])
            self.df = df.reset_index(inplace=False, drop=True)

            # Sensor data, depth and point
            self.meta_data = {
                "sensor": numsensor,
                "atmospheric_pressure" : atmospheric_pressure,
                "density": density,
                "Altitude": Altitude,
                "meta_Temperature": meta_Temperature,
                "pressure": pressure,
                "meta_conductivity": meta_conductivity,
                "salinity": salinity,
                "conductivity_specific": conductivity_specific,
                "instrument": "rbr",
                "depth": self.depthdict[self.folder],
                "point": self.folder
            }

            # Store variables
            self.time = df["Timestamp"]  # need to be checked
            self.depth = self.meta_data["depth"]
            self.PAR = df["PAR"]

            # Create vector quality for data quality flags at Level 1A
            self.PAR_qual = np.zeros(len(self.PAR))

        # Quality flags in Level 0

    def quality_flags_level0(self, save=False):
        ## Perform quality assurance tests on data
        ## Code quality: 0 = nothing to report, 1 = more investigation

        # First filter: at the beginning and the end of the vector

        # negative values = 0:
        ind = np.where(self.PAR < 0)[0]
        if len(ind) != 0:
            self.PAR[ind] = 0

        vec = np.copy(self.PAR)
        vec_qual = np.copy(self.PAR_qual)


        maintenance = pd.read_csv(self.root + "/scripts/PAR_maintenance.csv", sep=";")
        indanom = []  # index anomaly
        for a in range(0, len(maintenance["debut"])):
            ind = np.where(self.time >= maintenance["debut"][a])[0]
            ind1 = np.where(self.time <= maintenance["end"][a])[0]
            ind2 = np.intersect1d(ind, ind1)

            if len(ind2) != 0:
                indanom = np.append(indanom, ind2)
                indanom = np.unique(indanom,axis=0)

        # First filter:
        # remove maintenance dates:
        if len(indanom) != 0:
            vec_qual[indanom.astype(int)] = 1

        # nan values
        ind = np.where(np.isnan(vec)==True)[0]
        if len(ind) != 0:
            vec[ind] = 10000

        #Remove values > 10000 μmol s-1 m-2
        ind = np.where(vec>=10000)[0]
        if len(ind) !=0:
            vec_qual[ind] = 1

        # first n values and last n values
        # first n values
        n = 30
        y_sub = np.copy(vec[0:n])
        y_qual = np.zeros(len(y_sub))
        y_qual = rm_outliers.rm_IQR(y_sub, y_qual, factor=3)

        if sum(y_qual)>0:
            vec_qual[0:n][y_qual] = 1  # update vec_qual

        # Last n values
        y_sub = np.copy(vec[-n:])
        y_qual = np.zeros(len(y_sub))

        y_qual = rm_outliers.rm_IQR(y_sub, y_qual, factor=3)
        if sum(y_qual)>0:
            vec_qual[-n:][y_qual] = 1  # update vec_qual

        #Remove threshold values error
        vdepth = int(np.ravel(self.depth)[0])

        valid_range = {
            "50": [0, 2000],
            "250": [0, 2000],
            "500": [0, 1500],
            "1000": [0, 1000],
            "2000": [0, 30],
            "3000": [0, 10]}

        # Check if variable is within valid range
        indlow = np.where(vec<valid_range[str(vdepth)][0])
        indhigh = np.where(vec>=valid_range[str(vdepth)][1])
        vec_qual[indlow] = 1
        vec_qual[indhigh] = 1

        # Update vector PAR_qual
        self.PAR_qual = np.copy(vec_qual)

    def quality_flags_v2(self):
        """
        Simplified QC for PAR: range test + maintenance mask only.

        PAR is highly variable by nature (diurnal cycle, clouds, waves),
        so statistical filters are inappropriate. Only physically impossible
        values and maintenance periods are flagged.

        Flag bits (stored in PAR_qc integer array):
            bit 0 (1):  Range test         — hard physical limits (depth-specific)
            bit 1 (2):  Maintenance mask   — known deployment windows
        """
        import yaml

        # Load QC configuration
        config_obj = get_config()
        qc_config_path = os.path.join(str(config_obj.scripts_path), 'qc_config.yaml')
        with open(qc_config_path, 'r') as f:
            qc_config = yaml.safe_load(f)
        par_config = qc_config.get('PAR', {})

        # Prepare inputs
        values = np.copy(self.PAR).astype(float)
        # Clamp negative values to 0 (same as original)
        values[values < 0] = 0.0
        timestamps = np.asarray(self.time, dtype=float)
        depth_cm = int(np.ravel(self.depth)[0])

        # Load maintenance data
        try:
            maintenance_df = pd.read_csv(
                str(config_obj.maintenance_files['par_maintenance']), sep=";"
            )
        except Exception:
            maintenance_df = None

        # Run PAR tests (range + maintenance only)
        qc_flags, qual, flag_details = apply_par_tests(
            values, timestamps, depth_cm,
            maintenance_df, par_config
        )

        # Store results
        self.PAR_qc = qc_flags                        # bitfield (new)
        self.PAR_qual = qual.astype(float)             # 0/1 backward compatible
        self.PAR_flag_details = flag_details           # per-test masks

        # Print summary
        depth_str = f"{depth_cm}cm"
        logger.info("QC v2 - PAR at %s:", depth_str)
        summarize_flags(qc_flags)

class Temperature_data(ppmooring):
    # Level 1A files of this sensor: folder name pattern, variables stored, and the
    # prefix of the key used in self.outfolders
    level1A_folder_pattern = '*Temperature*'
    level1A_variables = ("Temp", "Temp_qual")
    level1A_folder_prefix = ""

    ## Class for the Temperature data

    def __init__(self):
        ## Set attributes for this class
        self.root = str(ppmooring.get_config().ppmooring_root)
        # Output folders
        self.outfolders = {
            "T_1": "/LexplorePPMooringTemperature_50cm",
            "T_2": "/LexplorePPMooringTemperature_250cm",
            "T_3": "/LexplorePPMooringTemperature_500cm",
            "T_4": "/LexplorePPMooringTemperature_750cm",
            "T_5": "/LexplorePPMooringTemperature_1000cm",
            "T_6": "/LexplorePPMooringTemperature_1250cm",
            "T_7": "/LexplorePPMooringTemperature_1500cm",
            "T_8": "/LexplorePPMooringTemperature_1750cm",
            "T_9": "/LexplorePPMooringTemperature_2000cm",
            "T_10": "/LexplorePPMooringTemperature_2250cm",
            "T_11": "/LexplorePPMooringTemperature_2500cm",
            "T_12": "/LexplorePPMooringTemperature_2750cm",
            "T_13": "/LexplorePPMooringTemperature_3000cm",
            "50cm": "/LexplorePPMooringTemperature_50cm",
            "250cm": "/LexplorePPMooringTemperature_250cm",
            "500cm": "/LexplorePPMooringTemperature_500cm",
            "750cm": "/LexplorePPMooringTemperature_750cm",
            "1000cm": "/LexplorePPMooringTemperature_1000cm",
            "1250cm": "/LexplorePPMooringTemperature_1250cm",
            "1500cm": "/LexplorePPMooringTemperature_1500cm",
            "1750cm": "/LexplorePPMooringTemperature_1750cm",
            "2000cm": "/LexplorePPMooringTemperature_2000cm",
            "2250cm": "/LexplorePPMooringTemperature_2250cm",
            "2500cm": "/LexplorePPMooringTemperature_2500cm",
            "2750cm": "/LexplorePPMooringTemperature_2750cm",
            "3000cm": "/LexplorePPMooringTemperature_3000cm"
        }
        self.Level2_outfolder = "/LexplorePPMooringTemperature"

        # Matrix depth - folder
        self.depthdict = {
            "T_1":50,
            "T_2":250,
            "T_3":500,
            "T_4":750,
            "T_5":1000,
            "T_6":1250,
            "T_7":1500,
            "T_8":1750,
            "T_9":2000,
            "T_10":2250,
            "T_11":2500,
            "T_12":2750,
            "T_13":3000
        }

        # General attributes
        self.att = {
            "institution": "EPFL",
            "source": "PPMooring DO",
            "references": "LéXPLORE commun instruments Viet Tran Khac <viet.tran-khac@inrae.fr>",
            "history": "See history on Renku",
            "conventions": "CF 1.7",
            "comment": "Water temperatures from PPMooring M2 on Lexplore Platform in Lake Geneva",
            "title": "Lexplore PPMooring Temperature"
        }

        # Dimension attributes
        self.dim_dict = {
            'time': {'dim_name':'time', 'dim_size':None}
        }

        # Coordinate variable attributes
        self.dim_var_dict = {
            'time': {'var_name':'time', 'dim':('time',), 'unit':'seconds since 1970-01-01 00:00:00', 'longname':'time'}
        }

        # Level 1A - Variable attributes
        self.L1A_dict = {
            'depth': {'var_name':'depth', 'unit':'cm', 'longname':'depth'},
            'Temp': {'var_name':'Temp', 'dim':('time',), 'unit':'degC', 'longname':'temperature'},
            'Temp_qual': {'var_name': 'Temp_qual', 'dim': ('time',),
                          'unit': '0 = nothing to report, 1 = more investigation',
                          'longname': 'temperature quality flag'}
        }

        # Level 1B - Quality flag attributes
        self.L1B_dict = {
            'Temp_qual': {'var_name':'Temp_qual', 'dim':('time',), 'unit':'0 = nothing to report, 1 = more investigation', 'longname':'temperature quality flag'}
        }

        # Level 2 - Variable attributes
        self.L2_dict = {
            'Temp': {'var_name':'Temp', 'dim':('time',), 'unit':'degC', 'longname':'temperature'}
        }

    def read_data_txt(self, infile):
        ## Read .txt data

        # Store infile variable
        self.infile = infile

        # Get input folder
        self.folder = self.infile.rsplit("/")[-3]

        # Import metadata and get indicies of metadata
        with open(self.infile, encoding="utf8", errors='ignore') as f:
            for idx, line in enumerate(f):
                if 'Serial' in line:
                    idx_serial = idx
                if 'Channel[1].calibration' in line:
                    idx_calibr = idx
                if 'NumberOfSamples' in line:
                    idx_meta = idx + 1
                    break
            f.seek(0)
            lines = f.readlines()
            meta_data = lines[0:idx_meta]

        # Extract meta-data
        numsensor = meta_data[idx_serial].rsplit("=")[1].rstrip("\n")
        numcalibration = meta_data[idx_calibr].rsplit("=")[1].rstrip("\n")

        # Read data
        # Handle encoding issues with temperature files using fallback encodings
        try:
            df = pd.read_csv(self.infile, sep=r'\s+', header=None, skiprows=idx_meta + 2, encoding='utf-8')
        except UnicodeDecodeError:
            try:
                df = pd.read_csv(self.infile, sep=r'\s+', header=None, skiprows=idx_meta + 2, encoding='latin1')
            except Exception:
                df = pd.read_csv(self.infile, sep=r'\s+', header=None, skiprows=idx_meta + 2, encoding='utf-8', errors='ignore')
        df.columns = ["Date", "Time", "Temp"]

        # Replace date and time by timestamp
        datehour = [df["Date"] + " " + df["Time"]]
        datehour = datehour[0].apply(lambda x: datetime.strptime(x, '%d-%b-%Y %H:%M:%S.%f'))
        df["Timestamp"] = [datetime.timestamp(s) for s in datehour]

        df = df[["Timestamp", "Temp"]]

        # Sort by date and reset index
        df = df.sort_values(by=['Timestamp'])
        self.df = df.reset_index(inplace=False, drop=True)

        # Sensor data, depth and point
        self.meta_data = {
            "sensor": numsensor,
            "calibration": numcalibration,
            "instrument": "rbr",
            "depth": self.depthdict[self.folder],
            "point": self.folder
        }

        # Store variables
        self.time = df["Timestamp"]  # need to be checked
        self.depth = self.meta_data["depth"]
        self.Temp = df["Temp"]

        # Create vector quality for data quality flags at Level 1A
        self.Temp_qual = np.zeros(len(self.Temp))

    def read_data_csv(self, infile):
        ## Read .txt data

        # Store infile variable
        self.infile = infile

        # Get input folder
        self.folder = self.infile.rsplit("/")[-3]

        # Import metadata and get indicies of metadata
        with open(self.infile, encoding="utf8", errors='ignore') as f:
            for idx, line in enumerate(f):
                if 'Source Device' in line:
                    idx_srcdev = idx
                    break
            f.seek(0)
            lines = f.readlines()
            meta_data = lines[0:8]

        # Extract meta-data
        numsensor = meta_data[idx_srcdev].rsplit("-")[-1].rstrip("\n")
        numcalibration = -99

        # Read data
        # Handle encoding issues with temperature files using fallback encodings
        try:
            df = pd.read_csv(self.infile, sep=",", header=None, skiprows=8, encoding='utf-8')
        except UnicodeDecodeError:
            try:
                df = pd.read_csv(self.infile, sep=",", header=None, skiprows=8, encoding='latin1')
            except Exception:
                df = pd.read_csv(self.infile, sep=",", header=None, skiprows=8, encoding='utf-8', errors='ignore')

        # Some files has ADC column, others do not have
        if len(df.columns) == 4:
            df = df.drop(columns=3)
        df.columns = ["Date", "Time", "Temp"]

        # Replace date and time by timestamp
        datehour = [df["Date"] + " " + df["Time"]]
        datehour = datehour[0].apply(lambda x: datetime.strptime(x, '%Y-%m-%d %H:%M:%S'))
        df["Timestamp"] = [datetime.timestamp(s) for s in datehour]

        df = df[["Timestamp", "Temp"]]

        # Sort by date and reset index
        df = df.sort_values(by=['Timestamp'])
        self.df = df.reset_index(inplace=False, drop=True)

        # Sensor data, depth and point
        self.meta_data = {
            "sensor": numsensor,
            "calibration": numcalibration,
            "instrument": "minidot",
            "depth": self.depthdict[self.folder],
            "point": self.folder
        }

        # Store variables
        self.time = df["Timestamp"]  # need to be checked
        self.depth = self.meta_data["depth"]
        self.Temp = df["Temp"]

        # Create vector quality for data quality flags at Level 1A
        self.Temp_qual = np.zeros(len(self.Temp))

    def quality_flags_v2(self):
        """
        Simplified QC for Temperature: range, maintenance, spike, gradient, flat-line.

        Flag bits (stored in Temp_qc integer array):
            bit 0 (1):  Range test         — depth-specific physical limits
            bit 1 (2):  Maintenance mask   — known deployment windows
            bit 3 (8):  Spike test         — short deviations returning to baseline
            bit 4 (16): Gradient test      — rate of change exceeds physical max
            bit 5 (32): Flat-line test     — consecutive identical values
        """
        import yaml

        # Load QC configuration
        config_obj = get_config()
        qc_config_path = os.path.join(str(config_obj.scripts_path), 'qc_config.yaml')
        with open(qc_config_path, 'r') as f:
            qc_config = yaml.safe_load(f)
        temp_config = qc_config.get('Temperature', {})

        # Prepare inputs
        values = np.copy(self.Temp).astype(float)
        timestamps = np.asarray(self.time, dtype=float)
        depth_cm = int(np.ravel(self.depth)[0])

        # Load maintenance data
        try:
            maintenance_df = pd.read_csv(
                str(config_obj.maintenance_files['par_maintenance']), sep=";"
            )
        except Exception:
            maintenance_df = None

        # Try self.maintenance if available (set by main_ppmooring.py)
        if maintenance_df is None and hasattr(self, 'maintenance'):
            maintenance_df = self.maintenance

        # Run Temperature tests
        qc_flags, qual, flag_details = apply_temperature_tests(
            values, timestamps, depth_cm,
            maintenance_df, temp_config
        )

        # Store results
        self.Temp_qc = qc_flags                        # bitfield (new)
        self.Temp_qual = qual.astype(float)             # 0/1 backward compatible
        self.Temp_flag_details = flag_details           # per-test masks

        # Print summary
        depth_str = f"{depth_cm}cm"
        logger.info("QC v2 - Temp at %s:", depth_str)
        summarize_flags(qc_flags)

class rm_outliers(ppmooring):

    def rm_IQR(vecX,vecX_qual,factor=3):
        """
            Remove values outside Interquartile Range (IQR) for timeseries data.

            Parameters:
                vecX (np.array): Data array to which to apply the quality assurance
                vecX_qual (np.array): An array of bools where True means non-trusted data before this filter
                factor (np.int): threshold for outlier labelling rule

            Returns:
                y_qual (np.array): An array of bools where True means non-trusted data after this outlier detection
        """
        y = np.copy(vecX)
        y_qual = np.copy(vecX_qual)
        q75 = np.quantile(y, 0.75)
        q25 = np.quantile(y, 0.25)
        IQR = q75 - q25
        if IQR != 0:
            outsup = q75 + factor * IQR
            outinf = q25 - factor * IQR

            idx1 = np.where(y >= outsup)[0]
            idx2 = np.where(y <= outinf)[0]

            idx0 = np.r_[idx1, idx2]
        else:
            idx0 = np.array([])

        if len(idx0) == 0:
            # print("No outliers detected - IQR")
            idx = idx0
        else:
            y_qual[idx0] = y_qual[idx0] + 1
            y_qual = np.array(y_qual,dtype=bool)

            #idx = np.where(y_qual >= 1)[0]
        return y_qual

