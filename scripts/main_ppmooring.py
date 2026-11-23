# -*- coding: utf-8 -*-
"""Process PPMooring sensor data: Level 0 -> Level 1A -> Level 1B -> Level 2.

Usage:
    python main_ppmooring.py                      all sensors
    python main_ppmooring.py --sensors DO PAR     selected sensors only
    python main_ppmooring.py --log-level DEBUG

A file that fails is logged and skipped; the other files are still processed. A failed
file is not recorded as processed, so the next run tries it again. The exit code is 1
when at least one file failed.

Files are processed in chronological order within each depth folder, to prevent
timestamp conflicts.
"""
# Notes to improve: metadata from Level0 was not stored in Level1A. Not necessary for the moment.
from __future__ import annotations

import time
_start_time = time.time()

import argparse
import logging
import os
import re
import sys
from dataclasses import dataclass
from datetime import datetime
from functools import partial
from pathlib import Path
from typing import Callable, Dict, Iterable, List, Optional, Sequence, Tuple

import pandas as pd

try:  # Use tqdm if available for fast progress bars
    from tqdm import tqdm  # type: ignore

    def iter_progress(iterable, desc=""):
        return tqdm(iterable, desc=desc, unit="file", dynamic_ncols=True, leave=False)
except Exception:  # Graceful fallback without external deps
    def iter_progress(iterable, desc=""):
        return iterable

# Make config.py / log_optimizer.py (parent folder) and the scripts folder importable
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from config import get_config
from log_optimizer import LogRegistryManager
from ppmooring import CalibrationCorrector, DO_data, PAR_data, Temperature_data

logger = logging.getLogger("main_ppmooring")

SENSOR_CHOICES = ("PAR", "DO", "Temperature")  # also the processing order
RUN_LOG_NAME = "pipeline_run.log"
CALIBRATION_FILE = "data/calibration_data/Calibration-sensore-OxyPAR-compilation.xlsx"


def dedupe_preserve_order(items: Iterable[str]) -> List[str]:
    seen = set()
    out = []
    for x in items:
        if x not in seen:
            seen.add(x)
            out.append(x)
    return out


class ChronologicalFileProcessor:
    """
    Manages chronological ordering of sensor files for processing.
    
    Extracts timestamps from filenames and file content to ensure proper
    temporal sequence during Level 0 -> Level 1A -> Level 1B -> Level 2 processing.
    """
    
    def __init__(self):
        self.config = get_config()
        
        # Filename date extraction patterns
        self.date_patterns = [
            # YYYY-MM-DD HHMMSSZ.txt format (e.g., 2018-09-17 155200Z.txt)
            r'(\d{4})-(\d{2})-(\d{2})\s+(\d{6})Z',  # 2018-09-17 155200Z
            # YYYY-MM-DD formats
            r'(\d{4})-(\d{2})-(\d{2})',             # 2018-09-17
            # YYYYMMDD formats  
            r'(\d{8})',           # 20181120 (8 digits)
            r'(\d{8})_\d{4}',     # 20181120_1039
            r'(\d{8})-',          # 20181120-
            r'_(\d{8})',          # _20181120
            r'-(\d{8})',          # -20181120
            r'(\d{4})(\d{2})(\d{2})'  # YYYYMMDD split
        ]
    
    def extract_date_from_filename(self, filename: str) -> Optional[datetime]:
        """
        Extract date from filename using various patterns.
        
        Args:
            filename: Name of the file
            
        Returns:
            datetime object if date found, None otherwise
        """
        for pattern in self.date_patterns:
            match = re.search(pattern, filename)
            if match:
                try:
                    groups = match.groups()
                    if len(groups) == 1:
                        date_str = groups[0]
                        if len(date_str) == 8:
                            return datetime.strptime(date_str, '%Y%m%d')
                    elif len(groups) == 3:
                        year, month, day = groups
                        return datetime(int(year), int(month), int(day))
                    elif len(groups) == 4:
                        # YYYY-MM-DD HHMMSSZ format (2018-09-17 155200Z)
                        year, month, day, time_str = groups
                        # Parse time: 155200 -> 15:52:00
                        if len(time_str) == 6:
                            hour = int(time_str[:2])
                            minute = int(time_str[2:4])
                            second = int(time_str[4:6])
                            return datetime(int(year), int(month), int(day), hour, minute, second)
                        else:
                            # If time parsing fails, use just the date
                            return datetime(int(year), int(month), int(day))
                except (ValueError, TypeError):
                    continue
        return None
    
    def extract_date_from_content(self, file_path: Path) -> Optional[datetime]:
        """
        Extract earliest timestamp from file content.
        
        Args:
            file_path: Path to the data file
            
        Returns:
            datetime of earliest data point, None if extraction fails
        """
        try:
            if file_path.suffix.lower() == '.csv':
                df = pd.read_csv(file_path, nrows=10)  # Read first few rows
                
                # Look for time columns
                time_columns = ['time', 'Time', 'datetime', 'timestamp', 'DATE_TIME']
                for col in time_columns:
                    if col in df.columns:
                        try:
                            # Try parsing as unix timestamp
                            first_time = pd.to_datetime(df[col].iloc[0], unit='s')
                            return first_time.to_pydatetime()
                        except:
                            try:
                                # Try parsing as datetime string
                                first_time = pd.to_datetime(df[col].iloc[0])
                                return first_time.to_pydatetime()
                            except:
                                continue
                                
            elif file_path.suffix.lower() in ['.xlsx', '.xls']:
                df = pd.read_excel(file_path, nrows=10)
                
                # Similar logic for Excel files
                time_columns = ['time', 'Time', 'datetime', 'timestamp', 'DATE_TIME']
                for col in time_columns:
                    if col in df.columns:
                        try:
                            first_time = pd.to_datetime(df[col].iloc[0])
                            return first_time.to_pydatetime()
                        except:
                            continue
                            
        except Exception as e:
            logger.warning("Could not extract date from content of %s: %s", file_path, e)
        
        return None
    
    def get_chronological_file_list_by_depth(self, sensor_type: str, level: str = "Level0") -> Dict[str, List[Tuple[Path, datetime, str]]]:
        """
        Get list of files for a sensor type organized by depth, each in chronological order.
        
        Args:
            sensor_type: 'DO', 'RBR_PAR', or 'Temperature'
            level: Data level to process ('Level0', 'Level1A', etc.)
            
        Returns:
            Dictionary where keys are depth folder names and values are 
            lists of (file_path, timestamp, sort_method) tuples in chronological order
        """
        sensor_dir = self.config.get_sensor_data_path(level, sensor_type)
        
        if not sensor_dir.exists():
            logger.warning("Sensor directory does not exist: %s", sensor_dir)
            return {}
        
        depth_files = {}
        
        # Find depth directories
        depth_dirs = []
        if sensor_type == "RBR_PAR":
            # PAR files have depth directories like RBR_PAR_1000cm
            depth_dirs = [d for d in sensor_dir.glob("*cm") if d.is_dir()]
        elif sensor_type == "Temperature":
            # Temperature has depth directories like T_1, T_10, etc.
            depth_dirs = [d for d in sensor_dir.glob("T_*") if d.is_dir()]
            # Also include Thermistor_string if it exists
            thermistor_dir = sensor_dir / "Thermistor_string"
            if thermistor_dir.exists() and thermistor_dir.is_dir():
                depth_dirs.append(thermistor_dir)
        else:
            # DO has depth directories like DO_50cm
            depth_dirs = [d for d in sensor_dir.glob("*cm") if d.is_dir()]
        
        # Handle files directly in the root sensor directory (not in depth subfolders)
        root_files = [f for f in sensor_dir.glob("*") if f.is_file() and f.suffix.lower() in ['.txt', '.csv', '.xlsx', '.xls']]
        if root_files:
            depth_dirs.append(sensor_dir)  # Add root directory as "depth"
        
        if not depth_dirs:
            logger.warning("No depth directories or files found in %s", sensor_dir)
            return {}
        
        logger.info("Found %d depth folders for %s", len(depth_dirs), sensor_type)
        
        for depth_dir in depth_dirs:
            depth_name = depth_dir.name if depth_dir != sensor_dir else "root"
            files_with_timestamps = []
            
            # Find all data files in this depth folder
            data_files = []
            if sensor_type == "RBR_PAR" and depth_dir != sensor_dir:
                # PAR files are in depth-specific subdirectories with conv folder
                conv_dir = depth_dir / "conv"
                if conv_dir.exists():
                    data_files.extend(conv_dir.glob("*"))
            elif depth_dir == sensor_dir:
                # Root directory - only get files directly in root, not in subdirectories
                data_files.extend(depth_dir.glob("*"))
            else:
                # Depth-specific folders - search recursively within this depth only
                data_files.extend(depth_dir.rglob("*"))
            
            # Filter for actual data files
            data_files = [f for f in data_files if f.is_file() and f.suffix.lower() in ['.txt', '.csv', '.xlsx', '.xls']]
            
            logger.info("  %s: Found %d files", depth_name, len(data_files))
            
            for file_path in data_files:
                # Try to extract date from filename first (faster)
                timestamp = self.extract_date_from_filename(file_path.name)
                sort_method = "filename"
                
                # If filename parsing failed, try content parsing
                if timestamp is None:
                    timestamp = self.extract_date_from_content(file_path)
                    sort_method = "content"
                
                # If still no timestamp, skip this file with error
                if timestamp is None:
                    logger.error("Could not extract timestamp from %s - skipping file "
                                 "(filename patterns not matched and no readable time data in file content)",
                                 file_path.name)
                    continue  # Skip this file entirely
                
                files_with_timestamps.append((file_path, timestamp, sort_method))
            
            # Sort by timestamp (chronological order) within this depth
            files_with_timestamps.sort(key=lambda x: x[1])
            
            if files_with_timestamps:
                depth_files[depth_name] = files_with_timestamps
        
        return depth_files
    
    def validate_chronological_order_by_depth(self, depth_files: Dict[str, List[Tuple[Path, datetime, str]]]) -> Dict[str, Dict[str, any]]:
        """
        Validate the chronological ordering and identify potential issues for each depth.
        
        Args:
            depth_files: Dictionary with depth names as keys and file lists as values
            
        Returns:
            Dictionary with depth names as keys and validation results as values
        """
        if not depth_files:
            return {}
        
        validation_results = {}
        
        for depth_name, files_list in depth_files.items():
            validation_results[depth_name] = self._validate_single_depth(files_list)
        
        return validation_results
    
    def _validate_single_depth(self, files_list: List[Tuple[Path, datetime, str]]) -> Dict[str, any]:
        """
        Validate chronological order for a single depth folder.
        
        Args:
            files_list: List of (file_path, timestamp, sort_method) tuples
            
        Returns:
            Dictionary with validation results and warnings
        """
        if not files_list:
            return {"valid": True, "warnings": [], "issues": []}
        
        warnings = []
        issues = []
        
        # Check for duplicated timestamps
        timestamps = [f[1] for f in files_list]
        duplicates = {}
        for i, ts in enumerate(timestamps):
            if timestamps.count(ts) > 1:
                if ts not in duplicates:
                    duplicates[ts] = []
                duplicates[ts].append(files_list[i][0].name)
        
        if duplicates:
            for ts, filenames in duplicates.items():
                issues.append(f"Duplicate timestamp {ts.strftime('%Y-%m-%d')}: {', '.join(filenames)}")
        
        # Check for large time gaps (potential missing data)
        large_gaps = []
        for i in range(len(files_list) - 1):
            current_ts = files_list[i][1]
            next_ts = files_list[i + 1][1]
            gap_days = (next_ts - current_ts).days
            
            if gap_days > 90:  # More than 3 months
                large_gaps.append(f"Large gap: {gap_days} days between {current_ts.strftime('%Y-%m-%d')} and {next_ts.strftime('%Y-%m-%d')}")
        
        if large_gaps:
            warnings.extend(large_gaps[:5])  # Limit to first 5 warnings
        
        # Check for very close timestamps (potential overlapping data)
        close_timestamps = []
        for i in range(len(files_list) - 1):
            current_ts = files_list[i][1]
            next_ts = files_list[i + 1][1]
            gap_hours = (next_ts - current_ts).total_seconds() / 3600
            
            if 0 < gap_hours < 24:  # Less than 24 hours apart
                close_timestamps.append(f"Close timestamps: {gap_hours:.1f} hours between {files_list[i][0].name} and {files_list[i+1][0].name}")
        
        if close_timestamps:
            warnings.extend(close_timestamps[:3])  # Limit to first 3 warnings
        
        return {
            "valid": len(issues) == 0,
            "warnings": warnings,
            "issues": issues,
            "total_files": len(files_list),
            "date_range": (files_list[0][1], files_list[-1][1]) if files_list else None,
            "sort_methods": {method: sum(1 for f in files_list if f[2] == method) for method in ["filename", "content"]}
        }
    
    def get_processing_order_for_sensor(self, sensor_type: str) -> Dict[str, any]:
        """
        Get recommended processing order for a specific sensor type, organized by depth.
        
        Args:
            sensor_type: 'DO', 'RBR_PAR', or 'Temperature'
            
        Returns:
            Dictionary with processing order and validation results organized by depth
        """
        depth_files = self.get_chronological_file_list_by_depth(sensor_type)
        validation_results = self.validate_chronological_order_by_depth(depth_files)
        
        # Create recommended order for each depth
        recommended_order_by_depth = {}
        for depth_name, files_list in depth_files.items():
            recommended_order_by_depth[depth_name] = [f[0] for f in files_list]  # Just the file paths
        
        return {
            "sensor_type": sensor_type,
            "depth_files": depth_files,
            "validation_by_depth": validation_results,
            "recommended_order_by_depth": recommended_order_by_depth
        }


class Run:
    """Shared state of one pipeline run."""

    def __init__(self) -> None:
        self.config = get_config()
        # Ensure all required directories exist
        self.config.ensure_directories(create_missing=True)
        # Optimized log system: which files are already processed
        self.log_manager = LogRegistryManager(str(self.config.log_path))
        # Chronological processor for proper temporal ordering
        self.chrono_processor = ChronologicalFileProcessor()
        self.calibration_corrector: Optional[CalibrationCorrector] = None
        self.directories = {
            "Level1A_dir": str(self.config.data_paths["Level1A"]),
            "Level1B_dir": str(self.config.data_paths["Level1B"]),
            "Level2_dir": str(self.config.data_paths["Level2"]),
        }
        # (sensor, level, file, error message) for every file that failed
        self.failed: List[Tuple[str, str, str, str]] = []


def setup_logging(level: str, log_dir: Path) -> None:
    """Send log messages to the console and to <log_dir>/pipeline_run.log."""
    os.makedirs(log_dir, exist_ok=True)
    formatter = logging.Formatter("%(asctime)s %(levelname)-7s %(name)s: %(message)s")
    console = logging.StreamHandler(sys.stdout)
    console.setFormatter(formatter)
    logfile = logging.FileHandler(os.path.join(str(log_dir), RUN_LOG_NAME), encoding="utf-8")
    logfile.setFormatter(formatter)
    root = logging.getLogger()
    root.setLevel(level)
    root.addHandler(console)
    root.addHandler(logfile)


def init_calibration(run: Run) -> None:
    """Load the DO calibration corrector, or continue without it."""
    calibration_file_path = run.config.ppmooring_root / CALIBRATION_FILE
    if calibration_file_path.exists():
        try:
            run.calibration_corrector = CalibrationCorrector(str(calibration_file_path))
            logger.info("Calibration corrector initialized: %d records loaded",
                        len(run.calibration_corrector.calibration_data))
        except Exception as e:
            run.calibration_corrector = None
            logger.warning("Failed to initialize calibration corrector: %s", e)
    else:
        run.calibration_corrector = None
        logger.warning("Calibration data not found: %s", calibration_file_path)


def delete_ds_store(run: Run) -> None:
    """Delete .DS_Store files under the data folder."""
    for root1, dirs, files in os.walk(run.config.data_root):
        for file in files:
            if file.endswith('.DS_Store'):
                os.remove(os.path.join(root1, file))


def list_level0_files(run: Run, sensor_type: str, label: str, remove_small: bool = False) -> List[str]:
    """Level 0 files of one sensor type, in chronological order within each depth."""
    logger.info("Getting %s files in chronological order...", label)
    by_depth = run.chrono_processor.get_chronological_file_list_by_depth(sensor_type, "Level0")

    # Extract just the file paths for processing from all depth directories
    files = []
    for depth_dir, file_list in by_depth.items():
        files.extend([str(file_info[0]) for file_info in file_list])

    # Deduplicate while preserving order (avoid duplicated run)
    files = dedupe_preserve_order(files)

    if remove_small:
        # Check and delete empty files (<= 1000 bytes)
        for infile in files.copy():
            if os.stat(infile).st_size <= 1000:
                os.remove(infile)
                files.remove(infile)

    logger.info("Processing %d %s files in chronological order", len(files), label)
    return files


def process_files(run: Run, files: Sequence[str], sensor: str, level: str,
                  process_one: Callable[[Run, str], Optional[bool]], desc: Optional[str] = None) -> None:
    """Run ``process_one(run, infile)`` on every file not processed yet.

    ``process_one`` returns False when it skipped the file on purpose; such a file is
    not recorded as processed. An exception is logged, the file is added to
    ``run.failed`` and the loop goes on with the next file.
    """
    for infile in (iter_progress(files, desc) if desc else files):
        infile = infile.replace(os.sep, '/')

        # Check if file already processed using optimized log system
        if run.log_manager.is_processed(infile, sensor, level):
            continue

        try:
            if process_one(run, infile) is False:
                continue
            # Mark file as successfully processed
            run.log_manager.mark_processed(infile, sensor, level)
        except Exception as e:
            logger.exception("Failed %s %s: %s", sensor, level, infile)
            run.failed.append((sensor, level, infile, f"{type(e).__name__}: {e}"))


def update_maintenance_csv(run: Run) -> None:
    """Rebuild PAR_maintenance.csv from the deployment boundaries of the calibration file.

    Each deployment boundary date gets a +/-12h maintenance window.
    """
    calibration_xlsx = run.config.ppmooring_root / CALIBRATION_FILE
    if not calibration_xlsx.exists():
        logger.warning("Calibration file not found: %s, keeping existing PAR_maintenance.csv",
                       calibration_xlsx)
        return

    cal_df = pd.read_excel(str(calibration_xlsx), sheet_name='Sheet1', header=2)
    cal_df['Beginning deployment'] = pd.to_datetime(cal_df['Beginning deployment'], errors='coerce')
    cal_df['End Deployment'] = pd.to_datetime(cal_df['End Deployment'], errors='coerce')

    # Collect all unique deployment boundary dates
    boundary_dates = sorted(set(
        cal_df['Beginning deployment'].dropna().tolist() +
        cal_df['End Deployment'].dropna().tolist()
    ))

    # Create maintenance windows: +/-12h around each boundary date
    half_day = pd.Timedelta(hours=12)
    windows = [(d - half_day, d + half_day) for d in boundary_dates]

    # Merge overlapping windows
    merged = [windows[0]]
    for start, end in windows[1:]:
        if start <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((start, end))

    # Convert to unix timestamps and save as debut;end CSV
    debut_ts = [w[0].timestamp() for w in merged]
    end_ts = [w[1].timestamp() for w in merged]
    df_maint = pd.DataFrame({"debut": debut_ts, "end": end_ts})

    filename = run.config.maintenance_files['par_maintenance']
    df_maint.to_csv(filename, sep=";", header=True, index=False)
    logger.info("Maintenance CSV updated: %d windows from %d deployment dates",
                len(merged), len(boundary_dates))


## Sensor descriptions ######################################

def _unsupported_temperature_depth(infile: str) -> bool:
    # Only T_1 to T_13 are processed; T_14 and beyond and Thermistor_string are skipped
    return any(f'/T_{i}/' in infile for i in range(14, 22)) or '/Thermistor_string/' in infile


@dataclass(frozen=True)
class SensorSpec:
    """What differs between sensors in the two processing stages."""
    name: str              # name on the command line and in messages
    log_key: str           # key of the processed-file logs
    cls: type              # sensor class of ppmooring.py
    level0_type: str       # folder under data/Level0
    readers: tuple         # (matches(infile), reader method), first match wins
    level1A_label: str     # label in Level 1A file names
    level2_label: str      # label in Level 2 file names
    remove_small: bool = False            # delete Level 0 files <= 1000 bytes before reading
    level0_qc: bool = False               # run quality_flags_level0() before saving Level 1A
    calibrate: bool = False               # apply the calibration correction before saving Level 1A
    skip_unknown: bool = False            # no reader matches: skip with a warning instead of failing
    skip_file: Optional[Callable[[str], bool]] = None  # files left out on purpose


SENSORS = {
    "PAR": SensorSpec(
        name="PAR", log_key="PAR", cls=PAR_data, level0_type="RBR_PAR",
        readers=((lambda f: True, PAR_data.read_data),),
        level1A_label="PAR", level2_label="PAR",
        remove_small=True, level0_qc=True,
    ),
    "DO": SensorSpec(
        name="DO", log_key="DO", cls=DO_data, level0_type="DO",
        readers=(
            (lambda f: "rbr" in f and "txt" in f, DO_data.read_data_rbr),
            (lambda f: "rbr" in f and "xlsx" in f, DO_data.read_data_rbr_excel),
            (lambda f: "exo" in f, DO_data.read_data_exo),
            (lambda f: "minidot" in f, DO_data.read_data_minidot),
        ),
        level1A_label="DO", level2_label="DO",
        remove_small=True, calibrate=True, skip_unknown=True,
    ),
    "Temperature": SensorSpec(
        name="Temperature", log_key="Temp", cls=Temperature_data, level0_type="Temperature",
        readers=(
            (lambda f: ".txt" in f, Temperature_data.read_data_txt),
            (lambda f: ".csv" in f, Temperature_data.read_data_csv),
        ),
        level1A_label="Temp", level2_label="Temp",
        skip_file=_unsupported_temperature_depth,
    ),
}


## Processing steps ######################################

def find_reader(spec: SensorSpec, infile: str) -> Optional[Callable]:
    for matches, reader in spec.readers:
        if matches(infile):
            return reader
    return None


def apply_calibration(run: Run, data: DO_data, infile: str) -> None:
    """Correct DO data with the calibration file. A failure here does not fail the file."""
    try:
        correction_result = run.calibration_corrector.apply_correction(data)
        if correction_result:
            logger.info("Calibration applied at %scm: %s period(s)",
                        correction_result['depth_cm'], correction_result['n_corrections'])
            for corr in correction_result['corrections']:
                logger.info("  Sensor %s: %+.2f%% (%s to %s): %s points",
                            corr['sensor'], corr['correction_percent'],
                            corr['period_start'], corr['period_end'], corr['n_points'])
    except Exception as e:
        logger.warning("Calibration correction failed for %s: %s", os.path.basename(infile), e)


def level0_to_1A(run: Run, infile: str, spec: SensorSpec) -> Optional[bool]:
    """Level 0 --> Level 1A for one file. Return False when the file is skipped on purpose."""
    if spec.skip_file is not None and spec.skip_file(infile):
        logger.info("Skipping file from unsupported depth: %s", infile)
        return False

    reader = find_reader(spec, infile)
    if reader is None:
        if spec.skip_unknown:
            logger.warning("Unrecognized %s file type, skipping: %s", spec.name, infile)
            return False
        raise ValueError(f"No reader for this {spec.name} file type")

    data = spec.cls()
    reader(data, infile)

    # Calibration correction BEFORE saving Level 1A
    if spec.calibrate and run.calibration_corrector is not None:
        apply_calibration(run, data, infile)

    # Set output name and folder for Level 1A
    data.set_output(run.directories["Level1A_dir"], "L1A", spec.level1A_label)
    if spec.level0_qc:
        data.quality_flags_level0()
    # Store Level 1A data to NetCDF and csv
    data.to_NetCDF("L1A")
    data.to_csv("L1A")


def level1A_to_1B_2(run: Run, infile: str, spec: SensorSpec) -> None:
    """Level 1A --> Level 1B --> Level 2 for one Level 1A NetCDF file."""
    data = spec.cls()
    data.read_data_Level1A_nc(infile)
    # Set output name and folder for Level 1B
    data.set_output(run.directories["Level1B_dir"], "L1B", "DO")

    # Level 1A --> Level 1B: QC v2 (independent tests, bit flags)
    data.quality_flags_v2()
    # Store Level 1B data to NetCDF and csv
    data.to_NetCDF("L1A", "L1B")
    data.to_csv("L1A", "L1B")

    # Level 1B --> Level 2: mask bad data, then store on the time grid
    data.mask_data()
    data.Level2_to_NetCDF(run.directories["Level2_dir"], spec.level2_label)


def process_level0(run: Run, spec: SensorSpec) -> None:
    """All Level 0 files of one sensor --> Level 1A."""
    files = list_level0_files(run, spec.level0_type, spec.name, remove_small=spec.remove_small)
    process_files(run, files, spec.log_key, '01A', partial(level0_to_1A, spec=spec),
                  desc=f"{spec.log_key} L0→L1A")


def process_level1(run: Run, spec: SensorSpec) -> None:
    """All Level 1A NetCDF files of one sensor --> Level 1B + Level 2."""
    lister = spec.cls()
    lister.get_files_to_read_Level1A_nc(run.directories["Level1A_dir"])
    logger.info("Processing %d %s L1A NetCDF files -> L1B + L2", len(lister.filestoread), spec.name)
    process_files(run, lister.filestoread, spec.log_key, '1A1B', partial(level1A_to_1B_2, spec=spec))


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Process PPMooring sensor data from Level 0 to Level 2.")
    parser.add_argument("--sensors", nargs="+", choices=SENSOR_CHOICES,
                        default=list(SENSOR_CHOICES), metavar="SENSOR",
                        help="sensors to process: PAR, DO, Temperature (default: all)")
    parser.add_argument("--log-level", default="INFO",
                        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
                        help="amount of detail in the log (default: INFO)")
    return parser.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)

    run = Run()
    setup_logging(args.log_level, run.config.log_path)
    logger.info("Sensors to process: %s", ", ".join(s for s in SENSOR_CHOICES if s in args.sensors))

    init_calibration(run)
    delete_ds_store(run)

    # PAR first: its Level 0 step runs before the maintenance file is rebuilt,
    # as it always did. The maintenance file is used by the QC of all sensors,
    # so it is rebuilt even when PAR is not selected.
    if "PAR" in args.sensors:
        process_level0(run, SENSORS["PAR"])
    update_maintenance_csv(run)
    if "PAR" in args.sensors:
        process_level1(run, SENSORS["PAR"])

    for name in ("DO", "Temperature"):
        if name in args.sensors:
            process_level0(run, SENSORS[name])
            process_level1(run, SENSORS[name])

    # Processing statistics
    stats = run.log_manager.get_statistics()
    logger.info("=== Processing Statistics ===")
    for sensor_level, data in stats.items():
        logger.info("%s: %s files processed", sensor_level, data['processed_count'])

    elapsed = time.time() - _start_time
    if run.failed:
        logger.error("=== %d file(s) failed ===", len(run.failed))
        for sensor, level, infile, error in run.failed:
            logger.error("%s %s: %s -> %s", sensor, level, infile, error)
        logger.error("Processing finished with %d failed file(s) in %.1f min (%.0f s)",
                     len(run.failed), elapsed / 60, elapsed)
        return 1

    logger.info("Processing completed successfully in %.1f min (%.0f s)", elapsed / 60, elapsed)
    return 0


if __name__ == "__main__":
    sys.exit(main())
