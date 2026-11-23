"""
PPMooring Configuration Management

This module provides centralized configuration for the ppmooring sensor data processing pipeline.
It automatically detects the project structure and provides cross-platform path management.
"""

import os
import logging
from pathlib import Path
from typing import Dict, Optional, Any
import yaml

logger = logging.getLogger(__name__)


class ConfigurationError(Exception):
    """Raised when configuration is invalid or missing."""
    pass


class Config:
    """
    Centralized configuration management for ppmooring pipeline.
    
    Automatically detects project root and provides validated paths for all data levels.
    Supports environment-specific overrides via environment variables.
    """
    
    def __init__(self, env: Optional[str] = None, base_path: Optional[Path] = None):
        """
        Initialize configuration.
        
        Args:
            env: Environment name ('development', 'production', 'testing')
            base_path: Override auto-detection of project root
        """
        self.env = env or os.getenv('PPMOORING_ENV', 'development')
        self.base_path = base_path or self._detect_base_path()
        self._validate_structure()
        
    def _detect_base_path(self) -> Path:
        """Auto-detect the project root directory."""
        # Start from current file location
        current = Path(__file__).parent.absolute()
        
        # Look for ppmooring directory structure
        while current != current.parent:  # Not at filesystem root
            if (current / 'ppmooring').exists():
                return current
            current = current.parent
            
        # Fallback: check if we're already inside ppmooring
        current = Path(__file__).parent.absolute()
        if current.name == 'ppmooring':
            return current.parent
            
        # Last resort: use environment variable
        if 'PPMOORING_ROOT' in os.environ:
            return Path(os.environ['PPMOORING_ROOT'])
            
        raise ConfigurationError(
            "Cannot detect ppmooring project root. "
            "Set PPMOORING_ROOT environment variable or ensure proper directory structure."
        )
    
    def _validate_structure(self):
        """Validate that required directories exist."""
        ppmooring_dir = self.base_path / 'ppmooring'
        
        if not ppmooring_dir.exists():
            raise ConfigurationError(
                f"ppmooring directory not found at {ppmooring_dir}. "
                f"Base path: {self.base_path}"
            )
            
        scripts_dir = ppmooring_dir / 'scripts'
        if not scripts_dir.exists():
            logger.warning(f"Scripts directory not found at {scripts_dir}")
    
    @property
    def ppmooring_root(self) -> Path:
        """Path to the ppmooring directory."""
        return self.base_path / 'ppmooring'
    
    @property
    def scripts_path(self) -> Path:
        """Path to the scripts directory."""
        return self.ppmooring_root / 'scripts'
    
    @property
    def data_root(self) -> Path:
        """Path to the data root directory."""
        data_dir = os.getenv('PPMOORING_DATA_DIR')
        if data_dir:
            return Path(data_dir)
        return self.ppmooring_root / 'data'
    
    @property
    def data_paths(self) -> Dict[str, Path]:
        """Dictionary of all data level paths."""
        base_data = self.data_root
        return {
            'Level0': base_data / 'Level0',
            'Level1A': base_data / 'Level1A',
            'Level1B': base_data / 'Level1B', 
            'Level2': base_data / 'Level2',
            'Level2a': base_data / 'Level2a'
        }
    
    @property
    def log_path(self) -> Path:
        """Path to log directory."""
        log_dir = os.getenv('PPMOORING_LOG_DIR')
        if log_dir:
            return Path(log_dir)
        return self.ppmooring_root / 'log'
    
    @property
    def maintenance_files(self) -> Dict[str, Path]:
        """Paths to maintenance and configuration files."""
        scripts = self.scripts_path
        return {
            'par_maintenance': scripts / 'PAR_maintenance.csv',
            'border_values': scripts / 'border_value_lower.csv',
            'out_of_water': scripts / 'out_of_water_periods.csv'
        }
    
    def get_sensor_data_path(self, level: str, sensor_type: str, depth: Optional[str] = None) -> Path:
        """
        Get path for specific sensor data.
        
        Args:
            level: Data level ('Level0', 'Level1A', etc.)
            sensor_type: Sensor type ('DO', 'RBR_PAR', 'Temperature')
            depth: Optional depth specification (e.g., '1000cm', '3000cm')
            
        Returns:
            Path to sensor data directory
        """
        base_path = self.data_paths[level]
        
        if level == 'Level0':
            return base_path / sensor_type
        else:
            # Level1A, Level1B, Level2 use different naming convention
            if sensor_type == 'RBR_PAR':
                sensor_name = f'LexplorePPMooringPAR_{depth}' if depth else 'LexplorePPMooringPAR'
            elif sensor_type == 'DO':
                sensor_name = f'LexplorePPMooringDO_{depth}' if depth else 'LexplorePPMooringDO'
            elif sensor_type == 'Temperature':
                sensor_name = f'LexplorePPMooringTemp_{depth}' if depth else 'LexplorePPMooringTemp'
            else:
                raise ConfigurationError(f"Unknown sensor type: {sensor_type}")
                
            return base_path / sensor_name
    
    def ensure_directories(self, create_missing: bool = True) -> bool:
        """
        Ensure all required directories exist.
        
        Args:
            create_missing: Create directories if they don't exist
            
        Returns:
            True if all directories exist/were created successfully
        """
        required_dirs = [
            self.data_root,
            self.log_path,
            *self.data_paths.values()
        ]
        
        success = True
        for directory in required_dirs:
            if not directory.exists():
                if create_missing:
                    try:
                        directory.mkdir(parents=True, exist_ok=True)
                        logger.info(f"Created directory: {directory}")
                    except OSError as e:
                        logger.error(f"Failed to create directory {directory}: {e}")
                        success = False
                else:
                    logger.warning(f"Required directory does not exist: {directory}")
                    success = False
                    
        return success
    
    def load_yaml_config(self, config_file: str = 'input_python.yaml') -> Dict[str, Any]:
        """
        Load YAML configuration file.
        
        Args:
            config_file: Name of YAML config file in scripts directory
            
        Returns:
            Dictionary of configuration values
        """
        config_path = self.scripts_path / config_file
        
        if not config_path.exists():
            logger.warning(f"Config file not found: {config_path}")
            return {}
            
        try:
            with open(config_path, 'r') as f:
                return yaml.safe_load(f) or {}
        except yaml.YAMLError as e:
            logger.error(f"Error parsing YAML config {config_path}: {e}")
            return {}
        except Exception as e:
            logger.error(f"Error reading config file {config_path}: {e}")
            return {}
    
    def __str__(self) -> str:
        """String representation of configuration."""
        return (
            f"Config(env={self.env}, "
            f"base_path={self.base_path}, "
            f"data_root={self.data_root})"
        )
    
    def __repr__(self) -> str:
        return self.__str__()


# Global configuration instance
_config_instance: Optional[Config] = None


def get_config(reload: bool = False, **kwargs) -> Config:
    """
    Get global configuration instance (singleton pattern).
    
    Args:
        reload: Force reload of configuration
        **kwargs: Arguments passed to Config constructor
        
    Returns:
        Global Config instance
    """
    global _config_instance
    
    if _config_instance is None or reload:
        _config_instance = Config(**kwargs)
        
    return _config_instance


def reset_config():
    """Reset global configuration instance (mainly for testing)."""
    global _config_instance
    _config_instance = None


# Convenience functions
def get_data_path(level: str) -> Path:
    """Get data path for specified level."""
    return get_config().data_paths[level]


def get_sensor_path(level: str, sensor_type: str, depth: Optional[str] = None) -> Path:
    """Get sensor-specific data path."""
    return get_config().get_sensor_data_path(level, sensor_type, depth)


def get_maintenance_file(file_type: str) -> Path:
    """Get maintenance file path."""
    return get_config().maintenance_files[file_type]