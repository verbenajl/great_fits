"""
Configuration loader for PCA analysis.

Handles loading and parsing of TOML configuration files and YAML mission metadata.
Provides a unified configuration interface for decomposition and correction steps.
"""

import os
from pathlib import Path
from typing import Any, Dict, Optional, Union
import logging

try:
    import tomllib  # Built-in for Python 3.11+
except ImportError:
    try:
        import tomli as tomllib  # Fallback for Python < 3.11
    except ImportError:
        tomllib = None

try:
    import yaml
except ImportError:
    yaml = None

from .errors import ConfigurationError, MissingConfigurationError, InvalidConfigurationError

logger = logging.getLogger(__name__)


class ConfigLoader:
    """
    Load and manage configuration for PCA decomposition and correction.
    
    Loads configuration from TOML files and mission metadata from YAML files.
    Provides access to configuration parameters with type checking and defaults.
    
    Attributes:
        config_dict (dict): Main configuration parameters
        mission_params (dict): Mission-specific parameters (flight metadata)
    """
    
    def __init__(self, config_file: Optional[Union[str, Path]] = None,
                 mission_file: Optional[Union[str, Path]] = None):
        """
        Initialize configuration loader.
        
        Parameters
        ----------
        config_file : str or Path, optional
            Path to TOML configuration file
        mission_file : str or Path, optional
            Path to YAML mission parameters file
        """
        self.config_dict = {}
        self.mission_params = {}
        
        if config_file:
            self.load_config(config_file)
        
        if mission_file:
            self.load_missions(mission_file)
    
    def load_config(self, config_file: Union[str, Path]) -> None:
        """
        Load configuration from TOML file.
        
        Parameters
        ----------
        config_file : str or Path
            Path to TOML configuration file
            
        Raises
        ------
        ConfigurationError
            If file doesn't exist or TOML parsing fails
        """
        config_file = Path(config_file)
        
        if not config_file.exists():
            raise ConfigurationError(f"Configuration file not found: {config_file}")
        
        if tomllib is None:
            raise ConfigurationError("toml library not available. Install tomli: pip install tomli")

        try:
            with open(config_file, "rb") as f:
                self.config_dict = tomllib.load(f)
            logger.info(f"Loaded configuration from {config_file}")
        except Exception as e:
            raise ConfigurationError(f"Failed to parse configuration file: {e}")
        except Exception as e:
            raise ConfigurationError(f"Error loading configuration: {e}")
    
    def load_missions(self, mission_file: Union[str, Path]) -> None:
        """
        Load mission-specific parameters from YAML file.
        
        Parameters
        ----------
        mission_file : str or Path
            Path to YAML mission parameters file
            
        Raises
        ------
        ConfigurationError
            If file doesn't exist or YAML parsing fails
        """
        mission_file = Path(mission_file)
        
        if not mission_file.exists():
            raise ConfigurationError(f"Mission file not found: {mission_file}")
        
        if yaml is None:
            raise ConfigurationError("pyyaml library not available. Install with: pip install pyyaml")
        
        try:
            with open(mission_file, 'r') as f:
                self.mission_params = yaml.safe_load(f)
            logger.info(f"Loaded mission parameters from {mission_file}")
        except yaml.YAMLError as e:
            raise ConfigurationError(f"Failed to parse mission file: {e}")
        except Exception as e:
            raise ConfigurationError(f"Error loading missions: {e}")
    
    def get(self, key: str, default: Any = None, required: bool = False) -> Any:
        """
        Get configuration parameter with optional default value.
        
        Parameters
        ----------
        key : str
            Configuration key (supports dot notation for nested keys, e.g., 'pca.decompose.number_components')
        default : Any, optional
            Default value if key not found
        required : bool, optional
            If True, raise error if key is missing
            
        Returns
        -------
        Any
            Configuration value
            
        Raises
        ------
        MissingConfigurationError
            If key is required but not found
        """
        value = self._get_nested(self.config_dict, key)
        
        if value is None:
            if required:
                raise MissingConfigurationError(f"Required configuration parameter missing: {key}")
            return default
        
        return value
    
    def get_mission(self, mission_id: str, key: str = None, default: Any = None) -> Any:
        """
        Get mission-specific parameter.
        
        Parameters
        ----------
        mission_id : str
            Mission identifier (e.g., '2017-02-10_GR_F373')
        key : str, optional
            Specific parameter key within mission. If None, returns all mission params.
        default : Any, optional
            Default value if not found
            
        Returns
        -------
        Any
            Mission parameter value or dict of all parameters
        """
        if mission_id not in self.mission_params:
            return default
        
        mission_data = self.mission_params[mission_id]
        
        if key is None:
            return mission_data
        
        return mission_data.get(key, default)
    
    @staticmethod
    def _get_nested(d: dict, key: str) -> Any:
        """
        Get value from nested dictionary using dot notation.
        
        Parameters
        ----------
        d : dict
            Dictionary to search
        key : str
            Key with dot notation (e.g., 'pca.decompose.number_components')
            
        Returns
        -------
        Any
            Value from nested dict, or None if not found
        """
        keys = key.split('.')
        value = d
        
        for k in keys:
            if isinstance(value, dict):
                value = value.get(k)
            else:
                return None
        
        return value
    
    def get_line_parameters(self, mission_id: str) -> Dict[str, int]:
        """
        Get line masking parameters for a specific mission.
        
        Parameters
        ----------
        mission_id : str
            Mission identifier
            
        Returns
        -------
        dict
            Dictionary with 'telluric_line_center' and 'telluric_line_width'
            
        Raises
        ------
        ConfigurationError
            If line parameters are missing for the mission
        """
        mission_data = self.get_mission(mission_id)
        
        if mission_data is None:
            raise ConfigurationError(f"Mission not found in configuration: {mission_id}")
        
        try:
            center = mission_data['telluric_line_center']
            width = mission_data['telluric_line_width']
            return {'center': center, 'width': width}
        except KeyError as e:
            raise ConfigurationError(
                f"Missing line parameter {e} for mission {mission_id}"
            )
    
    def get_drop_filters(self, mission_id: str) -> Dict[str, Any]:
        """
        Get data quality drop filters for a mission.
        
        Parameters
        ----------
        mission_id : str
            Mission identifier
            
        Returns
        -------
        dict
            Drop filters (telescopes to exclude, bad scans, etc.)
            Empty dict if no filters defined
        """
        mission_data = self.get_mission(mission_id)
        
        if mission_data is None:
            return {}
        
        return mission_data.get('drop', {})
    
    def validate(self) -> bool:
        """
        Validate that essential configuration is present.
        
        Returns
        -------
        bool
            True if configuration is valid
            
        Raises
        ------
        MissingConfigurationError
            If essential parameters are missing
        """
        # Check for at least basic structure
        if not self.config_dict:
            raise MissingConfigurationError("Configuration is empty")
        
        logger.info("Configuration validation passed")
        return True
    
    def load_correction_config(self, mission_id: str):
        """
        Load correction configuration for Phase 3.
        
        Parameters
        ----------
        mission_id : str
            Mission identifier
            
        Returns
        -------
        CorrectionConfig
            Configuration object for correction module
            
        Raises
        ------
        ConfigurationError
            If configuration is invalid
        """
        # Import here to avoid circular imports
        from .correct import CorrectionConfig
        
        # Get correction parameters from config
        n_components = self.get('pca.decompose.number_components', default=5)
        science_source = self.get('pca.correct.science_source', default='M51CENTER')
        preserve_line = self.get('pca.correct.preserve_science_line', default=True)
        mask_width = self.get('pca.correct.mask_width', default=20)
        fit_window = self.get('pca.correct.fit_window', default=512)
        regularization = self.get('pca.correct.regularization', default=0.0)
        
        config = CorrectionConfig(
            n_components=n_components,
            science_source=science_source,
            preserve_science_line=preserve_line,
            mask_width=mask_width,
            fit_window=fit_window,
            regularization=regularization,
        )
        
        logger.info(f"Loaded correction config: {n_components} components, "
                   f"source={science_source}")
        
        return config
    
    def __repr__(self) -> str:
        """String representation of configuration."""
        return (
            f"ConfigLoader("
            f"config_keys={list(self.config_dict.keys())}, "
            f"missions={len(self.mission_params)} missions"
            f")"
        )


def create_default_config_template() -> str:
    """
    Create a template TOML configuration file.
    
    Returns
    -------
    str
        TOML template for PCA configuration
    """
    template = """
# PCA Decomposition Configuration
# ================================

[pca.common]
# Enable/disable PCA analysis
enabled = true

# Source spectra for PCA decomposition
pca_source = "SKYCHOPDIFF"

# Line window to ignore during analysis
# Format: [channel_start, channel_end]
line_window = [450, 550]

# Tag for identifying decomposition run
tag = "decomposition_v1"

# Frequency-space smoothing of components
smoothing_kernel_size = 3

# Window for rolling RMS noise estimate
rolling_noise_window = 11

# Type of decomposition: "PCA", "SparsePCA", or "ICA"
decomposition_type = "PCA"

# Normalize by std dev before PCA
do_scale = false

[pca.decompose]
# Number of principal components to derive
number_components = 5

# Add SKY-DIFF data in addition to SKYCHOPDIFF
add_sky_diff = false

# Z-score cutoff for outlier detection (false to disable)
noise_cutoff = false

# Scramble spectra (experimental)
scramble = false

[pca.correct]
# Min variance explained ratio to use component
cutoff = false

# Directory for diagnostic plots
output_folder = "pca_plots"

# Write component spectra to output
export_components = false

# Max noise ratio threshold for component selection (default: 15)
global_noise_ratio_cutoff = 15

# Override automatic noise ratio cutoff
noise_ratio_cutoff = false

# Gaussian blur kernel size for line detection
line_kernel_size = 61

# Write variance explained ratios to file
dump_exp_variance = false
"""
    return template
