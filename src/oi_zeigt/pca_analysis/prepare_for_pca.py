"""
Prepare FITS data for PCA analysis.

This module filters FITS data to include only specific sources and objects,
and fills telluric line regions with Gaussian noise.
"""

import logging
import os
from pathlib import Path
from typing import Optional

import numpy as np
import yaml
from astropy.io import fits

logger = logging.getLogger(__name__)


def load_mission_parameters(mission_id: str) -> dict:
    """
    Load mission-specific parameters from YAML file.
    
    Parameters
    ----------
    mission_id : str
        Mission identifier
    
    Returns
    -------
    dict
        Mission parameters including telluric_line_center and telluric_line_width
    """
    yaml_file = Path(__file__).parent / 'mission_id_parameters.yml'
    
    if not yaml_file.exists():
        logger.warning(f"Mission parameters file not found: {yaml_file}")
        return {}
    
    try:
        with open(yaml_file, 'r') as f:
            config = yaml.safe_load(f)
        
        # Find matching mission ID (may be a substring)
        for key in config:
            if key and mission_id in key:
                logger.debug(f"Found mission config for '{mission_id}': {key}")
                return config[key] or {}
        
        logger.warning(f"Mission '{mission_id}' not found in parameters file")
        return {}
    except Exception as e:
        logger.warning(f"Error loading mission parameters: {e}")
        return {}


def get_telluric_indices(mission_id: str, velocity_axis_kms: np.ndarray) -> Optional[np.ndarray]:
    """
    Get channel indices for telluric line regions based on mission parameters.
    
    Parameters
    ----------
    mission_id : str
        Mission identifier
    velocity_axis_kms : ndarray
        Velocity axis in km/s
    
    Returns
    -------
    ndarray or None
        Boolean array of channel indices to fill with noise
    """
    params = load_mission_parameters(mission_id)
    
    if not params or 'telluric_line_center' not in params:
        return None
    
    center_km_s = params['telluric_line_center']
    width_km_s = params.get('telluric_line_width', 30)
    
    # Create mask for velocities within the line region
    v_min = center_km_s - width_km_s / 2.0
    v_max = center_km_s + width_km_s / 2.0
    
    telluric_mask = (velocity_axis_kms >= v_min) & (velocity_axis_kms <= v_max)
    
    if np.any(telluric_mask):
        indices = np.where(telluric_mask)[0]
        logger.info(f"Telluric mask for {mission_id}: center={center_km_s} km/s, width={width_km_s} km/s")
        logger.info(f"  Channels {indices[0]}-{indices[-1]} to be filled with noise")
        return telluric_mask
    
    return None


def fill_telluric_with_noise(fits_file: str, output_fits: str, 
                            pca_source: str = "SKYCHOPDIFF",
                            object_filter: str = "M51CENTER") -> None:
    """
    Filter FITS file and fill telluric lines with Gaussian noise.
    
    Parameters
    ----------
    fits_file : str
        Input FITS file path
    output_fits : str
        Output FITS file path
    pca_source : str
        Source to filter for (default: SKYCHOPDIFF)
    object_filter : str
        Object substring to filter for (default: M51CENTER)
    """
    logger.info(f"Opening {fits_file}...")
    
    with fits.open(fits_file) as hdul:
        header = hdul[1].header
        data = hdul[1].data
        
        # Get velocity axis from data columns
        n_channels = len(data['SPECTRUM'][0])
        
        try:
            # Try to get VELOCITY_AXIS column first (velocity per channel)
            if 'VELOCITY_AXIS' in data.dtype.names:
                velocity_axis_ms = data['VELOCITY_AXIS'][0]  # m/s per channel
                velocity_axis_kms = velocity_axis_ms / 1000.0  # Convert to km/s
                logger.info(f"Velocity axis from VELOCITY_AXIS column: {len(velocity_axis_kms)} channels, {velocity_axis_kms[0]:.1f} to {velocity_axis_kms[-1]:.1f} km/s")
            # Fall back to constructing from VELOCITY and DELTAV
            elif ('VELOCITY' in data.dtype.names and 'DELTAV' in data.dtype.names):
                velo_ref = float(data['VELOCITY'][0])
                deltav = float(data['DELTAV'][0])
                
                # Get reference pixel from header
                crpix1_spec = float(header.get('CRPIX1', 1.0))
                
                # Create velocity axis
                channel_indices = np.arange(n_channels, dtype=np.float64)
                velocity_axis = velo_ref + (channel_indices - (crpix1_spec - 1.0)) * deltav
                velocity_axis_kms = velocity_axis / 1000.0  # Convert m/s to km/s
                logger.info(f"Velocity axis from VELOCITY/DELTAV: {len(velocity_axis_kms)} channels, {velocity_axis_kms[0]:.1f} to {velocity_axis_kms[-1]:.1f} km/s")
            else:
                raise KeyError("VELOCITY_AXIS, VELOCITY, or DELTAV columns missing")
        except (AttributeError, KeyError, TypeError, ValueError) as e:
            logger.warning(f"Could not construct velocity axis: {e}")
            velocity_axis_kms = np.arange(n_channels)
        
        # Filter data
        logger.info(f"Filtering to OBJECT = '{pca_source}' or '{object_filter}'...")
        
        # Get masks for filtering
        objects = np.array([s.strip() for s in data['OBJECT']])
        mission_ids = np.array([s.strip() for s in data['MISSION_ID']])
        
        # Keep spectra where OBJECT equals either pca_source or object_filter
        combined_mask = (objects == pca_source) | (objects == object_filter)
        
        n_original = len(data)
        n_filtered = np.sum(combined_mask)
        logger.info(f"Filtered from {n_original} to {n_filtered} spectra")
        
        # Filter all data
        filtered_data = data[combined_mask]
        filtered_mission_ids = mission_ids[combined_mask]
        
        # Fill telluric lines with Gaussian noise
        logger.info("Filling telluric lines with Gaussian noise...")
        
        n_filled = 0
        current_mission_id = None
        for idx, spectrum in enumerate(filtered_data['SPECTRUM']):
            mission_id = filtered_mission_ids[idx]
            
            # Log when we start processing a new mission_id
            if mission_id != current_mission_id:
                logger.info(f"  Processing {mission_id}")
                current_mission_id = mission_id
            
            # Get velocity axis for this specific spectrum
            if 'VELOCITY_AXIS' in filtered_data.dtype.names:
                spectrum_velocity_ms = filtered_data['VELOCITY_AXIS'][idx]
                spectrum_velocity_kms = spectrum_velocity_ms / 1000.0  # Convert m/s to km/s
            else:
                # Fallback to shared velocity axis if per-spectrum not available
                spectrum_velocity_kms = velocity_axis_kms
            
            telluric_mask = get_telluric_indices(mission_id, spectrum_velocity_kms)
            
            if telluric_mask is not None and np.any(telluric_mask):
                # Get noise level from surrounding channels, ignoring NaNs
                non_telluric = ~telluric_mask
                if np.any(non_telluric):
                    # Use nanstd to ignore NaN values
                    non_telluric_values = spectrum[non_telluric]
                    noise_level = np.nanstd(non_telluric_values)
                    
                    # If all non-telluric values are NaN, use the full spectrum
                    if np.isnan(noise_level):
                        noise_level = np.nanstd(spectrum)
                else:
                    noise_level = np.nanstd(spectrum)
                
                # Only fill if we have a valid noise level
                if not np.isnan(noise_level) and noise_level > 0:
                    # Fill with Gaussian noise
                    noise = np.random.normal(0, noise_level, np.sum(telluric_mask))
                    spectrum[telluric_mask] = noise
                    n_filled += 1
                else:
                    logger.warning(f"Could not determine noise level for spectrum {idx}, skipping telluric fill")
        
        logger.info(f"Filled telluric lines in {n_filled} spectra")
        
        # Create output FITS
        logger.info(f"Writing output to {output_fits}...")
        
        os.makedirs(os.path.dirname(output_fits) or '.', exist_ok=True)
        
        # Create new HDU with filtered data
        primary_hdu = fits.PrimaryHDU(header=hdul[0].header)
        table_hdu = fits.BinTableHDU(filtered_data, header=header)
        
        hdul_out = fits.HDUList([primary_hdu, table_hdu])
        hdul_out.writeto(output_fits, overwrite=True)
        
        logger.info(f"✓ Wrote {n_filtered} spectra to {output_fits}")


def prepare_for_pca(fits_file: Optional[str] = None, 
                   output_fits: Optional[str] = None,
                   config: Optional[str] = None,
                   pca_source: Optional[str] = None,
                   object_filter: Optional[str] = None) -> None:
    """
    Main entry point for prepare_for_pca functionality.
    
    Loads configuration from config.toml if provided, using defaults for missing parameters.
    
    Parameters
    ----------
    fits_file : str, optional
        Input FITS file (overrides config)
    output_fits : str, optional
        Output FITS file (overrides config)
    config : str, optional
        Path to config.toml file
    pca_source : str, optional
        PCA source to filter for (overrides config)
    object_filter : str, optional
        Object substring to filter for (overrides config)
    """
    # Configure logging to display output
    if not logger.handlers:
        handler = logging.StreamHandler()
        handler.setLevel(logging.INFO)
        formatter = logging.Formatter('%(message)s')
        handler.setFormatter(formatter)
        logger.addHandler(handler)
        logger.setLevel(logging.INFO)
    
    # Load configuration
    if config or (not fits_file and not output_fits):
        try:
            import tomllib
        except ModuleNotFoundError:
            import tomli as tomllib
        
        config_path = config or "config.toml"
        logger.info(f"Loading configuration from {config_path}")
        
        with open(config_path, 'rb') as f:
            cfg = tomllib.load(f)
        
        # Get input file
        if not fits_file:
            output_cfg = cfg.get('output', {})
            if 'reduced_fits' in output_cfg:
                fits_file = output_cfg['reduced_fits']
                logger.info(f"Using reduced_fits from config: {fits_file}")
            else:
                raise ValueError("Could not find input file. Specify --fits or ensure output.reduced_fits in config.toml")
        
        # Get output file
        if not output_fits:
            output_cfg = cfg.get('output', {})
            if 'prepared_for_pca' in output_cfg:
                output_fits = output_cfg['prepared_for_pca']
                logger.info(f"Using prepared_for_pca from config: {output_fits}")
            else:
                raise ValueError("Could not find output file. Specify --output or ensure output.prepared_for_pca in config.toml")
        
        # Get PCA source
        if not pca_source:
            pca_cfg = cfg.get('pca', {})
            pca_source = pca_cfg.get('pca_source', 'SKYCHOPDIFF')
        
        # Get object filter
        if not object_filter:
            parameters_cfg = cfg.get('parameters', {})
            object_filter = parameters_cfg.get('object', 'M51CENTER')
    
    # Validate required parameters
    if not fits_file:
        raise ValueError("Input FITS file not specified")
    if not output_fits:
        raise ValueError("Output FITS file not specified")
    if not pca_source:
        pca_source = 'SKYCHOPDIFF'
    if not object_filter:
        object_filter = 'M51CENTER'
    
    logger.info(f"PCA source filter: {pca_source}")
    logger.info(f"Object filter: {object_filter}")
    
    # Process the file
    fill_telluric_with_noise(fits_file, output_fits, pca_source, object_filter)
