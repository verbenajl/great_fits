"""
PCA-based spectral correction for FITS data.

This module applies PCA decomposition results to correct spectra in FITS files.
It follows the same workflow as pca_correct.py but adapted for FITS binary tables.

The correction workflow:
1. Load decomposed PCA components from pickle file
2. For each spectrum, fit PCA components using noise-ratio weighting
3. Subtract selected components based on noise thresholds
4. Generate diagnostic plots showing before/after correction

IMPORTANT ALIGNMENT NOTE:
- PCA components are eigenvectors stored in channel space (0-699 for 700 channels)
- Decomposition masks out telluric lines during component generation (km/s -> channels)
- When applying correction, spectra from the SAME FITS file have the SAME velocity axis
- This ensures automatic alignment by velocity - no manual conversion needed
- Components and spectra both use identical channel indexing from the FITS file
"""

import os
import sys
import pickle
import logging
import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure
import numpy.ma as ma
from astropy.io import fits
from pathlib import Path
import warnings
import yaml

from .config import ConfigLoader

try:
    import cv2 as cv
except ImportError:
    cv = None

logger = logging.getLogger(__name__)

# Configure logging if not already configured
if not logging.getLogger().hasHandlers():
    logging.basicConfig(
        level=logging.INFO,
        format='%(levelname)s:%(name)s:%(message)s'
    )


def find_science_lines(spectra, kernel_size=51, cutoff_std=2.0, smoothing_kernel=None, 
                       velocity_axis_kms=None, velocity_window_kms=None):
    """
    Detect science lines in spectra using 2D OpenCV image processing.
    
    Treats the 2D array of spectra as a 2D image and applies Gaussian blur
    and threshold detection. This is similar to pca_utilities.find_lines().
    
    If velocity_window_kms is provided, detection is restricted to channels
    within that velocity range, masking out other regions before processing.
    
    Parameters
    ----------
    spectra : ndarray
        2D array of spectra (n_spectra, n_channels)
    kernel_size : int
        Kernel size for Gaussian blur (must be odd, default: 51)
    cutoff_std : float
        Threshold in units of standard deviation (default: 2.0)
    smoothing_kernel : int, optional
        Smoothing kernel size (not used, for API compatibility)
    velocity_axis_kms : ndarray, optional
        Velocity axis in km/s for each channel (shape: n_channels)
    velocity_window_kms : tuple, optional
        (v_min, v_max) velocity range in km/s to constrain detection
    
    Returns
    -------
    ndarray or None
        2D boolean mask (n_spectra, n_channels) indicating detected lines
        Returns None if no lines detected in any spectrum
    """
    if cv is None:
        logger.warning("OpenCV not available, skipping line detection")
        return None
    
    spectra = np.asarray(spectra)
    if len(spectra.shape) == 1:
        # Duplicate 10 times so the 2D Gaussian blur behaves the same as in
        # the legacy find_lines() call in pca_correct.py (lines 1238-1242).
        spectra = np.vstack([spectra] * 10)
    
    n_spectra, n_channels = spectra.shape
    logger.debug(f"Line detection: input shape {n_spectra}x{n_channels}, kernel_size={kernel_size}, cutoff_std={cutoff_std}")
    
    # Ensure kernel_size is positive odd integer
    if kernel_size <= 0 or kernel_size % 2 == 0:
        kernel_size = max(1, kernel_size - 1)
    if kernel_size < 3:
        kernel_size = 3
    
    logger.debug(f"  Using kernel_size={kernel_size}")
    
    # Determine velocity window mask if provided
    velocity_mask = np.ones(n_channels, dtype=bool)
    if velocity_axis_kms is not None and velocity_window_kms is not None:
        v_min, v_max = velocity_window_kms
        velocity_mask = (velocity_axis_kms >= v_min) & (velocity_axis_kms <= v_max)
        n_channels_in_window = np.sum(velocity_mask)
        logger.debug(f"  Velocity window: {v_min}-{v_max} km/s, {n_channels_in_window}/{n_channels} channels in window")
    
    # Use FULL spectra for detection (don't mask before OpenCV processing)
    spectra_for_detection = spectra
    
    # CRITICAL: Reject spectra with zero or negligible signal
    # Calculate spectrum amplitude for each row
    spectrum_amplitudes = np.max(np.abs(spectra_for_detection), axis=1)
    global_max_amplitude = np.max(spectrum_amplitudes)
    
    # Find spectra with negligible signal (< 1e-10 of the maximum)
    amplitude_threshold = 1e-10 * global_max_amplitude
    if amplitude_threshold == 0:
        amplitude_threshold = 1e-10  # Fallback if all spectra are zero
    
    good_spectrum_mask = spectrum_amplitudes > amplitude_threshold
    logger.debug(f"  Amplitude threshold: {amplitude_threshold:.2e}")
    logger.debug(f"  Good spectra: {np.sum(good_spectrum_mask)}/{n_spectra}")
    
    # Normalize spectra to 8-bit range for OpenCV
    data_min = np.min(spectra_for_detection)
    data_max = np.max(spectra_for_detection)
    
    logger.debug(f"  Data range: [{data_min:.6f}, {data_max:.6f}]")
    
    if data_max > data_min:
        data_8uc = 255 * (spectra_for_detection - data_min) / (data_max - data_min)
    else:
        data_8uc = np.zeros_like(spectra_for_detection)
    
    img = data_8uc.astype(np.uint8)
    img = cv.normalize(img, None, 0, 100, cv.NORM_MINMAX)
    
    logger.debug(f"  Normalized image range: [{np.min(img)}, {np.max(img)}]")
    
    # Apply Gaussian blur
    gray = cv.GaussianBlur(img, (kernel_size, kernel_size), 0)
    
    # Calculate initial threshold
    valid_pixels = gray[np.where(gray != 0)]
    if len(valid_pixels) > 0:
        mean = valid_pixels.mean()
        std = valid_pixels.std()
    else:
        mean = 0
        std = 1
    
    threshold_val = mean + cutoff_std * std
    logger.debug(f"  Initial: mean={mean:.2f}, std={std:.2f}, threshold={threshold_val:.2f}")
    
    _, threshold = cv.threshold(gray, threshold_val, 1, cv.THRESH_BINARY)
    
    initial_detections = np.sum(threshold > 0)
    logger.debug(f"  Initial detections: {initial_detections} pixels")
    
    # Iterative refinement: recalculate threshold excluding detected regions.
    # new_gray is created once and accumulates zeroed-out regions across iterations,
    # matching the legacy find_lines() behaviour in pca_utilities.py.
    new_mean = 0
    new_std = 0
    new_gray = gray.copy()

    while (not np.isclose(new_mean, mean).all()) and (not np.isclose(new_std, std).all()):
        new_mean = mean
        new_std = std

        # Accumulate detected regions into new_gray (do not reset each iteration)
        new_gray[np.where(threshold == 1)] = 0

        # Recalculate statistics on remaining data
        valid_pixels = new_gray[np.where(new_gray != 0)]
        if len(valid_pixels) > 0:
            mean = valid_pixels.mean()
            std = valid_pixels.std()
        else:
            mean = 0
            std = 1

        # Recalculate threshold against the original blurred image
        threshold_val = mean + cutoff_std * std
        _, threshold = cv.threshold(gray, threshold_val, 1, cv.THRESH_BINARY)
    
    final_detections = np.sum(threshold > 0)
    logger.debug(f"  Final detections: {final_detections} pixels")
    
    # Convert to boolean mask
    mask = threshold.astype(bool)
    
    # CRITICAL: Zero out detections in spectra with negligible signal
    # These spectra should not have any lines detected
    mask[~good_spectrum_mask, :] = False
    logger.debug(f"  Zeroed out {np.sum(~good_spectrum_mask)} spectra with negligible signal")
    
    # Apply velocity window constraint: zero out channels outside the window (AFTER detection)
    if not np.all(velocity_mask):
        mask[:, ~velocity_mask] = False
        logger.debug(f"  Applying velocity window filter to detection results")
    
    # Log per-spectrum and per-channel statistics
    n_spectra_with_lines = np.sum(np.any(mask, axis=1))
    line_channels = np.any(mask, axis=0)
    n_channels_with_lines = np.sum(line_channels)
    
    logger.debug(f"  Result: {n_spectra_with_lines}/{n_spectra} spectra with lines, {n_channels_with_lines}/{n_channels} channels detected")
    
    if n_channels_with_lines > 0:
        line_indices = np.where(line_channels)[0]
        logger.debug(f"    Channel range: {line_indices[0]}-{line_indices[-1]}")
        if velocity_axis_kms is not None:
            v_min_detected = velocity_axis_kms[line_indices[0]]
            v_max_detected = velocity_axis_kms[line_indices[-1]]
            logger.debug(f"    Velocity range: {v_min_detected:.2f}-{v_max_detected:.2f} km/s")
    
    # Return None if no lines detected in any spectrum (to match original behavior)
    return mask if n_spectra_with_lines > 0 else None


def get_velocity_window_channels(velocity_axis, window_km_s, tolerance=1.0):
    """
    Convert velocity window to channel indices.
    
    Parameters
    ----------
    velocity_axis : ndarray
        Velocity axis in km/s
    window_km_s : tuple or list
        (v_min, v_max) in km/s
    tolerance : float
        Tolerance for matching (default: 1.0 km/s)
    
    Returns
    -------
    tuple
        (channel_min, channel_max) indices
    """
    if velocity_axis is None or window_km_s is None:
        return None
    
    v_min, v_max = window_km_s
    
    # Find closest channel indices
    ch_min = np.argmin(np.abs(velocity_axis - v_min))
    ch_max = np.argmin(np.abs(velocity_axis - v_max))
    
    return (min(ch_min, ch_max), max(ch_min, ch_max))


def load_mission_parameters(mission_id):
    """
    Load mission-specific parameters from YAML file.
    
    Parameters
    ----------
    mission_id : str
        Mission identifier (e.g., '2017-02-01_GR_F367')
    
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


def get_telluric_line_mask(mission_id, velocity_axis_kms, n_spectra):
    """
    Create a mask for telluric lines based on mission-specific parameters.
    
    Parameters
    ----------
    mission_id : str
        Mission identifier
    velocity_axis_kms : ndarray
        Velocity axis in km/s (already converted from m/s)
    n_spectra : int
        Number of spectra
    
    Returns
    -------
    ndarray or None
        Boolean mask array of shape (n_spectra, len(velocity_axis_kms)) where True = line region
    """
    params = load_mission_parameters(mission_id)
    
    if not params or 'telluric_line_center' not in params:
        return None
    
    center_km_s = params['telluric_line_center']
    width_km_s = params.get('telluric_line_width', 30)
    
    # velocity_axis_kms is already in km/s
    # Create mask for velocities within the line region
    v_min = center_km_s - width_km_s / 2.0
    v_max = center_km_s + width_km_s / 2.0
    
    telluric_mask = (velocity_axis_kms >= v_min) & (velocity_axis_kms <= v_max)
    
    # Expand to all spectra
    mask = np.zeros((n_spectra, len(velocity_axis_kms)), dtype=bool)
    mask[:, telluric_mask] = True
    
    logger.info(f"Telluric line mask: center={center_km_s} km/s, width={width_km_s} km/s")
    logger.info(f"  Channels {np.where(telluric_mask)[0][0]}-{np.where(telluric_mask)[0][-1]} masked")
    
    return mask


class PCACorrector:
    """Apply PCA-based correction to spectra in FITS files."""
    
    def __init__(self, decomposition_pkl, config=None, 
                 line_kernel_size=51, line_cutoff_std=2.0, 
                 smoothing_kernel_size=None, line_window_velocities=None):
        """
        Initialize PCA corrector with decomposition results.
        
        Parameters
        ----------
        decomposition_pkl : str
            Path to pickle file containing PCA decomposition results
        config : dict, optional
            Configuration parameters
        line_kernel_size : int
            Kernel size for line detection (default: 51, must be odd)
        line_cutoff_std : float
            Standard deviation cutoff for line detection (default: 2.0)
        smoothing_kernel_size : int, optional
            Smoothing kernel size for line detection refinement
        line_window_velocities : tuple, optional
            (v_min, v_max) velocity range in km/s to constrain line detection
        """
        self.config = config or {}
        self.line_kernel_size = line_kernel_size
        self.line_cutoff_std = line_cutoff_std
        self.smoothing_kernel_size = smoothing_kernel_size
        self.line_window_velocities = line_window_velocities
        
        self.decomposition = self._load_decomposition(decomposition_pkl)
        self.pca = self.decomposition.pca if hasattr(self.decomposition, 'pca') else None
        
        # Handle both dict format (from decomposition.py) and object format
        if self.pca is None:
            if isinstance(self.decomposition, dict):
                # Dict format: extract components directly
                if 'components' in self.decomposition and 'explained_variance_ratio' in self.decomposition:
                    self.components = self.decomposition['components']
                    self.explained_variance_ratio = self.decomposition['explained_variance_ratio']
                else:
                    raise ValueError("Decomposition dict missing 'components' or 'explained_variance_ratio'")
            else:
                raise ValueError("No PCA object found in decomposition")
        else:
            # Object format: extract from sklearn PCA object
            self.components = self.pca.components_
            self.explained_variance_ratio = self.pca.explained_variance_ratio_
        
        logger.info(f"✓ Loaded PCA with {len(self.components)} components")
        logger.info(f"  Variance explained: {np.sum(self.explained_variance_ratio)*100:.2f}%")
        logger.info(f"  Line detection kernel size: {line_kernel_size}")
        logger.info(f"  Line detection cutoff: {line_cutoff_std} sigma")
        if line_window_velocities:
            logger.info(f"  Line detection window: {line_window_velocities[0]}-{line_window_velocities[1]} km/s")
    
    def _load_decomposition(self, pkl_path):
        """Load decomposition from pickle file."""
        if not os.path.exists(pkl_path):
            raise FileNotFoundError(f"Decomposition file not found: {pkl_path}")
        
        file_size = os.path.getsize(pkl_path)
        logger.info(f"DEBUG: Loading decomposition from {pkl_path} (size: {file_size} bytes)")
        
        with open(pkl_path, 'rb') as f:
            decomposition = pickle.load(f)
        
        # Log the keys in the loaded decomposition dict for debugging
        if isinstance(decomposition, dict):
            logger.info(f"DEBUG: Decomposition dict keys: {list(decomposition.keys())}")
            if 'reference_spectrum_indices' in decomposition:
                indices = decomposition.get('reference_spectrum_indices')
                if indices is not None:
                    logger.info(f"DEBUG: reference_spectrum_indices found, shape={indices.shape}, dtype={indices.dtype}")
                else:
                    logger.info(f"DEBUG: reference_spectrum_indices key exists but value is None")
            else:
                logger.info(f"DEBUG: reference_spectrum_indices key NOT in dict")
                # Also check file size to see if it's a new or old pickle
                logger.info(f"DEBUG: Pickle file size: {file_size} bytes")
        
        return decomposition
    
    def load_reference_spectra_from_fits(self, fits_file, spectrum_indices, hdu_index=1, spectrum_col='SPECTRUM'):
        """
        Load reference (SKYCHOPDIFF) spectra from FITS file using stored indices.
        
        Parameters
        ----------
        fits_file : str
            Path to FITS file
        spectrum_indices : np.ndarray
            Array of indices pointing to SKYCHOPDIFF spectra in the FITS file
        hdu_index : int
            HDU index (default: 1)
        spectrum_col : str
            Name of spectrum column (default: 'SPECTRUM')
        
        Returns
        -------
        np.ndarray or None
            Reference spectra array, or None if indices not available
        """
        if spectrum_indices is None or len(spectrum_indices) == 0:
            logger.warning("spectrum_indices is None or empty")
            return None
        
        try:
            logger.info(f"Loading {len(spectrum_indices)} reference SKYCHOPDIFF spectra from {fits_file}")
            with fits.open(fits_file) as hdul:
                hdu = hdul[hdu_index]
                data = hdu.data
                
                # Load spectra at the specified indices
                reference_spectra = data[spectrum_col][spectrum_indices]
                logger.info(f"✓ Successfully loaded {len(reference_spectra)} reference SKYCHOPDIFF spectra")
                return np.array(reference_spectra)
        except Exception as e:
            logger.warning(f"Failed to load reference spectra: {e}")
            return None
    
    @staticmethod
    def load_decompositions_from_dir(decomp_dir):
        """
        Load all decomposition files from a directory.
        
        Looks for files matching pattern: decomposition_*_components.pkl
        Returns dict mapping (mission_id, telescop) -> decomposition
        Falls back to mapping mission_id -> decomposition for single-mission files
        
        Parameters
        ----------
        decomp_dir : str
            Directory containing decomposition files
        
        Returns
        -------
        dict
            Dict mapping (mission_id, telescop) or mission_id -> PCACorrector instance
        """
        import glob
        
        decomp_dir = Path(decomp_dir)
        if not decomp_dir.exists():
            logger.warning(f"Decomposition directory not found: {decomp_dir}")
            return {}
        
        decompositions = {}
        decomp_files = sorted(glob.glob(str(decomp_dir / "decomposition_*_components.pkl")))
        
        logger.info(f"Found {len(decomp_files)} decomposition files in {decomp_dir}")
        
        if len(decomp_files) == 0:
            logger.warning(f"No decomposition files found matching pattern: decomposition_*_components.pkl")
            return {}
        
        for decomp_file in decomp_files:
            try:
                # Extract mission_id and telescope from filename
                # Format: decomposition_{YYYY-MM-DD}_{GR/HR}_{F###}_{LFAH/LFAV}_{PX##}_{S}_components.pkl
                # OR: decomposition_{YYYYMMDD}_{mission_id}_components.pkl
                # Example: decomposition_2016-05-12_GR_F296_LFAH_PX00_S_components.pkl
                filename = os.path.basename(decomp_file)
                base = filename.replace('decomposition_', '').replace('_components.pkl', '')
                
                # The date is at the START and contains dashes: YYYY-MM-DD
                # After the date comes: GR/HR_F###_LFAH/LFAV_PX##_S
                flight_date = None
                mission_id = None
                telescope = None
                
                # Try to match YYYY-MM-DD at the beginning
                import re
                match = re.match(r'^(\d{4})-(\d{2})-(\d{2})_(.+)$', base)
                if match:
                    year, month, day = match.group(1), match.group(2), match.group(3)
                    flight_date = f"{year}{month}{day}"
                    date_str = f"{year}-{month}-{day}"  # Keep with dashes for mission_id
                    remaining = match.group(4)
                    
                    # Now parse: GR_F296_LFAH_PX00_S
                    parts = remaining.split('_')
                    
                    # Look for telescope pattern: LFAH_PX##_S or LFAV_PX##_S
                    # The mission_id is everything before the telescope code, INCLUDING the date
                    for i in range(len(parts) - 1, -1, -1):
                        if parts[i] in ['LFAH', 'LFAI', 'LFAV', 'LFBI', 'LFBH']:
                            # Found telescope type
                            if i + 2 < len(parts) and parts[i+1].startswith('PX') and parts[i+2] == 'S':
                                # Reconstruct mission_id with the date at the front
                                mission_id = date_str + '_' + '_'.join(parts[:i])
                                telescope = '_'.join(parts[i:i+3])  # LFAH_PX00_S
                                break
                
                # Fallback for old YYYYMMDD format
                if flight_date is None:
                    parts = base.split('_')
                    if len(parts) >= 2 and len(parts[-1]) >= 7 and parts[-1][:7].isdigit():
                        flight_date = parts[-1]
                        mission_id = '_'.join(parts[:-1])
                
                if flight_date is None:
                    # No valid date found, skip
                    logger.warning(f"  Could not parse date from {filename}")
                    continue
                
                # Load the decomposition
                with open(decomp_file, 'rb') as f:
                    decomposition = pickle.load(f)
                
                # Create a PCACorrector-like object
                corrector = type('DecompositionData', (), {})()
                corrector.components = decomposition.get('components', np.array([]))
                corrector.explained_variance_ratio = decomposition.get('explained_variance_ratio', np.array([]))
                corrector.mean_spectrum = decomposition.get('mean_spectrum', None)
                corrector.metadata = decomposition.get('metadata', {})
                corrector.reference_spectrum_indices = decomposition.get('reference_spectrum_indices', None)
                
                # Store with appropriate key
                if mission_id and telescope:
                    key = (mission_id, telescope)
                    decompositions[key] = corrector
                    logger.info(f"  Loaded {mission_id}/{telescope}: {len(corrector.components)} components")
                elif mission_id:
                    # Fall back to old format (mission_id only)
                    decompositions[mission_id] = corrector
                    logger.info(f"  Loaded {mission_id}: {len(corrector.components)} components")
                else:
                    # No mission_id could be parsed, skip
                    logger.warning(f"  Could not parse mission_id from {filename}")
                    continue
                
            except Exception as e:
                logger.warning(f"  Failed to load {decomp_file}: {e}")
        
        return decompositions
    
    def get_decomposition_for_mission(self, mission_id, telescope, mission_decompositions):
        """
        Get the right decomposition for a given mission and telescope.
        
        Tries exact match (mission_id, telescope) first, then tries base telescope name,
        then falls back to mission_id only.
        
        Parameters
        ----------
        mission_id : str
            Mission identifier (e.g., "2016-05-12_GR_F296")
        telescope : str
            Telescope identifier (e.g., LFAH_PX00_S)
        mission_decompositions : dict
            Dict mapping (mission_id, telescope) or mission_id -> decomposition
        
        Returns
        -------
        PCACorrector or dict
            Decomposition for this mission/telescope, or None if not found
        """
        # Try exact match first: (mission_id, telescope)
        key_tuple = (mission_id, telescope)
        if key_tuple in mission_decompositions:
            return mission_decompositions[key_tuple]
        
        # Try with base telescope name (extract first 4 chars: LFAH, LFAV, etc.)
        if len(telescope) > 4 and telescope[:4] in ['LFAH', 'LFAI', 'LFAV', 'LFBI', 'LFBH', 'PRISM']:
            base_telescope = telescope[:4]
            key_tuple_base = (mission_id, base_telescope)
            if key_tuple_base in mission_decompositions:
                logger.debug(f"Found decomposition for {mission_id}/{base_telescope} (from {telescope})")
                return mission_decompositions[key_tuple_base]
        
        # Fall back to mission_id only
        if mission_id in mission_decompositions:
            return mission_decompositions[mission_id]
        
        logger.warning(f"No decomposition found for {mission_id}/{telescope}, using default")
        return None
    
    def get_noise_ratio(self, spectrum, component, scaled_spectrum=None, smoothing_kernel_size=None):
        """
        Calculate noise ratio: noise in spectrum / noise in component.
        
        This metric measures the reliability of a component for this spectrum:
        - noise_ratio = std(spectrum) / std(scaled_component)
        - Higher values = component explains less signal, more noise-like
        - Lower values = component explains signal well
        
        Optionally applies box-car smoothing to distinguish between noise and signal,
        following the approach in pca_correct.py lines 939-956.
        
        Parameters
        ----------
        spectrum : ndarray
            Original spectrum
        component : ndarray
            PCA component
        scaled_spectrum : ndarray, optional
            Scaled spectrum (if different from original)
        smoothing_kernel_size : int, optional
            Kernel size for box-car smoothing. If None, no smoothing applied.
        
        Returns
        -------
        float
            Noise ratio (higher = less reliable component)
        """
        if scaled_spectrum is None:
            scaled_spectrum = spectrum
        
        # Calculate scaled component
        coeff = np.dot(component, scaled_spectrum)
        scaled_comp = coeff * component
        
        # Calculate spectrum std (unsmoothed)
        spectrum_std = np.nanstd(scaled_spectrum)
        component_std = np.nanstd(scaled_comp)
        
        # Apply smoothing if kernel size provided (like original pca_correct.py)
        # CRITICAL: Handle NaNs properly - convolve() does NOT ignore NaNs!
        if smoothing_kernel_size and smoothing_kernel_size > 0:
            # Replace NaNs with interpolated values before convolution
            nan_mask = np.isnan(scaled_spectrum)
            if np.any(nan_mask):
                # Interpolate NaNs before smoothing
                valid_idx = np.where(~nan_mask)[0]
                if len(valid_idx) > 0:
                    spectrum_for_smooth = scaled_spectrum.copy()
                    spectrum_for_smooth[nan_mask] = np.interp(
                        np.where(nan_mask)[0], 
                        valid_idx, 
                        scaled_spectrum[valid_idx]
                    )
                else:
                    # All NaN - can't smooth
                    spectrum_for_smooth = scaled_spectrum
            else:
                spectrum_for_smooth = scaled_spectrum
            
            kernel = np.ones(smoothing_kernel_size) / smoothing_kernel_size
            spectrum_smoothed = np.convolve(spectrum_for_smooth, kernel, mode="same")
            spectrum_std = np.nanstd(spectrum_smoothed)
        
        # Handle edge cases that could produce NaN or inf
        if np.isnan(spectrum_std) or np.isnan(component_std):
            # If either std is NaN, we have bad/masked data (all NaN)
            return np.inf
        
        # Prevent division by zero with a minimum threshold
        epsilon = 1e-15
        if component_std < epsilon:
            # Component has negligible amplitude
            return 1e10
        
        if spectrum_std < epsilon:
            # Spectrum is essentially flat/zero
            return 1e5
        
        noise_ratio = spectrum_std / component_std
        
        # Ensure the result is finite
        if not np.isfinite(noise_ratio):
            return 1e10
        
        return noise_ratio
    
    def fit_coefficients(self, spectrum, good_channels=None):
        """
        Fit PCA coefficients to a spectrum.
        
        Parameters
        ----------
        spectrum : ndarray
            Input spectrum
        good_channels : ndarray, optional
            Boolean mask of good channels (False = bad)
        
        Returns
        -------
        ndarray
            Coefficients for each component
        """
        if good_channels is None:
            good_channels = np.ones(len(spectrum), dtype=bool)
        
        # Fit coefficients using only good channels
        # Slice both components and spectrum to good channels only
        if not np.all(good_channels):
            # Only use good channels for fitting
            coeff = np.dot(self.components[:, good_channels], spectrum[good_channels])
        else:
            # All channels are good
            coeff = np.dot(self.components, spectrum)
        
        return coeff
    
    def apply_correction(self, spectrum, good_channels=None, good_channels_for_subtraction=None,
                         cutoff_variance=None, cutoff_noise_ratio=None, verbose=False, smoothing_kernel_size=None):
        """
        Apply PCA correction to a single spectrum.
        
        NOTE ON VELOCITY ALIGNMENT:
        The components are eigenvectors in channel space (0-699). The spectra being
        corrected come from the same FITS file, so they use the SAME velocity axis.
        This means channel indices automatically correspond to the same velocities
        in both components and spectra - no velocity-to-channel conversion needed here.
        
        The telluric line masking (if any) was already applied during decomposition,
        where mission-specific km/s parameters were converted to channel indices.
        
        Parameters
        ----------
        spectrum : ndarray
            Input spectrum to correct
        good_channels : ndarray, optional
            Boolean mask of good channels for FITTING (should include science lines and telluric)
        good_channels_for_subtraction : ndarray, optional
            Boolean mask of channels where components should be SUBTRACTED
            (excludes science lines and telluric). If None, uses good_channels.
        cutoff_variance : float, optional
            Only use components explaining > this fraction of variance
        cutoff_noise_ratio : float, optional
            Skip components where noise_ratio > this threshold
        verbose : bool
            Print debug information
        smoothing_kernel_size : int, optional
            Kernel size for box-car smoothing in noise_ratio calculation.
            If None, no smoothing applied.
        
        Returns
        -------
        ndarray
            Corrected spectrum
        dict
            Correction details (coefficients, used components, etc.)
        """
        if good_channels is None:
            good_channels = np.ones(len(spectrum), dtype=bool)
        
        # If not specified, use the same channels for subtraction as for fitting
        if good_channels_for_subtraction is None:
            good_channels_for_subtraction = good_channels.copy()
        
        # Fit coefficients
        coeff = self.fit_coefficients(spectrum, good_channels)
        
        if verbose:
            logger.info(f"  Fitted coefficients: {coeff}")
            logger.info(f"  Cutoff variance: {cutoff_variance}, Cutoff noise ratio: {cutoff_noise_ratio}")
            logger.info(f"  Explained variance ratios: {self.explained_variance_ratio}")
        
        corrected = spectrum.copy()
        used_components = []
        component_info = {}
        
        # Subtract components based on criteria
        for i, (comp, var_ratio) in enumerate(zip(self.components, 
                                                   self.explained_variance_ratio)):
            # Check variance cutoff
            if cutoff_variance is not None and var_ratio < cutoff_variance:
                if verbose:
                    logger.info(f"  Skip comp {i}: variance {var_ratio:.4f} < {cutoff_variance}")
                component_info[i] = {
                    'used': False,
                    'reason': 'variance_cutoff',
                    'coeff': coeff[i],
                    'variance': var_ratio
                }
                continue
            
            # Calculate noise ratio (with optional smoothing)
            scaled_comp = coeff[i] * comp
            noise_ratio = self.get_noise_ratio(spectrum, comp, smoothing_kernel_size=smoothing_kernel_size)
            
            if verbose:
                logger.info(f"  Comp {i}: coeff={coeff[i]:.6f}, var_ratio={var_ratio:.6f}, "
                           f"noise_ratio={noise_ratio:.6f}, scaled_comp_range=[{scaled_comp.min():.6f}, {scaled_comp.max():.6f}]")
            
            # Skip if noise ratio is NaN (bad data, masked/corrupted spectrum)
            # Note: Very large noise ratios (1e5+) are now handled in cutoff_noise_ratio check below
            if np.isnan(noise_ratio):
                if verbose:
                    logger.info(f"    -> SKIPPED (noise_ratio is NaN, bad/masked data)")
                component_info[i] = {
                    'used': False,
                    'reason': 'bad_noise_ratio',
                    'coeff': coeff[i],
                    'noise_ratio': noise_ratio
                }
                continue
            
            if cutoff_noise_ratio is not None and noise_ratio > cutoff_noise_ratio:
                if verbose:
                    logger.info(f"    -> SKIPPED (noise_ratio {noise_ratio:.4f} > {cutoff_noise_ratio})")
                component_info[i] = {
                    'used': False,
                    'reason': 'noise_ratio_cutoff',
                    'coeff': coeff[i],
                    'noise_ratio': noise_ratio
                }
                continue
            
            # Subtract component from subtraction channels only (excludes science lines, telluric)
            if verbose:
                logger.info(f"    -> APPLYING: subtracting from {np.sum(good_channels_for_subtraction)} channels")
            corrected[good_channels_for_subtraction] -= scaled_comp[good_channels_for_subtraction]
            used_components.append(i)
            
            if verbose:
                logger.info(f"      After subtraction: corrected range=[{corrected.min():.6f}, {corrected.max():.6f}]")
            
            component_info[i] = {
                'used': True,
                'coeff': coeff[i],
                'variance': var_ratio,
                'noise_ratio': noise_ratio
            }
        
        details = {
            'coefficients': coeff,
            'n_components_used': len(used_components),
            'used_components': used_components,
            'component_info': component_info,
            'n_components_used': len(used_components)
        }
        
        return corrected, details
    
    def correct_fits_file(self, input_fits, output_fits, cutoff_variance=None,
                         cutoff_noise_ratio=None, hdu_index=1, spectrum_col='SPECTRUM',
                         object_filter=None, overwrite=False, generate_plots=False,
                         output_dir='output/pca_corrected', config_window=None,
                         detect_science_lines=True, scan_filter=None, subscan_filter=None,
                         telescope_filter=None, mission_id_filter=None, mission_decompositions=None):
        """
        Correct all spectra in a FITS file.
        
        Parameters
        ----------
        input_fits : str
            Input FITS file path
        output_fits : str
            Output FITS file path (corrected spectra)
        cutoff_variance : float, optional
            Only use components explaining > this fraction of variance
        cutoff_noise_ratio : float, optional
            Skip components where noise_ratio > this threshold
        hdu_index : int
            HDU index containing spectrum table (default: 1)
        spectrum_col : str
            Name of spectrum column (default: 'SPECTRUM')
        object_filter : str, optional
            Only correct spectra matching this OBJECT value
        overwrite : bool
            Overwrite output file if exists
        generate_plots : bool
            Generate diagnostic plots per scan/telescope
        output_dir : str
            Output directory for plots
        config_window : tuple, optional
            (v_min, v_max) in km/s - fallback window if line detection fails
        detect_science_lines : bool
            Detect science lines using OpenCV (default: True)
        scan_filter : int, optional
            Only correct spectra from this SCAN (for testing)
        subscan_filter : int, optional
            Only correct spectra from this SUBSCAN (for testing)
        telescope_filter : str, optional
            Only correct spectra from this TELESCOP (for testing)
        mission_id_filter : str, optional
            Only correct spectra from this MISSION_ID (for testing)
        mission_decompositions : dict, optional
            Dictionary mapping (mission_id, telescope) or mission_id -> decomposition
            If provided, uses per-mission/telescope decompositions instead of self.components
        
        Returns
        -------
        dict
            Correction statistics
        """
        logger.info(f"Opening {input_fits}...")
        with fits.open(input_fits) as hdul:
            hdu = hdul[hdu_index]
            data = hdu.data.copy()
            header = hdu.header.copy()
            
            # Get mission_id from the data table (it's a column, not a header keyword)
            # Use the first spectrum's mission_id
            try:
                mission_id = data['MISSION_ID'][0].strip()
            except (KeyError, IndexError):
                mission_id = "UNKNOWN"
            
            n_spectra = len(data)
            logger.info(f"✓ Loaded {n_spectra} spectra from {spectrum_col}")
            
            # Filter by object if requested
            if object_filter:
                mask = np.array([object_filter.strip() in s.strip() for s in data['OBJECT']])
                indices = np.where(mask)[0]
                logger.info(f"  Filtering to {len(indices)} spectra containing '{object_filter}'")
            else:
                indices = np.arange(n_spectra)
            
            # Additional filtering for scan/subscan/telescope (for testing)
            if scan_filter is not None:
                mask = data['SCAN'][indices] == scan_filter
                indices = indices[mask]
                logger.info(f"  Filtering to {len(indices)} spectra with SCAN={scan_filter}")
            
            if subscan_filter is not None:
                mask = data['SUBSCAN'][indices] == subscan_filter
                indices = indices[mask]
                logger.info(f"  Filtering to {len(indices)} spectra with SUBSCAN={subscan_filter}")
            
            if telescope_filter is not None:
                def to_string(x):
                    if isinstance(x, bytes):
                        return x.decode('utf-8').strip()
                    return str(x).strip()
                
                telescopes = np.array([to_string(data['TELESCOP'][i]) for i in indices])
                mask = telescopes == telescope_filter.strip()
                indices = indices[mask]
                logger.info(f"  Filtering to {len(indices)} spectra with TELESCOP='{telescope_filter}'")
            
            if mission_id_filter is not None:
                def to_string(x):
                    if isinstance(x, bytes):
                        return x.decode('utf-8').strip()
                    return str(x).strip()
                
                mission_ids_filtered = np.array([to_string(data['MISSION_ID'][i]) for i in indices])
                mask = mission_ids_filtered == mission_id_filter.strip()
                indices = indices[mask]
                logger.info(f"  Filtering to {len(indices)} spectra with MISSION_ID='{mission_id_filter}'")
                
                # Debug: show what mission_ids are in the filtered data
                if len(indices) > 0:
                    filtered_mission_ids = np.unique(mission_ids_filtered[mask])
                    logger.info(f"  Actual mission_ids in filtered data: {filtered_mission_ids}")
                    if mission_decompositions:
                        logger.info(f"  Available decompositions: {list(mission_decompositions.keys())[:3]}... ({len(mission_decompositions)} total)")
                else:
                    logger.warning(f"  WARNING: No spectra found with MISSION_ID='{mission_id_filter}'")
                    # Show what mission_ids ARE available
                    unique_mission_ids = np.unique(mission_ids_filtered)
                    logger.info(f"  Available mission_ids in data: {unique_mission_ids}")
            
            # CRITICAL: Filter out spectra with all zeros or negligible signal
            # These spectra cannot be corrected anyway and contaminate line detection
            logger.info("Filtering out spectra with negligible signal...")
            n_before_zero_filter = len(indices)
            
            # Calculate amplitude (max absolute value) for each spectrum, ignoring NaN
            spectrum_amplitudes = np.array([
                np.nanmax(np.abs(data[spectrum_col][i]))
                for i in indices
            ])
            
            # Find global maximum amplitude
            valid_amplitudes = spectrum_amplitudes[~np.isnan(spectrum_amplitudes)]
            if len(valid_amplitudes) > 0:
                global_max_amplitude = np.max(valid_amplitudes)
            else:
                global_max_amplitude = 1.0
            
            # Filter: keep only spectra with amplitude > 1e-10 * global_max
            # This filters out spectra that are all zeros or all NaN
            amplitude_threshold = 1e-10 * global_max_amplitude
            good_amplitude_mask = spectrum_amplitudes > amplitude_threshold
            
            indices = indices[good_amplitude_mask]
            n_after_zero_filter = len(indices)
            
            if n_after_zero_filter < n_before_zero_filter:
                logger.info(f"  Removed {n_before_zero_filter - n_after_zero_filter} spectra with negligible signal ({n_after_zero_filter} remaining)")
            else:
                logger.info(f"  All {n_after_zero_filter} spectra have sufficient signal")
            
            # Get velocity axis if available
            velocity_axis = None
            if 'VELOCITY_AXIS' in hdu.columns.names:
                velocity_axis = data['VELOCITY_AXIS'][0]  # Same for all spectra
                velocity_axis_kms = velocity_axis / 1000.0  # Convert m/s to km/s
                logger.info(f"  Using VELOCITY_AXIS: {len(velocity_axis)} channels")
            else:
                velocity_axis_kms = None
            
            # Prepare for iterative line detection (following original pca_correct.py)
            # Convert string columns to consistent format once
            def to_string(x):
                if isinstance(x, bytes):
                    return x.decode('utf-8').strip()
                return str(x).strip()
            
            mission_ids = np.array([to_string(data['MISSION_ID'][i]) for i in indices])
            telescopes = np.array([to_string(data['TELESCOP'][i]) for i in indices])
            scans = np.array([to_string(x) for x in data['SCAN'][indices]])
            subscans = np.array([to_string(x) for x in data['SUBSCAN'][indices]])
            
            # Get unique mission/telescope/scan/subscan combinations (needed for iterative line detection)
            unique_groups = np.unique(
                np.column_stack((mission_ids, telescopes, scans, subscans)),
                axis=0
            )
            
            # Also get unique mission/telescope/scan combinations for plotting
            unique_groups_for_plotting = np.unique(
                np.column_stack((mission_ids, telescopes, scans)),
                axis=0
            )
            
            logger.info(f"  Found {len(unique_groups)} unique mission/telescope/scan/subscan combinations (for line detection)")
            logger.info(f"  Will generate {len(unique_groups_for_plotting)} plots (one per mission/telescope/scan)")
            logger.info(f"  Total spectra to process: {len(indices)}")
            
            # Load telluric line mask from mission parameters
            telluric_line_mask = None
            if velocity_axis_kms is not None and mission_id != "UNKNOWN":
                try:
                    telluric_line_mask = get_telluric_line_mask(mission_id, velocity_axis_kms, len(indices))
                    if telluric_line_mask is not None:
                        telluric_channels = np.any(telluric_line_mask, axis=0)
                        n_telluric_channels = np.sum(telluric_channels)
                        logger.info(f"  Loaded telluric line mask: {n_telluric_channels} channels masked")
                except Exception as e:
                    logger.debug(f"  Could not load telluric mask: {e}")
            
            # Initialize corrected spectra and plots arrays
            corrected_spectra = data[spectrum_col].copy()
            original_spectra_for_plotting = data[spectrum_col].copy()
            correction_details = {}
            stats = {
                'total': len(indices),
                'corrected': 0,
                'failed': 0,
                'components_used': [],
                'lines_detected': False
            }
            
            # ============================================================================
            # ITERATIVE LINE DETECTION FOLLOWING ORIGINAL pca_correct.py APPROACH
            # ============================================================================
            # Step 1: First pass correction WITHOUT line mask to get preliminary corrected spectra
            logger.info("STEP 1: First pass correction (without line detection mask)")
            
            prelim_corrected_spectra = corrected_spectra.copy()
            detected_lines_mask = None
            
            # Log if we're using per-mission decompositions
            if mission_decompositions:
                logger.info(f"Using per-mission decompositions ({len(mission_decompositions)} available)")
                logger.info(f"  Sample keys: {list(mission_decompositions.keys())[:3]}")
            else:
                logger.info("Using single default decomposition for all spectra")
            
            for spec_idx, idx in enumerate(indices):
                try:
                    spectrum = data[spectrum_col][idx]
                    
                    # Select correct decomposition for this spectrum if using per-mission decompositions
                    if mission_decompositions:
                        spec_mission_id = data['MISSION_ID'][idx].strip()
                        spec_telescope = data['TELESCOP'][idx].strip()
                        decomp = self.get_decomposition_for_mission(spec_mission_id, spec_telescope, mission_decompositions)
                        if decomp is None:
                            if spec_idx < 3:  # Log only first 3 misses
                                logger.warning(f"Spectrum {idx}: No decomposition found for {spec_mission_id}/{spec_telescope} - keeping original (uncorrected)")
                            # Keep original spectrum uncorrected instead of skipping
                            prelim_corrected_spectra[idx] = spectrum
                            correction_details[idx] = {'n_components_used': 0, 'skipped': True, 'reason': 'no_decomposition'}
                            stats['failed'] += 1
                            continue
                        saved_components = self.components
                        saved_variance_ratio = self.explained_variance_ratio
                        self.components = decomp.components
                        self.explained_variance_ratio = decomp.explained_variance_ratio
                    
                    # Check for bad channels only (no line mask yet)
                    bad_channels = np.isnan(spectrum) | np.isinf(spectrum) | (spectrum == 0)
                    good_channels_for_fitting = ~bad_channels
                    
                    if not np.any(good_channels_for_fitting):
                        correction_details[idx] = {'n_components_used': 0, 'skipped': True}
                        if mission_decompositions:
                            self.components = saved_components
                            self.explained_variance_ratio = saved_variance_ratio
                        continue
                    
                    # First pass: correct without line mask
                    corrected, details = self.apply_correction(
                        spectrum,
                        good_channels=good_channels_for_fitting,
                        good_channels_for_subtraction=good_channels_for_fitting,
                        cutoff_variance=cutoff_variance,
                        cutoff_noise_ratio=cutoff_noise_ratio,
                        verbose=False,
                        smoothing_kernel_size=self.smoothing_kernel_size
                    )
                    
                    prelim_corrected_spectra[idx] = corrected
                    correction_details[idx] = details
                    
                    if mission_decompositions:
                        self.components = saved_components
                        self.explained_variance_ratio = saved_variance_ratio
                
                except Exception as e:
                    logger.debug(f"First pass correction failed for spectrum {idx}: {e}")
                    correction_details[idx] = {'n_components_used': 0, 'skipped': True}
            
            # Log coefficient statistics after first pass
            all_coefficients = []
            for idx in indices:
                if idx in correction_details and 'coefficients' in correction_details[idx]:
                    all_coefficients.append(correction_details[idx]['coefficients'])
            
            if all_coefficients:
                all_coefficients = np.array(all_coefficients)
                logger.info(f"COEFFICIENT STATISTICS (first pass, {len(all_coefficients)} spectra):")
                for comp_idx in range(all_coefficients.shape[1]):
                    comp_coeff = all_coefficients[:, comp_idx]
                    logger.info(f"  Component {comp_idx}: "
                              f"min={np.min(comp_coeff):.4e}, "
                              f"max={np.max(comp_coeff):.4e}, "
                              f"mean={np.mean(comp_coeff):.4e}, "
                              f"std={np.std(comp_coeff):.4e}, "
                              f"median={np.median(comp_coeff):.4e}")
            
            # Step 2: Detect lines FROM the corrected spectra (following original pca_correct.py)
            logger.info("STEP 2: Detecting lines from corrected spectra (iterative refinement)")
            
            detected_lines_mask = np.zeros((len(indices), len(velocity_axis_kms)), dtype=bool)
            group_thresholds = {}
            
            # Detect lines from each group's corrected spectra (3 iterations like original)
            for iteration in range(3):
                logger.info(f"  Iteration {iteration + 1}/3 of line detection")
                
                for mission_id_group, telescope_group, scan_id, subscan_id in unique_groups:
                    # Find spectra matching this mission/telescope/scan/subscan combination
                    group_mask = (
                        (mission_ids == mission_id_group) &
                        (telescopes == telescope_group) &
                        (scans == scan_id) &
                        (subscans == subscan_id)
                    )
                    group_indices = np.where(group_mask)[0]
                    if len(group_indices) == 0:
                        continue
                    
                    # Get CORRECTED spectra from first pass for line detection
                    group_corrected_spectra = prelim_corrected_spectra[indices[group_indices]]
                    
                    try:
                        # Detect lines in the corrected spectra
                        line_window_kms = None
                        if hasattr(self, 'line_window_velocities') and self.line_window_velocities:
                            line_window_kms = tuple(self.line_window_velocities)
                        
                        # Match legacy behaviour: raise cutoff_std to 3 on the
                        # final iteration (pca_correct.py line 1252-1253).
                        this_cutoff_std = 3 if iteration == 2 else self.line_cutoff_std
                        group_threshold = find_science_lines(
                            group_corrected_spectra,
                            kernel_size=self.line_kernel_size,
                            cutoff_std=this_cutoff_std,
                            smoothing_kernel=self.smoothing_kernel_size,
                            velocity_axis_kms=velocity_axis_kms,
                            velocity_window_kms=line_window_kms
                        )
                        
                        if group_threshold is not None:
                            # Save threshold for visualization (only on last iteration)
                            if iteration == 2:
                                group_key = (mission_id_group, telescope_group, scan_id, subscan_id)
                                group_thresholds[group_key] = group_threshold
                            
                            # Apply detected lines to mask
                            for group_spectrum_position, global_spectrum_idx in enumerate(group_indices):
                                line_mask = group_threshold[group_spectrum_position]
                                if np.any(line_mask):
                                    local_idx = np.where(indices == global_spectrum_idx)[0][0]
                                    detected_lines_mask[local_idx] = line_mask
                    
                    except Exception as e:
                        logger.debug(f"Line detection failed for {mission_id_group}/{telescope_group}/scan{scan_id}/subscan{subscan_id}: {e}")
                
                # If last iteration, skip re-correction
                if iteration < 2:
                    # Re-correct with the current detected lines mask
                    for spec_idx, idx in enumerate(indices):
                        try:
                            spectrum = data[spectrum_col][idx]
                            
                            if mission_decompositions:
                                spec_mission_id = data['MISSION_ID'][idx].strip()
                                spec_telescope = data['TELESCOP'][idx].strip()
                                decomp = self.get_decomposition_for_mission(spec_mission_id, spec_telescope, mission_decompositions)
                                if decomp is None:
                                    continue
                                saved_components = self.components
                                saved_variance_ratio = self.explained_variance_ratio
                                self.components = decomp.components
                                self.explained_variance_ratio = decomp.explained_variance_ratio
                            
                            bad_channels = np.isnan(spectrum) | np.isinf(spectrum) | (spectrum == 0)
                            good_channels_for_fitting = ~bad_channels
                            good_channels_for_subtraction = good_channels_for_fitting.copy()
                            
                            # Exclude detected lines from subtraction (but use for fitting - like original)
                            line_regions = detected_lines_mask[spec_idx]
                            good_channels_for_subtraction = good_channels_for_subtraction & ~line_regions
                            
                            if not np.any(good_channels_for_fitting):
                                if mission_decompositions:
                                    self.components = saved_components
                                    self.explained_variance_ratio = saved_variance_ratio
                                continue
                            
                            # Re-correct with detected lines excluded from subtraction
                            corrected, details = self.apply_correction(
                                spectrum,
                                good_channels=good_channels_for_fitting,
                                good_channels_for_subtraction=good_channels_for_subtraction,
                                cutoff_variance=cutoff_variance,
                                cutoff_noise_ratio=cutoff_noise_ratio,
                                verbose=False,
                                smoothing_kernel_size=self.smoothing_kernel_size
                            )
                            
                            prelim_corrected_spectra[idx] = corrected
                            
                            if mission_decompositions:
                                self.components = saved_components
                                self.explained_variance_ratio = saved_variance_ratio
                        
                        except Exception as e:
                            logger.debug(f"Re-correction iteration {iteration + 1} failed for spectrum {idx}: {e}")
            
            # Step 3: Final correction with detected lines excluded from fitting
            logger.info("STEP 3: Final correction with detected lines excluded from fitting")
            
            for spec_idx, idx in enumerate(indices):
                try:
                    spectrum = data[spectrum_col][idx]
                    
                    # Select correct decomposition for this spectrum if using per-mission decompositions
                    if mission_decompositions:
                        spec_mission_id = data['MISSION_ID'][idx].strip()
                        spec_telescope = data['TELESCOP'][idx].strip()
                        decomp = self.get_decomposition_for_mission(spec_mission_id, spec_telescope, mission_decompositions)
                        if decomp is None:
                            logger.debug(f"Spectrum {idx}: no decomposition for {spec_mission_id}/{spec_telescope}, keeping original (uncorrected)")
                            corrected_spectra[idx] = spectrum
                            stats['failed'] += 1
                            continue
                        # Temporarily swap components for this spectrum
                        saved_components = self.components
                        saved_variance_ratio = self.explained_variance_ratio
                        self.components = decomp.components
                        self.explained_variance_ratio = decomp.explained_variance_ratio
                    
                    # Check for bad channels
                    bad_channels = np.isnan(spectrum) | np.isinf(spectrum) | (spectrum == 0)
                    good_channels_for_fitting = ~bad_channels
                    good_channels_for_subtraction = good_channels_for_fitting.copy()
                    
                    n_good = np.sum(good_channels_for_fitting)
                    n_line = 0
                    
                    # Exclude detected lines from FITTING (following original pca_correct.py)
                    if detected_lines_mask is not None and np.any(detected_lines_mask):
                        line_regions = detected_lines_mask[spec_idx]
                        good_channels_for_fitting = good_channels_for_fitting & ~line_regions
                        good_channels_for_subtraction = good_channels_for_subtraction & ~line_regions
                        n_line = np.sum(line_regions)
                        stats['lines_detected'] = True
                    
                    # Also exclude telluric lines
                    if telluric_line_mask is not None:
                        telluric_regions = telluric_line_mask[min(spec_idx, len(telluric_line_mask)-1)]
                        good_channels_for_fitting = good_channels_for_fitting & ~telluric_regions
                        good_channels_for_subtraction = good_channels_for_subtraction & ~telluric_regions

                    
                    if not np.any(good_channels_for_fitting):
                        logger.warning(f"Spectrum {idx}: no good channels for fitting (had {n_good} good, {n_line if detected_lines_mask is not None else 0} lines), skipping")
                        corrected_spectra[idx] = spectrum
                        # Record as skipped but still track it
                        correction_details[idx] = {'n_components_used': 0, 'skipped': True}
                        stats['failed'] += 1
                        # Restore original decomposition if using per-mission
                        if mission_decompositions:
                            self.components = saved_components
                            self.explained_variance_ratio = saved_variance_ratio
                        continue
                    
                    # Final correction with detected lines excluded from fitting
                    corrected, details = self.apply_correction(
                        spectrum, 
                        good_channels=good_channels_for_fitting,
                        good_channels_for_subtraction=good_channels_for_subtraction,
                        cutoff_variance=cutoff_variance,
                        cutoff_noise_ratio=cutoff_noise_ratio,
                        verbose=(spec_idx < 3),  # Log first 3 spectra
                        smoothing_kernel_size=self.smoothing_kernel_size
                    )
                    
                    corrected_spectra[idx] = corrected
                    correction_details[idx] = details
                    
                    stats['corrected'] += 1
                    if details['n_components_used'] > 0:
                        stats['components_used'].append(details['n_components_used'])
                    
                    # Log first few spectra to track what's happening
                    if spec_idx < 3:
                        logger.info(f"Spectrum {idx}: {details['n_components_used']} components used, "
                                  f"correction_avg={np.mean(np.abs(corrected - spectrum)):.4e}")
                    
                    # Restore original decomposition if using per-mission
                    if mission_decompositions:
                        self.components = saved_components
                        self.explained_variance_ratio = saved_variance_ratio
                
                except Exception as e:
                    logger.error(f"Spectrum {idx}: ERROR: {str(e)}", exc_info=True)
                    corrected_spectra[idx] = data[spectrum_col][idx]
                    stats['failed'] += 1
                    # Restore original decomposition if using per-mission
                    if mission_decompositions:
                        self.components = saved_components
                        self.explained_variance_ratio = saved_variance_ratio
            
            # Update data with corrected spectra
            data[spectrum_col] = corrected_spectra
            
            # Preserve calibration objects (TSYS, TAU_SIG) even when filtering
            # These are needed for downstream analysis and pairing with M51CENTER
            if object_filter:
                # Get indices of calibration objects to preserve
                calibration_mask = np.array([
                    s.strip() in ('TSYS', 'TAU_SIG') 
                    for s in data['OBJECT']
                ])
                calibration_indices = np.where(calibration_mask)[0]
                
                # Combine filtered science object indices with calibration indices
                all_indices_to_keep = np.concatenate([indices, calibration_indices])
                all_indices_to_keep = np.sort(all_indices_to_keep)
                
                logger.info(f"Filters applied - preserving {len(indices)} corrected {object_filter} + {len(calibration_indices)} calibration spectra (TSYS/TAU_SIG)")
                filtered_data = data[all_indices_to_keep]
                
                # CRITICAL: Recalculate TSYS_INDEX and TAU_SIG_INDEX for the output file
                # These indices must point to rows in the OUTPUT file, not the input file
                # Create a mapping from old row numbers (in input) to new row numbers (in output)
                row_mapping = {old_idx: new_idx for new_idx, old_idx in enumerate(all_indices_to_keep)}
                
                # Recalculate indices
                if 'TSYS_INDEX' in filtered_data.columns.names:
                    old_tsys_indices = filtered_data['TSYS_INDEX'].copy()
                    new_tsys_indices = np.full_like(old_tsys_indices, -1, dtype=np.int32)
                    
                    for new_row_idx in range(len(filtered_data)):
                        old_row_idx = all_indices_to_keep[new_row_idx]
                        old_tsys_ref = old_tsys_indices[new_row_idx]
                        
                        # If this row had a valid TSYS_INDEX in the input, remap it to output
                        if old_tsys_ref >= 0 and old_tsys_ref in row_mapping:
                            new_tsys_indices[new_row_idx] = row_mapping[old_tsys_ref]
                    
                    filtered_data['TSYS_INDEX'] = new_tsys_indices
                    n_remapped = np.sum(new_tsys_indices >= 0)
                    logger.debug(f"Remapped TSYS_INDEX: {n_remapped} entries now point to correct rows in output file")
                
                if 'TAU_SIG_INDEX' in filtered_data.columns.names:
                    old_tau_indices = filtered_data['TAU_SIG_INDEX'].copy()
                    new_tau_indices = np.full_like(old_tau_indices, -1, dtype=np.int32)
                    
                    for new_row_idx in range(len(filtered_data)):
                        old_row_idx = all_indices_to_keep[new_row_idx]
                        old_tau_ref = old_tau_indices[new_row_idx]
                        
                        # If this row had a valid TAU_SIG_INDEX in the input, remap it to output
                        if old_tau_ref >= 0 and old_tau_ref in row_mapping:
                            new_tau_indices[new_row_idx] = row_mapping[old_tau_ref]
                    
                    filtered_data['TAU_SIG_INDEX'] = new_tau_indices
                    n_remapped = np.sum(new_tau_indices >= 0)
                    logger.debug(f"Remapped TAU_SIG_INDEX: {n_remapped} entries now point to correct rows in output file")
                
                # Preserve ALL header keywords before table restructuring
                header_keywords_to_preserve = dict(header)
                
                output_table = fits.BinTableHDU(filtered_data, header=header)
                
                # Restore any header keywords that were lost during table restructuring
                for key, value in header_keywords_to_preserve.items():
                    if key not in output_table.header:
                        logger.debug(f"Restoring header keyword {key} that was lost during table restructuring")
                        output_table.header[key] = value
                        
            elif scan_filter is not None or subscan_filter is not None or telescope_filter or mission_id_filter:
                logger.info(f"Filters applied - writing only {len(indices)} filtered spectra to output")
                filtered_data = data[indices]
                
                # Preserve ALL header keywords before table restructuring
                header_keywords_to_preserve = dict(header)
                
                output_table = fits.BinTableHDU(filtered_data, header=header)
                
                # Restore any header keywords that were lost during table restructuring
                for key, value in header_keywords_to_preserve.items():
                    if key not in output_table.header:
                        logger.debug(f"Restoring header keyword {key} that was lost during table restructuring")
                        output_table.header[key] = value
                        
            else:
                logger.info(f"No filters applied - writing all {len(data)} spectra to output")
                
                # Preserve ALL header keywords before table restructuring
                header_keywords_to_preserve = dict(header)
                
                output_table = fits.BinTableHDU(data, header=header)
                
                # Restore any header keywords that were lost during table restructuring
                for key, value in header_keywords_to_preserve.items():
                    if key not in output_table.header:
                        logger.debug(f"Restoring header keyword {key} that was lost during table restructuring")
                        output_table.header[key] = value
            
            # Write output FITS
            os.makedirs(os.path.dirname(output_fits) or '.', exist_ok=True)
            hdul_out = fits.HDUList([hdul[0], output_table])
            hdul_out.writeto(output_fits, overwrite=overwrite)
            
            logger.info(f"✓ Saved corrected spectra to {output_fits}")
            logger.info(f"  Total: {stats['total']} spectra processed")
            logger.info(f"  Corrected: {stats['corrected']} spectra")
            if stats['failed'] > 0:
                logger.info(f"  Kept uncorrected: {stats['failed']} spectra (no decomposition or other issues)")
            if stats['components_used']:
                logger.info(f"  Components used: {np.mean(stats['components_used']):.1f} ± "
                           f"{np.std(stats['components_used']):.1f}")
            
            # Generate plots if requested
            # Use only the corrected spectra for plotting (filtered indices)
            if generate_plots:
                # IMPORTANT: Use all filtered spectra for plotting, including those that may not have been corrected
                # (e.g., due to bad channels). The plot function can show all data and note which was corrected.
                # We DON'T filter again here because the filtering already happened via correction_details tracking.
                
                # Create a filtered data table with spectra from filtered indices
                filtered_data = data[indices]
                filtered_original = original_spectra_for_plotting[indices]
                filtered_corrected = corrected_spectra[indices]
                filtered_details = correction_details
                
                # Note: detected_lines_mask and telluric_line_mask are already sized for the filtered spectra
                # (they were initialized as np.zeros((len(indices), n_channels))), so we pass them directly
                # without further filtering
                
                # Generate plots with all filtered spectra
                self._generate_plots(filtered_data, filtered_original, filtered_corrected, filtered_details,
                                   velocity_axis, mission_id, output_dir,
                                   science_line_mask=detected_lines_mask,
                                   telluric_line_mask=telluric_line_mask,
                                   input_fits=input_fits,
                                   group_thresholds=group_thresholds,
                                   global_indices=indices,
                                   mission_decompositions=mission_decompositions)
            
            return stats
    
    def _generate_plots(self, data, original_spectra, corrected_spectra, correction_details,
                       velocity_axis, mission_id, output_dir, science_line_mask=None, telluric_line_mask=None,
                       input_fits=None, group_thresholds=None, global_indices=None, mission_decompositions=None):
        """
        Generate diagnostic plots matching pca_correct.py structure.
        
        Layout (5 columns):
        COL1: Mean original/corrected, then components
        COL2: (top) spectra used to derive components (waterfall)
              (bottom) spectra before correction (heatmap)
        COL3: (top) difference original vs corrected
              (bottom) spectra after correction
        COL4: (top) masked array used for fitting
              (bottom) coefficient values heatmap
        COL5: Example spectrum correction with component breakdown
        """
        os.makedirs(output_dir, exist_ok=True)
        logger.info(f"Generating diagnostic plots...")
        
        # Load SKYCHOPDIFF reference spectra from the same input FITS file for waterfall plots
        skychopdiff_spectra = {}
        if mission_decompositions and input_fits:
            try:
                from astropy.io import fits as pyfits
                with pyfits.open(input_fits) as hdul:
                    # Load all SKYCHOPDIFF spectra and their metadata
                    hdu = hdul[1]
                    all_data = hdu.data
                    all_object = np.array([s.strip() for s in all_data['OBJECT']])
                    all_mission = np.array([s.strip() for s in all_data['MISSION_ID']])
                    all_telescope = np.array([s.strip() for s in all_data['TELESCOP']])
                    
                    # Filter to only SKYCHOPDIFF spectra
                    sky_mask = all_object == 'SKYCHOPDIFF'
                    sky_spectra = all_data['SPECTRUM'][sky_mask]
                    sky_mission = all_mission[sky_mask]
                    sky_telescope = all_telescope[sky_mask]
                    
                    # Group by mission/telescope
                    for key in mission_decompositions.keys():
                        # key is (mission_id, telescope)
                        mission_str, telescope_str = key
                        combo_mask = (sky_mission == mission_str) & (sky_telescope == telescope_str)
                        if np.any(combo_mask):
                            skychopdiff_spectra[key] = sky_spectra[combo_mask]
                            logger.info(f"Loaded {len(sky_spectra[combo_mask])} SKYCHOPDIFF spectra for {mission_str}/{telescope_str}")
            except Exception as e:
                logger.warning(f"Failed to load SKYCHOPDIFF spectra for plotting: {e}")
        
        # Debug logging for masks
        logger.info(f"_generate_plots called with science_line_mask: {science_line_mask is not None}, telluric_line_mask: {telluric_line_mask is not None}")
        if science_line_mask is not None:
            logger.info(f"  science_line_mask shape: {science_line_mask.shape}, any detections: {np.any(science_line_mask)}")
        if telluric_line_mask is not None:
            logger.info(f"  telluric_line_mask shape: {telluric_line_mask.shape}")
        
        # Get unique scan/subscan/telescope/mission combinations
        scans = data['SCAN']
        subscans = data['SUBSCAN']
        telescopes = data['TELESCOP']
        missions = np.array([data['MISSION_ID'][i].strip() for i in range(len(data))])
        
        # Group by (mission_id, telescope, scan) like original pca_correct.py
        # This creates one plot per (mission_id, telescope, scan) group
        # ALL subscans within that group are ACCUMULATED and shown together
        unique_mission_tel_scan = sorted(list(set(zip(missions, telescopes, scans))))
        
        logger.info(f"_generate_plots received {len(data)} spectra")
        logger.info(f"Grouping by (mission_id, telescope, scan): {len(unique_mission_tel_scan)} groups")
        logger.info(f"correction_details contains {len(correction_details)} entries")
        logger.info(f"group_thresholds passed to _generate_plots: {len(group_thresholds) if group_thresholds else 0} groups")
        if group_thresholds:
            logger.info(f"  First few threshold keys: {list(group_thresholds.keys())[:3]}")
        
        total_plots = 0
        for mission_id, telescope, scan in unique_mission_tel_scan:
            # Get all spectra for this (mission_id, telescope, scan) combination
            # This includes ALL subscans for this scan
            group_mask = (missions == mission_id) & (telescopes == telescope) & (scans == scan)
            group_indices_in_filtered = np.where(group_mask)[0]
            
            if len(group_indices_in_filtered) == 0:
                continue
            
            # Convert filtered indices back to global indices
            if global_indices is not None:
                group_indices_global = global_indices[group_indices_in_filtered]
            else:
                group_indices_global = group_indices_in_filtered
            
            # Get all spectra for this group (all subscans accumulated)
            original_spectra_subset = original_spectra[group_indices_in_filtered]
            corrected_spectra_subset = corrected_spectra[group_indices_in_filtered]
            
            # Safety check: ensure we have spectra
            if len(original_spectra_subset) == 0 or len(original_spectra_subset[0]) == 0:
                logger.warning(f"mission_id={mission_id}, telescope={telescope}, scan={scan}: No spectral data (skipping plot)")
                continue
            
            # Get indices of spectra with correction details (were actually corrected)
            # Now check global indices against correction_details
            usable_global_indices = [global_idx for global_idx in group_indices_global
                                    if global_idx in correction_details]
            
            if len(usable_global_indices) == 0:
                logger.info(f"mission_id={mission_id}, telescope={telescope}, scan={scan}: No corrected spectra in group (skipping plot)")
                continue
            
            # CRITICAL: Filter subset arrays to ONLY include successfully corrected spectra
            # Create a mapping from global indices back to positions in group_indices_global
            usable_indices_in_group = []
            for i, global_idx in enumerate(group_indices_global):
                if global_idx in correction_details:
                    usable_indices_in_group.append(i)
            
            # Now filter the subset arrays using these indices (which are positions in group_indices_in_filtered)
            original_spectra_subset = original_spectra_subset[usable_indices_in_group]
            corrected_spectra_subset = corrected_spectra_subset[usable_indices_in_group]
            
            # Additional safety check: ensure spectra have data
            if original_spectra_subset.ndim < 2 or original_spectra_subset.shape[0] == 0 or original_spectra_subset.shape[1] == 0:
                logger.warning(f"mission_id={mission_id}, telescope={telescope}, scan={scan}: Invalid spectrum shape {original_spectra_subset.shape} (skipping plot)")
                continue
            
            logger.info(f"Creating plot {total_plots + 1} for mission_id={mission_id}, telescope={telescope}, scan={scan}: {len(usable_indices_in_group)} usable spectra out of {len(group_indices_in_filtered)}")
            total_plots += 1
            
            # Setup x_axis
            if velocity_axis is not None and len(velocity_axis) > 0:
                x_axis = velocity_axis / 1000.0  # Convert m/s to km/s
                x_label = 'Velocity (km/s)'
            else:
                # Safe to access now since we've checked shape above
                n_channels = original_spectra_subset.shape[1]
                x_axis = np.arange(n_channels)
                x_label = 'Channel'
            
            # Safety check: ensure x_axis has elements
            if len(x_axis) == 0:
                logger.warning(f"mission_id={mission_id}, telescope={telescope}, scan={scan}: Empty x_axis (skipping plot)")
                continue
            
            # Create figure (25x10 to match pca_correct.py)
            fig = Figure(figsize=(25, 10))
            canvas = FigureCanvasAgg(fig)
            
            plot_width = 0.17
            padding = 0.02
            
            # ==================================================================
            # COLUMN 1: ORIGINAL/CORRECTED MEAN + COMPONENTS (LEFT SIDE)
            # ==================================================================
            original_mean = np.mean(original_spectra_subset, axis=0)
            corrected_mean = np.mean(corrected_spectra_subset, axis=0)
            
            # Handle NaN/Inf in mean spectra
            original_mean = np.nan_to_num(original_mean, nan=0.0, posinf=0.0, neginf=0.0)
            corrected_mean = np.nan_to_num(corrected_mean, nan=0.0, posinf=0.0, neginf=0.0)
            
            min_val = min(np.min(original_mean), np.min(corrected_mean)) * 1.1
            max_val = max(np.max(original_mean), np.max(corrected_mean)) * 1.1
            
            # Skip this plot if axis limits are invalid
            if not (np.isfinite(min_val) and np.isfinite(max_val)):
                logger.warning(f"Skipping plot: invalid axis limits for mission_id={mission_id}, telescope={telescope}, scan={scan}")
                continue
            
            # Get the correct decomposition for this specific mission/telescope
            plot_decomposition = self.get_decomposition_for_mission(mission_id, telescope, mission_decompositions)
            
            # Use components from the correct mission decomposition, or fall back to self.components
            if plot_decomposition:
                components_to_plot = plot_decomposition.components
                variance_ratio_to_plot = plot_decomposition.explained_variance_ratio
            else:
                components_to_plot = self.components
                variance_ratio_to_plot = self.explained_variance_ratio
            
            n_components_show = len(components_to_plot)
            spectral_height = 0.85 / (2 + n_components_show)  # 2 for orig/corr + n for components
            
            # Original mean
            ax1 = fig.add_axes([padding, 0.90 - spectral_height, plot_width, spectral_height])
            ax1.plot(x_axis, original_mean, 'r-', lw=1.5)
            ax1.set_title('Original Mean', fontsize=10)
            ax1.set_ylim(min_val, max_val)
            ax1.grid(True, alpha=0.3)
            ax1.tick_params(labelsize=8)
            
            # Corrected mean
            ax2 = fig.add_axes([padding, 0.90 - 2*spectral_height, plot_width, spectral_height])
            ax2.plot(x_axis, original_mean, 'r-', lw=1, alpha=0.5, label='Orig')
            ax2.plot(x_axis, corrected_mean, 'g-', lw=1.5, label='Corr')
            ax2.set_title('Corrected Mean', fontsize=10)
            ax2.set_ylim(min_val, max_val)
            ax2.legend(fontsize=8, loc='upper right')
            ax2.grid(True, alpha=0.3)
            ax2.tick_params(labelsize=8)
            
            # Components
            for i in range(n_components_show):
                axi = fig.add_axes([padding, 0.90 - (2+i+1)*spectral_height, plot_width, spectral_height])
                if len(x_axis) > 0 and len(components_to_plot[i]) == len(x_axis):
                    axi.plot(x_axis, components_to_plot[i], 'k-', lw=1)
                elif len(x_axis) > 0:
                    axi.plot(components_to_plot[i], 'k-', lw=1)
                else:
                    axi.plot(components_to_plot[i], 'k-', lw=1)
                label = f"Component {i+1}\n(Var: {variance_ratio_to_plot[i]*100:.1f}%)"
                axi.text(0.03, 0.9, label, transform=axi.transAxes, fontsize=8, va='top')
                axi.grid(True, alpha=0.3)
                axi.tick_params(labelsize=8)
                if i < n_components_show - 1:
                    axi.set_xticklabels([])
            
            # ==================================================================
            # COLUMN 2: WATERFALL + BEFORE HEATMAP (TOP-MID LEFT)
            # ==================================================================
            ax_waterfall = fig.add_axes([2*padding + plot_width, 0.55, plot_width, 0.35])
            
            # Get the correct decomposition for this specific mission/telescope
            plot_decomposition = self.get_decomposition_for_mission(mission_id, telescope, mission_decompositions)
            
            # Display SKYCHOPDIFF spectra used for this mission/telescope decomposition
            mission_telescope_key = (mission_id, telescope)
            if mission_telescope_key in skychopdiff_spectra:
                to_display = skychopdiff_spectra[mission_telescope_key]
                waterfall_title = f'SKYCHOPDIFF Spectra - {mission_id}/{telescope}\n({len(to_display)} spectra used for decomposition)'
            else:
                # Fallback: show all spectra from current scan
                to_display = original_spectra_subset
                waterfall_title = f'Spectra - {mission_id}/{telescope}\n({len(to_display)} spectra from this scan)'
            
            n_display_waterfall = len(to_display)  # Display all SKYCHOPDIFF spectra
            
            # Safety check: ensure x_axis is valid before using it
            if len(x_axis) > 0:
                im_waterfall = ax_waterfall.imshow(to_display[:n_display_waterfall], aspect='auto',
                                                  extent=[x_axis[0], x_axis[-1], n_display_waterfall, 0],
                                                  interpolation='nearest', cmap='viridis')
            else:
                logger.warning(f"Skipping waterfall plot: x_axis is empty")
                im_waterfall = ax_waterfall.imshow(to_display[:n_display_waterfall], aspect='auto',
                                                  interpolation='nearest', cmap='viridis')
            ax_waterfall.set_title(waterfall_title, fontsize=10)
            ax_waterfall.set_xlabel(x_label, fontsize=9)
            ax_waterfall.set_ylabel('Spectrum #', fontsize=9)
            ax_waterfall.set_ylim(n_display_waterfall, 0)  # Top to bottom
            fig.colorbar(im_waterfall, ax=ax_waterfall, pad=0.02)
            ax_waterfall.tick_params(labelsize=8)
            
            # Calculate n_display early - all plots will use this for consistent sizing
            # Display only the usable processed spectra (those with correction details)
            n_display = len(usable_global_indices)
            
            # Before correction heatmap with line detection contours overlay
            threshold_2d = None
            n_sampled_for_threshold = None
            
            if group_thresholds is not None and len(group_thresholds) > 0:
                # Convert scan/telescope to strings for consistent comparison
                # (handles both byte strings and regular strings from FITS)
                def to_string(x):
                    if isinstance(x, bytes):
                        return x.decode('utf-8').strip()
                    return str(x).strip()
                
                scan_str = to_string(scan)
                telescope_str = to_string(telescope)
                mission_str = to_string(mission_id)
                
                # Consolidate thresholds across ALL subscans for this (mission, telescope, scan) group
                # group_thresholds keys are: (mission_id, telescope, scan, subscan)
                # We want to collect all thresholds that match (mission, telescope, scan)
                found_thresholds = []
                for group_key, thresh in group_thresholds.items():
                    group_mission_id, group_telescope, group_scan, group_subscan = group_key
                    # Match on mission, scan, and telescope
                    # Convert group values to strings too to handle np.str_ and bytes
                    if (to_string(group_mission_id) == mission_str and 
                        to_string(group_scan) == scan_str and 
                        to_string(group_telescope) == telescope_str):
                        found_thresholds.append((group_subscan, thresh))
                
                # If we found thresholds, vertically stack them (accumulating across subscans)
                if found_thresholds:
                    logger.info(f"Found {len(found_thresholds)} threshold sets for mission_id={mission_str}, scan={scan_str}, telescope={telescope_str}")
                    thresholds_to_stack = [thresh for _, thresh in sorted(found_thresholds)]
                    try:
                        threshold_2d = np.vstack(thresholds_to_stack)
                        n_sampled_for_threshold = threshold_2d.shape[0]
                        logger.info(f"Consolidated thresholds shape: {threshold_2d.shape}")
                    except Exception as e:
                        logger.debug(f"Could not consolidate thresholds: {e}")
                        threshold_2d = None
                        n_sampled_for_threshold = None
                else:
                    logger.debug(f"No thresholds found for mission_id={mission_str}, scan={scan_str}, telescope={telescope_str}")
            
            ax_before = fig.add_axes([2*padding + plot_width, 0.05, plot_width, 0.45])
            
            # Display original spectra as heatmap
            before_display = original_spectra_subset[:n_display].copy()
            # CRITICAL: extent MUST match the actual data shape being displayed
            # Shape of before_display: (n_display, n_channels)
            # extent format: [left, right, bottom, top] where bottom and top are in data coordinates
            if len(x_axis) > 0:
                im_before = ax_before.imshow(before_display, aspect='auto',
                                            extent=[x_axis[0], x_axis[-1], n_display, 0],
                                            interpolation='nearest', cmap='viridis')
            else:
                im_before = ax_before.imshow(before_display, aspect='auto',
                                            interpolation='nearest', cmap='viridis')

            
            # Draw line detection contours as overlay if available
            # Use group_thresholds (2D detection result) for visualization,
            # even if science_line_mask is None (which happens when no lines applied to individual spectra)
            if threshold_2d is not None and n_sampled_for_threshold is not None:
                logger.info(f"Attempting to draw contours for plot (threshold_2d shape: {threshold_2d.shape})")
                try:
                    # Use actual 2D thresholds from group_thresholds if available
                    # This shows the actual line detections from OpenCV
                    logger.debug(f"Drawing contours: threshold_2d shape={threshold_2d.shape}, n_sampled={n_sampled_for_threshold}")
                    from matplotlib.patches import Rectangle
                    
                    # Use threshold directly for the rows that were sampled
                    # Don't interpolate - just overlay it on the first n_sampled rows
                    threshold_2d_resized = threshold_2d[:min(n_sampled_for_threshold, n_display)].astype(bool)
                    
                    logger.debug(f"Drawing contours for {threshold_2d_resized.shape[0]} spectra with threshold")
                    
                    # Draw the 2D threshold mask as colored rectangles
                    # For each row that has any True values, draw rectangles at the boundaries
                    has_line_contours = False
                    rect_count = 0
                    for spectrum_idx in range(threshold_2d_resized.shape[0]):
                        row = threshold_2d_resized[spectrum_idx]
                        if np.any(row):
                            has_line_contours = True
                            # Find contiguous regions in this row
                            diff = np.diff(row.astype(int))
                            # +1 means False to True transition, -1 means True to False
                            starts = np.where(diff == 1)[0] + 1
                            ends = np.where(diff == -1)[0]
                            
                            # Handle edge cases
                            if row[0]:
                                starts = np.concatenate([[0], starts])
                            if row[-1]:
                                ends = np.concatenate([ends, [len(row)]])
                            
                            # Draw rectangles for each detected region
                            for start, end in zip(starts, ends):
                                # Convert channel indices to velocity/frequency coordinates
                                if len(x_axis) > 0:
                                    x_start = np.interp(start, [0, len(row)], [x_axis[0], x_axis[-1]])
                                    x_end = np.interp(end, [0, len(row)], [x_axis[0], x_axis[-1]])
                                    x_width = x_end - x_start
                                else:
                                    # Fall back to channel indices if x_axis is empty
                                    x_start = start
                                    x_end = end
                                    x_width = x_end - x_start
                                
                                # Draw rectangle: (x, y, width, height)
                                # y is spectrum_idx (0 at top), height is 1 spectrum
                                # Note: imshow with extent puts y=0 at top and increases downward
                                rect = Rectangle(
                                    (x_start, spectrum_idx),  # lower left corner
                                    x_width,  # width
                                    1.0,  # height (1 spectrum)
                                    linewidth=1.5,
                                    edgecolor='yellow',
                                    facecolor='yellow',
                                    alpha=0.3
                                )
                                ax_before.add_patch(rect)
                                rect_count += 1
                    
                    logger.info(f"Contours drawn: {rect_count} rectangles, has_line_contours={has_line_contours}")
                except Exception as e:
                    import traceback
                    logger.debug(f"Error drawing line contours: {e}")
                    logger.debug(traceback.format_exc())
                    has_line_contours = False
            else:
                logger.debug(f"Skipping contours: science_line_mask is {science_line_mask is None}")
                has_line_contours = False
            
            # Add annotation showing line detection status
            title_text = f'Before Correction\n(first {n_display} specs)'
            if has_line_contours and science_line_mask is not None and len(science_line_mask) > 0:
                try:
                    line_channels = np.any(science_line_mask, axis=0)
                    n_line_channels = np.sum(line_channels)
                    title_text += f'\nYellow contours = {n_line_channels} detected line channels'
                except (IndexError, ValueError) as e:
                    logger.debug(f"Could not annotate line channels: {e}")
                    # Use threshold_2d instead if available
                    if threshold_2d is not None:
                        line_channels = np.any(threshold_2d, axis=0)
                        n_line_channels = np.sum(line_channels)
                        title_text += f'\nYellow contours = {n_line_channels} detected line channels'
            
            ax_before.set_title(title_text, fontsize=10)
            ax_before.set_xlabel(x_label, fontsize=9)
            ax_before.set_ylabel('Spectrum #', fontsize=9)
            # IMPORTANT: Explicitly set y-axis limits to match extent
            ax_before.set_ylim(n_display, 0)  # Top to bottom (inverted for imshow)
            fig.colorbar(im_before, ax=ax_before, pad=0.02)
            ax_before.tick_params(labelsize=8)
            
            # ==================================================================
            # COLUMN 3: DIFFERENCE + AFTER HEATMAP (TOP-MID)
            # ==================================================================
            # Difference heatmap
            ax_diff = fig.add_axes([3*padding + 2*plot_width, 0.55, plot_width, 0.35])
            difference = original_spectra_subset[:n_display] - corrected_spectra_subset[:n_display]
            if len(x_axis) > 0:
                im_diff = ax_diff.imshow(difference, aspect='auto',
                                        extent=[x_axis[0], x_axis[-1], n_display, 0],
                                        interpolation='nearest', cmap='RdBu_r')
            else:
                im_diff = ax_diff.imshow(difference, aspect='auto',
                                        interpolation='nearest', cmap='RdBu_r')
            ax_diff.set_title(f'Correction Applied\n(Original - Corrected)', fontsize=10)
            ax_diff.set_xlabel(x_label, fontsize=9)
            ax_diff.set_ylabel('Spectrum #', fontsize=9)
            ax_diff.set_ylim(n_display, 0)  # Top to bottom
            fig.colorbar(im_diff, ax=ax_diff, pad=0.02)
            ax_diff.tick_params(labelsize=8)
            
            # After correction heatmap
            ax_after = fig.add_axes([3*padding + 2*plot_width, 0.05, plot_width, 0.45])
            if len(x_axis) > 0:
                im_after = ax_after.imshow(corrected_spectra_subset[:n_display], aspect='auto',
                                          extent=[x_axis[0], x_axis[-1], n_display, 0],
                                          interpolation='nearest', cmap='viridis')
            else:
                im_after = ax_after.imshow(corrected_spectra_subset[:n_display], aspect='auto',
                                          interpolation='nearest', cmap='viridis')
            ax_after.set_title(f'After Correction\n(first {n_display} specs)', fontsize=10)
            ax_after.set_xlabel(x_label, fontsize=9)
            ax_after.set_ylabel('Spectrum #', fontsize=9)
            ax_after.set_ylim(n_display, 0)  # Top to bottom
            fig.colorbar(im_after, ax=ax_after, pad=0.02)
            ax_after.tick_params(labelsize=8)
            
            # ==================================================================
            # COLUMN 4: MASKED ARRAY + COEFFICIENTS HEATMAP (TOP-RIGHT)
            # ==================================================================
            # Masked array visualization: show science spectra with masks applied
            # This is the same science spectra as "Before Correction" but with masks shown
            ax_masked = fig.add_axes([4*padding + 3*plot_width, 0.55, plot_width, 0.35])
            
            # Use the same science spectra as the Before Correction plot
            masked_science_spectra = original_spectra_subset[:n_display].copy()
            
            # ====================================================================
            # IMPORTANT: Apply masks PER-SPECTRUM, not as union!
            # This matches the original pca_correct.py behavior where each spectrum
            # has its own apply_mask = bad_channels OR detected_lines[spec_idx]
            # ====================================================================
            
            # 1. Apply BAD CHANNELS uniformly (vertical lines - same for all spectra)
            # These are channels that are NaN/Inf/zero in the original data
            bad_channel_mask = np.isnan(original_spectra_subset[:n_display]) | \
                               np.isinf(original_spectra_subset[:n_display]) | \
                               (original_spectra_subset[:n_display] == 0)
            # Get which channels are bad in ANY spectrum (will show as vertical lines)
            bad_channels_any = np.any(bad_channel_mask, axis=0)
            masked_science_spectra[:, bad_channels_any] = 0
            
            # 2. Apply TELLURIC CHANNELS uniformly (vertical lines - same for all spectra)
            if telluric_line_mask is not None:
                # Get channels that are masked in ANY spectrum (uniform across all)
                telluric_masked_channels = np.any(telluric_line_mask, axis=0)
                masked_science_spectra[:, telluric_masked_channels] = 0
            
            # 3. Apply DETECTED LINE CHANNELS PER-SPECTRUM (different for each row)
            # This is the key difference: detected lines vary by spectrum
            if threshold_2d is not None and len(threshold_2d) > 0:
                # Apply each spectrum's detected lines to that row only
                for spec_idx in range(min(len(threshold_2d), n_display)):
                    detected_channels = threshold_2d[spec_idx].astype(bool)
                    if np.any(detected_channels):
                        masked_science_spectra[spec_idx, detected_channels] = 0
            
            im_masked = ax_masked.imshow(masked_science_spectra, aspect='auto',
                                        extent=[x_axis[0], x_axis[-1], n_display, 0] if len(x_axis) > 0 else None,
                                        interpolation='nearest', cmap='viridis')
            
            # Draw the same line detection contours as in the Before Correction plot
            if threshold_2d is not None and n_sampled_for_threshold is not None:
                try:
                    from matplotlib.patches import Rectangle
                    threshold_2d_resized = threshold_2d[:min(n_sampled_for_threshold, n_display)].astype(bool)
                    rect_count = 0
                    for spectrum_idx in range(threshold_2d_resized.shape[0]):
                        row = threshold_2d_resized[spectrum_idx]
                        if np.any(row):
                            # Find contiguous regions in this row
                            diff = np.diff(row.astype(int))
                            starts = np.where(diff == 1)[0] + 1
                            ends = np.where(diff == -1)[0]
                            
                            # Handle edge cases
                            if row[0]:
                                starts = np.concatenate([[0], starts])
                            if row[-1]:
                                ends = np.concatenate([ends, [len(row)]])
                            
                            # Draw rectangles for each detected region
                            for start, end in zip(starts, ends):
                                # Convert channel indices to velocity/frequency coordinates
                                if len(x_axis) > 0:
                                    x_start = np.interp(start, [0, len(row)], [x_axis[0], x_axis[-1]])
                                    x_end = np.interp(end, [0, len(row)], [x_axis[0], x_axis[-1]])
                                    x_width = x_end - x_start
                                else:
                                    x_start = start
                                    x_end = end
                                    x_width = x_end - x_start
                                
                                # Draw rectangle: (x, y, width, height)
                                rect = Rectangle(
                                    (x_start, spectrum_idx),  # lower left corner
                                    x_width,  # width
                                    1.0,  # height (1 spectrum)
                                    linewidth=1.5,
                                    edgecolor='yellow',
                                    facecolor='yellow',
                                    alpha=0.3
                                )
                                ax_masked.add_patch(rect)
                                rect_count += 1
                except Exception as e:
                    logger.debug(f"Error drawing contours on masked plot: {e}")
            
            ax_masked.set_title(f'Spectra Used for Fitting\n({n_display} specs with masks)', fontsize=10)
            ax_masked.set_xlabel(x_label, fontsize=9)
            ax_masked.set_ylabel('Spectrum #', fontsize=9)
            ax_masked.set_ylim(n_display, 0)  # Top to bottom - same height as Before Correction
            fig.colorbar(im_masked, ax=ax_masked, pad=0.02)
            ax_masked.tick_params(labelsize=8)
            
            # Coefficients heatmap
            # Plot coefficient magnitudes heatmap
            ax_coeff = fig.add_axes([4*padding + 3*plot_width, 0.05, plot_width, 0.45])
            all_coefficients = []
            for idx in usable_global_indices:
                if idx in correction_details and 'coefficients' in correction_details[idx]:
                    # Get the coefficients for this spectrum
                    all_coefficients.append(correction_details[idx]['coefficients'])
                else:
                    # Zero coefficients for spectra without correction details
                    all_coefficients.append(np.zeros(len(self.components)))
            
            all_coefficients = np.array(all_coefficients)[:n_display]
            coeff_magnitudes = np.abs(all_coefficients)
            
            # Display with component number on X-axis, spectrum number on Y-axis
            im_coeff = ax_coeff.imshow(coeff_magnitudes, aspect='auto',
                                      extent=[0, len(self.components), n_display, 0],
                                      interpolation='nearest', cmap='viridis')
            ax_coeff.set_title(f'Coefficient Magnitudes', fontsize=10)
            ax_coeff.set_xlabel('Component #', fontsize=9)
            ax_coeff.set_ylabel('Spectrum #', fontsize=9)
            fig.colorbar(im_coeff, ax=ax_coeff, pad=0.02)
            ax_coeff.tick_params(labelsize=8)
            
            # ==================================================================
            # COLUMN 5: EXAMPLE SPECTRUM BREAKDOWN (RIGHT SIDE)
            # ==================================================================
            # Find example spectrum (median correction)
            # IMPORTANT: example_idx must be within usable_global_indices, not all spectra
            # Only use indices that actually have correction details
            valid_indices_for_example = [
                i for i in range(len(original_spectra_subset))
                if i < len(usable_global_indices)  # Ensure index exists in usable list
            ]
            
            if len(valid_indices_for_example) == 0:
                # If no valid indices, use the first one
                logger.debug(f"No valid indices for example spectrum, using index 0")
                example_idx = 0
                example_global_idx = usable_global_indices[0] if len(usable_global_indices) > 0 else 0
            else:
                correction_amounts = np.array([
                    np.mean(np.abs(original_spectra_subset[i] - corrected_spectra_subset[i]))
                    for i in valid_indices_for_example
                ])
                
                if np.any(correction_amounts > 0):
                    # Get median index in valid list
                    median_pos = len(correction_amounts) // 2
                    example_idx_in_valid = np.argsort(correction_amounts)[median_pos]
                    example_idx = valid_indices_for_example[example_idx_in_valid]
                else:
                    example_idx = valid_indices_for_example[0]
                
                # Get the global index
                if example_idx < len(usable_global_indices):
                    example_global_idx = usable_global_indices[example_idx]
                else:
                    example_global_idx = usable_global_indices[0] if len(usable_global_indices) > 0 else 0
            
            original_ex = original_spectra_subset[example_idx]
            corrected_ex = corrected_spectra_subset[example_idx]
            
            # Get component info for this spectrum
            if example_global_idx in correction_details:
                details = correction_details[example_global_idx]
                comp_info = details.get('component_info', {})
                coeffs = details.get('coefficients', np.zeros(len(self.components)))
            else:
                comp_info = {}
                coeffs = np.zeros(len(self.components))
            
            # Plot correction for each component (show all components)
            n_comp_show = len(self.components)
            # Allocate space evenly for all components (no separate title, use title area for first component)
            spec_height = 0.90 / n_comp_show
            
            for comp_idx in range(n_comp_show):
                # First component gets title above it
                y_pos = 0.90 - (comp_idx + 1)*spec_height
                
                ax_ex_comp = fig.add_axes([5*padding + 4*plot_width, 
                                          y_pos, 
                                          plot_width, spec_height])
                
                # Plot original, corrected, and component contribution
                ax_ex_comp.plot(x_axis, original_ex, 'k-', lw=1.5, alpha=0.8, label='Orig', color="black", zorder=1)
                ax_ex_comp.plot(x_axis, corrected_ex, color='gray', lw=0.5, alpha=0.5, label='Corr', zorder=2)
                
                scaled_comp = coeffs[comp_idx] * self.components[comp_idx]
                ax_ex_comp.plot(x_axis, scaled_comp, 'b-', lw=1, label='Comp')
                
                info = comp_info.get(comp_idx, {})
                used = info.get('used', False)
                coeff_val = coeffs[comp_idx]
                noise_ratio = info.get('noise_ratio', 0)
                
                # Create label
                label_text = f"C{comp_idx+1}: coeff={coeff_val:.2f}, NR={noise_ratio:.1f}"
                if not used:
                    label_text += " (NOT USED)"
                    ax_ex_comp.text(0.03, 0.85, label_text, transform=ax_ex_comp.transAxes,
                                  fontsize=7, va='top', color='red', fontweight='bold')
                else:
                    ax_ex_comp.text(0.03, 0.85, label_text, transform=ax_ex_comp.transAxes,
                                  fontsize=7, va='top', color='green', fontweight='bold')
                
                # Add title above first component
                if comp_idx == 0:
                    ax_ex_comp.text(0.5, 1.15, f'Example Spectrum #{example_global_idx}\n(Median Correction)',
                                  transform=ax_ex_comp.transAxes, ha='center', va='bottom',
                                  fontsize=9, fontweight='bold')
                
                ax_ex_comp.grid(True, alpha=0.3)
                ax_ex_comp.tick_params(labelsize=7)
                
                if comp_idx == 0:
                    ax_ex_comp.legend(fontsize=7, loc='upper right')
                
                if comp_idx == n_comp_show - 1:
                    ax_ex_comp.set_xlabel(x_label, fontsize=8)
            
            fig.suptitle(f'{mission_id} | Scan {scan} | Telescope {telescope} | {len(usable_global_indices)} Spectra', 
                        fontsize=12, fontweight='bold')
            
            # Save plot - one plot per (mission_id, telescope, scan) group
            plot_file = os.path.join(output_dir, 
                                    f'pca_correction_{mission_id}_scan{scan}_{telescope}.png')
            fig.savefig(plot_file, dpi=100, bbox_inches='tight')
            plt.close(fig)
        
        logger.info(f"✓ Generated {total_plots} diagnostic plots")


def correct_fits_file(input_fits, output_fits, decomposition_pkl,
                     cutoff_variance=None, cutoff_noise_ratio=None,
                     object_filter=None, generate_plots=False,
                     output_dir='output/pca_corrected', overwrite=False,
                     scan_filter=None, subscan_filter=None, telescope_filter=None):
    """
    Convenience function to correct FITS file with PCA.
    
    Parameters
    ----------
    input_fits : str
        Input FITS file
    output_fits : str
        Output FITS file (corrected)
    decomposition_pkl : str
        Path to PCA decomposition pickle file
    cutoff_variance : float, optional
        Only use components explaining > this fraction of variance
    cutoff_noise_ratio : float, optional
        Skip components where noise_ratio > this threshold
    object_filter : str, optional
        Only correct spectra matching this OBJECT
    generate_plots : bool
        Generate diagnostic plots
    output_dir : str
        Output directory for plots
    overwrite : bool
        Overwrite output file if exists
    scan_filter : int, optional
        Only correct spectra from this SCAN (for testing)
    subscan_filter : int, optional
        Only correct spectra from this SUBSCAN (for testing)
    telescope_filter : str, optional
        Only correct spectra from this TELESCOP (for testing)
    
    Returns
    -------
    dict
        Correction statistics
    """
    corrector = PCACorrector(decomposition_pkl)
    stats = corrector.correct_fits_file(
        input_fits, output_fits,
        cutoff_variance=cutoff_variance,
        cutoff_noise_ratio=cutoff_noise_ratio,
        object_filter=object_filter,
        generate_plots=generate_plots,
        output_dir=output_dir,
        overwrite=overwrite,
        scan_filter=scan_filter,
        subscan_filter=subscan_filter,
        telescope_filter=telescope_filter
    )
    return stats


def main_correct_cli():
    """
    Command-line interface for PCA-based spectral correction.
    
    Usage:
        pca_correct --input reduced_data.fits --decomposition decomp.pkl --output corrected.fits
    """
    parser = argparse.ArgumentParser(
        description='Apply PCA-based correction to spectral data in FITS files',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Correct using decomposition from same mission
  pca_correct --input reduced_data.fits --decomposition decomp.pkl \\
              --output corrected.fits --plot

  # Correct specific object with cutoffs
  pca_correct --input reduced_data.fits \\
              --decomposition decomp.pkl \\
              --output corrected.fits \\
              --object M51CENTER \\
              --variance-cutoff 0.01 \\
              --noise-ratio-cutoff 3.0 \\
              --plot
        """
    )
    
    # Configuration file
    parser.add_argument(
        '--config',
        default='config.toml',
        help='Configuration file (TOML format, default: config.toml)'
    )
    
    # Input/output arguments
    parser.add_argument(
        '--input', '-i',
        default=None,
        help='Input FITS file with spectra to correct (default: from config [output][prepared_for_pca] or [output][reduced_fits])'
    )
    parser.add_argument(
        '--decomposition', '-d',
        default=None,
        help='Path to PCA decomposition pickle file (default: ./output/pca_components/decomposition_*.pkl)'
    )
    parser.add_argument(
        '--output', '-o',
        default=None,
        help='Output FITS file (corrected spectra) (default: from config [output][pcad_fits])'
    )
    
    # Correction parameters
    parser.add_argument(
        '--variance-cutoff',
        type=float,
        default=None,
        help='Only use components explaining > this fraction of variance (0-1)'
    )
    parser.add_argument(
        '--noise-ratio-cutoff',
        type=float,
        default=None,
        help='Skip components where noise_ratio > this threshold'
    )
    parser.add_argument(
        '--object',
        type=str,
        default=None,
        help='Only correct spectra matching this OBJECT substring (default: from config [parameters][object])'
    )
    parser.add_argument(
        '--hdu',
        type=int,
        default=1,
        help='HDU index containing spectral data (default: 1)'
    )
    parser.add_argument(
        '--spectrum-column',
        default='SPECTRUM',
        help='Name of spectrum column (default: SPECTRUM)'
    )
    
    # Line detection options
    parser.add_argument(
        '--detect-science-lines',
        action='store_true',
        default=True,
        help='Detect science lines using OpenCV (default: True)'
    )
    parser.add_argument(
        '--no-line-detection',
        action='store_false',
        dest='detect_science_lines',
        help='Disable science line detection'
    )
    parser.add_argument(
        '--line-kernel-size',
        type=int,
        default=51,
        help='Kernel size for line detection (must be odd, default: 51)'
    )
    parser.add_argument(
        '--line-cutoff-std',
        type=float,
        default=2.0,
        help='Cutoff in sigma for line detection (default: 2.0)'
    )
    parser.add_argument(
        '--smoothing-kernel',
        type=int,
        default=None,
        help='Smoothing kernel size for line detection refinement (optional)'
    )
    parser.add_argument(
        '--config-window',
        type=float,
        nargs=2,
        default=None,
        metavar=('V_MIN', 'V_MAX'),
        help='Fallback velocity window in km/s if line detection fails (e.g., 450 500)'
    )
    
    # Testing/filtering options
    parser.add_argument(
        '--scan',
        type=int,
        default=None,
        help='Only process spectra from this SCAN (for testing specific combinations)'
    )
    parser.add_argument(
        '--subscan',
        type=int,
        default=None,
        help='Only process spectra from this SUBSCAN (for testing specific combinations)'
    )
    parser.add_argument(
        '--telescope',
        type=str,
        default=None,
        help='Only process spectra from this TELESCOP (for testing specific combinations)'
    )
    parser.add_argument(
        '--mission-id',
        type=str,
        default=None,
        help='Only process spectra from this MISSION_ID (e.g., 2016-05-18_GR_F298)'
    )
    
    # Output options (plot-dir default will be set from config)
    parser.add_argument(
        '--plot',
        action='store_true',
        help='Generate diagnostic plots per scan/telescope'
    )
    parser.add_argument(
        '--plot-dir',
        default=None,  # Will be set from config below
        help='Directory for output plots (default: from config [output][pca_plots_dir])'
    )
    parser.add_argument(
        '--no-overwrite',
        action='store_true',
        help='Do NOT overwrite output file if exists (default: overwrite enabled)'
    )
    
    # Logging
    parser.add_argument(
        '--verbose', '-v',
        action='store_true',
        help='Verbose output'
    )
    parser.add_argument(
        '--debug',
        action='store_true',
        help='Debug output'
    )
    
    args = parser.parse_args()
    
    # Load config to get defaults
    config_plot_dir = 'output/pca_corrected'  # Fallback default
    config_input_file = None
    config_object = None
    config_output_file = None
    config_line_window = None
    config_variance_cutoff = None
    config_noise_ratio_cutoff = None
    config_line_kernel_size = 51
    config_line_cutoff_std = 2.0
    config_smoothing_kernel = None
    
    try:
        config_path = Path(args.config)
        if config_path.exists():
            logger.info(f"Loading configuration from {args.config}")
            config = ConfigLoader(args.config)
            output_config = config.get('output', {})
            config_plot_dir = output_config.get('pca_plots_dir', 'output/pca_corrected')
            # Try prepared_for_pca first, then fall back to reduced_fits
            config_input_file = output_config.get('prepared_for_pca') or output_config.get('reduced_fits')
            config_output_file = output_config.get('pcad_fits')
            
            # Get object filter from [parameters] section
            params_config = config.get('parameters', {})
            config_object = params_config.get('object')
            
            # Get line_window from [reduction] section for science line fallback
            reduction_config = config.get('reduction', {})
            config_line_window = reduction_config.get('line_window')
            
            # Get PCA correction parameters from [pca] section
            pca_config = config.get('pca', {})
            
            # Read cutoff parameters (can be False, float, or int)
            cutoff_val = pca_config.get('cutoff', False)
            if cutoff_val and cutoff_val is not False:
                try:
                    config_variance_cutoff = float(cutoff_val)
                except (ValueError, TypeError):
                    config_variance_cutoff = None
            
            noise_ratio_val = pca_config.get('noise_ratio_cutoff', False)
            if noise_ratio_val and noise_ratio_val is not False:
                try:
                    config_noise_ratio_cutoff = float(noise_ratio_val)
                except (ValueError, TypeError):
                    config_noise_ratio_cutoff = None
            
            # NOTE: global_noise_ratio_cutoff in the original pca_correct.py triggers a sophisticated
            # Gaussian-based component selection algorithm, NOT a direct threshold cutoff.
            # For now, we don't use it automatically. To enable component filtering, users should
            # explicitly set 'noise_ratio_cutoff' (direct threshold) or provide --noise-ratio-cutoff CLI arg.
            # global_noise_ratio = pca_config.get('global_noise_ratio_cutoff', False)
            # For now, leave config_noise_ratio_cutoff as None to use all components by default.
            
            # Read line detection parameters
            line_kernel_val = pca_config.get('line_kernel_size', 51)
            if line_kernel_val:
                try:
                    config_line_kernel_size = int(line_kernel_val)
                except (ValueError, TypeError):
                    config_line_kernel_size = 51
            
            line_cutoff_val = pca_config.get('line_cutoff_std', 2.0)
            if line_cutoff_val:
                try:
                    config_line_cutoff_std = float(line_cutoff_val)
                except (ValueError, TypeError):
                    config_line_cutoff_std = 2.0
            
            smoothing_val = pca_config.get('smoothing_kernel_size', None)
            if smoothing_val:
                try:
                    config_smoothing_kernel = int(smoothing_val)
                except (ValueError, TypeError):
                    config_smoothing_kernel = None
            
            logger.debug(f"PCA plots directory from config: {config_plot_dir}")
            if config_input_file:
                logger.debug(f"Input FITS file from config: {config_input_file}")
            if config_object:
                logger.debug(f"Object filter from config: {config_object}")
            if config_output_file:
                logger.debug(f"Output FITS file from config: {config_output_file}")
            if config_line_window:
                logger.debug(f"Line window from config: {config_line_window}")
            if config_variance_cutoff is not None:
                logger.debug(f"Variance cutoff from config: {config_variance_cutoff}")
            if config_noise_ratio_cutoff is not None:
                logger.debug(f"Noise ratio cutoff from config: {config_noise_ratio_cutoff}")
            logger.debug(f"Line kernel size from config: {config_line_kernel_size}")
            logger.debug(f"Line cutoff std from config: {config_line_cutoff_std}")
            if config_smoothing_kernel is not None:
                logger.debug(f"Smoothing kernel from config: {config_smoothing_kernel}")
        else:
            logger.debug(f"Config file not found: {args.config}, using defaults")
    except Exception as e:
        logger.warning(f"Could not load config file: {e}, using defaults")
    
    # Determine input file: explicit CLI arg > config [output][prepared_for_pca] > config [output][reduced_fits]
    if args.input is not None:
        input_file = args.input
        logger.info(f"✓ Input file = {input_file} (from command line)")
    elif config_input_file:
        input_file = config_input_file
        if config_input_file == config.get('output', {}).get('prepared_for_pca'):
            logger.info(f"✓ Input file = {input_file} (from config [output][prepared_for_pca])")
        else:
            logger.info(f"✓ Input file = {input_file} (from config [output][reduced_fits])")
    else:
        logger.error("No input file specified in configuration or command line")
        logger.error("Specify via:")
        logger.error("  1. Command line: pca_correct --input /path/to/file.fits ...")
        logger.error("  2. Config [output][prepared_for_pca]")
        logger.error("  3. Config [output][reduced_fits]")
        sys.exit(1)
    
    # Determine object filter: explicit CLI arg > config [parameters][object]
    if args.object is not None:
        object_filter = args.object
        logger.info(f"✓ Object filter = {object_filter} (from command line)")
    elif config_object:
        object_filter = config_object
        logger.info(f"✓ Object filter = {object_filter} (from config [parameters][object])")
    else:
        object_filter = None
        logger.info("⚠ No object filter specified - will correct ALL objects in input file")
    
    # Determine decomposition file: explicit CLI arg > ./output/pca_components/
    # First, try to find per-mission decompositions
    mission_decompositions = {}
    decomposition_file = None
    
    if args.decomposition is not None:
        # Single decomposition file specified on CLI
        decomposition_file = args.decomposition
        logger.info(f"✓ Decomposition = {decomposition_file} (from command line)")
    else:
        # Look for per-mission decomposition files in default location
        decomp_dir = Path("./output/pca_components")
        if decomp_dir.exists():
            # Try to load per-mission decompositions
            mission_decompositions = PCACorrector.load_decompositions_from_dir(decomp_dir)
            
            if mission_decompositions:
                logger.info(f"✓ Loaded {len(mission_decompositions)} per-mission decompositions")
            else:
                # Fall back to single decomposition file
                import glob
                decomp_files = sorted(glob.glob(str(decomp_dir / "decomposition_*_components.pkl")))
                if decomp_files:
                    decomposition_file = decomp_files[-1]  # Use most recent
                    logger.info(f"✓ Decomposition = {decomposition_file} (from ./output/pca_components)")
                else:
                    logger.error("No decomposition files found in ./output/pca_components/")
                    logger.error("Specify via:")
                    logger.error("  1. Command line: pca_correct --decomposition /path/to/decomposition.pkl ...")
                    logger.error("  2. Place decomposition files in ./output/pca_components/decomposition_*.pkl")
                    sys.exit(1)
        else:
            logger.error("No decomposition file specified and default location not found")
            logger.error("Specify via:")
            logger.error("  1. Command line: pca_correct --decomposition /path/to/decomposition.pkl ...")
            logger.error("  2. Place decomposition files in ./output/pca_components/decomposition_*.pkl")
            sys.exit(1)
    
    # Determine output file: explicit CLI arg > config [output][pcad_fits]
    if args.output is not None:
        output_file = args.output
        logger.info(f"✓ Output file = {output_file} (from command line)")
    elif config_output_file:
        output_file = config_output_file
        logger.info(f"✓ Output file = {output_file} (from config [output][pcad_fits])")
    else:
        logger.error("No output file specified in configuration or command line")
        logger.error("Specify via:")
        logger.error("  1. Command line: pca_correct --output /path/to/output.fits ...")
        logger.error("  2. Config [output][pcad_fits]")
        sys.exit(1)
    
    # Set plot directory: explicit CLI arg > config > fallback
    if args.plot_dir is None:
        args.plot_dir = config_plot_dir
    
    # Set logging level
    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    elif args.verbose:
        logging.getLogger().setLevel(logging.INFO)
    
    try:
        # Check input files exist
        input_path = Path(input_file)
        
        if not input_path.exists():
            logger.error(f"Input FITS file not found: {input_file}")
            sys.exit(1)
        
        # Check decomposition files
        if decomposition_file:
            decomp_path = Path(decomposition_file)
            if not decomp_path.exists():
                logger.error(f"Decomposition file not found: {decomposition_file}")
                sys.exit(1)
            logger.info(f"Decomposition: {decomposition_file}")
        elif mission_decompositions:
            logger.info(f"Using {len(mission_decompositions)} per-mission decompositions")
        else:
            logger.error("No decomposition file or per-mission decompositions found")
            sys.exit(1)
        
        logger.info(f"Output file: {output_file}")
        
        # Create corrector and apply correction
        logger.info("Initializing PCA corrector...")

        # Resolve final parameter values: CLI args take precedence over config file.
        # This must happen BEFORE constructing PCACorrector so the corrector receives
        # the correct values (previously these merges were done after construction).
        cutoff_variance = args.variance_cutoff if args.variance_cutoff is not None else config_variance_cutoff
        cutoff_noise_ratio = args.noise_ratio_cutoff if args.noise_ratio_cutoff is not None else config_noise_ratio_cutoff
        line_kernel_size = args.line_kernel_size if args.line_kernel_size != 51 else config_line_kernel_size
        line_cutoff_std = args.line_cutoff_std if args.line_cutoff_std != 2.0 else config_line_cutoff_std
        smoothing_kernel = args.smoothing_kernel if args.smoothing_kernel is not None else config_smoothing_kernel

        # Get line_window from [reduction] section (already extracted above)
        line_window_velocities = tuple(config_line_window) if config_line_window else None

        if decomposition_file:
            # Use single decomposition file
            corrector = PCACorrector(
                str(decomp_path),
                line_kernel_size=line_kernel_size,
                line_cutoff_std=line_cutoff_std,
                smoothing_kernel_size=smoothing_kernel,
                line_window_velocities=line_window_velocities
            )
            mission_decompositions_to_use = None
        else:
            # Use first mission's decomposition as default, but will use per-mission in correct_fits_file
            first_mission = list(mission_decompositions.keys())[0]
            # Create a proper PCACorrector instance with the first mission's decomposition
            # We'll create a temporary pickle file in memory to load it
            import tempfile
            temp_decomp = {
                'components': mission_decompositions[first_mission].components,
                'explained_variance_ratio': mission_decompositions[first_mission].explained_variance_ratio,
                'mean_spectrum': mission_decompositions[first_mission].mean_spectrum,
                'metadata': mission_decompositions[first_mission].metadata,
                'config': mission_decompositions[first_mission].metadata.get('config', {}),
                'spectrum_metadata': [],
                'reference_spectrum_indices': mission_decompositions[first_mission].reference_spectrum_indices
            }
            with tempfile.NamedTemporaryFile(mode='wb', suffix='.pkl', delete=False) as tf:
                pickle.dump(temp_decomp, tf)
                temp_decomp_path = tf.name

            try:
                corrector = PCACorrector(
                    temp_decomp_path,
                    line_kernel_size=line_kernel_size,
                    line_cutoff_std=line_cutoff_std,
                    smoothing_kernel_size=smoothing_kernel,
                    line_window_velocities=line_window_velocities
                )
            finally:
                # Clean up temp file
                os.unlink(temp_decomp_path)
            mission_decompositions_to_use = mission_decompositions

        logger.info("Applying PCA correction...")
        
        logger.info(f"Correction parameters:")
        logger.info(f"  Variance cutoff: {cutoff_variance}")
        logger.info(f"  Noise ratio cutoff: {cutoff_noise_ratio}")
        logger.info(f"  Line kernel size: {line_kernel_size}")
        logger.info(f"  Line cutoff std: {line_cutoff_std}")
        logger.info(f"  Smoothing kernel: {smoothing_kernel}")
        
        stats = corrector.correct_fits_file(
            input_fits=str(input_path),
            output_fits=output_file,
            cutoff_variance=cutoff_variance,
            cutoff_noise_ratio=cutoff_noise_ratio,
            hdu_index=args.hdu,
            spectrum_col=args.spectrum_column,
            object_filter=object_filter,
            overwrite=not args.no_overwrite,
            generate_plots=args.plot,
            output_dir=args.plot_dir,
            config_window=config_line_window if config_line_window else (tuple(args.config_window) if args.config_window else None),
            detect_science_lines=args.detect_science_lines,
            scan_filter=args.scan,
            subscan_filter=args.subscan,
            telescope_filter=args.telescope,
            mission_id_filter=args.mission_id,
            mission_decompositions=mission_decompositions_to_use
        )
        # Print summary
        logger.info("\n=== Correction Summary ===")
        logger.info(f"Total spectra: {stats['total']}")
        logger.info(f"Successfully corrected: {stats['corrected']}")
        logger.info(f"Failed: {stats['failed']}")
        
        if stats['components_used']:
            avg_comp = np.mean(stats['components_used'])
            std_comp = np.std(stats['components_used'])
            logger.info(f"Components used: {avg_comp:.1f} ± {std_comp:.1f}")
        
        logger.info(f"\n✓ Corrected FITS saved to: {output_file}")
        if args.plot:
            logger.info(f"✓ Diagnostic plots saved to: {args.plot_dir}")
        
    except Exception as e:
        logger.error(f"Error: {e}")
        if args.debug:
            import traceback
            traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main_correct_cli()
