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


def find_science_lines(spectra, kernel_size=51, cutoff_std=2.0, smoothing_kernel=None):
    """
    Detect science lines in spectra using OpenCV edge detection.
    
    This uses Gaussian blur and threshold detection to identify emission/absorption lines.
    Follows the same algorithm as pca_utilities.find_lines().
    
    Parameters
    ----------
    spectra : ndarray
        2D array of spectra (n_spectra, n_channels)
    kernel_size : int
        Kernel size for Gaussian blur (must be odd, default: 51)
    cutoff_std : float
        Threshold in units of standard deviation (default: 2.0)
    smoothing_kernel : int, optional
        Smoothing kernel size (separate from detection kernel)
    
    Returns
    -------
    ndarray
        2D boolean mask (n_spectra, n_channels) indicating detected lines
    """
    if cv is None:
        logger.warning("OpenCV not available, skipping line detection")
        return None
    
    # Ensure kernel_size is positive odd integer
    if kernel_size <= 0 or kernel_size % 2 == 0:
        kernel_size = max(1, kernel_size - 1)
    
    spectra = np.asarray(spectra)
    if len(spectra.shape) == 1:
        spectra = spectra.reshape(1, -1)
    
    # Normalize spectra to 8-bit range for OpenCV
    data_min = np.min(spectra)
    data_max = np.max(spectra)
    
    if data_max > data_min:
        data_8uc = 255 * (spectra - data_min) / (data_max - data_min)
    else:
        data_8uc = np.zeros_like(spectra)
    
    img = data_8uc.astype(np.uint8)
    img = cv.normalize(img, None, 0, 100, cv.NORM_MINMAX)
    
    # Apply Gaussian blur
    gray = cv.GaussianBlur(img, (kernel_size, kernel_size), 0)
    
    # Detect threshold
    mean = gray[np.where(gray != 0)].mean() if np.any(gray != 0) else 0
    std = gray[np.where(gray != 0)].std() if np.any(gray != 0) else 1
    threshold_val = mean + cutoff_std * std
    
    _, threshold = cv.threshold(
        gray, threshold_val, 1, cv.THRESH_BINARY
    )
    
    initial_detections = np.sum(threshold > 0)
    
    # Iterative refinement
    new_mean = 0
    new_std = 0
    iterations = 0
    max_iterations = 10
    
    while (not np.isclose(new_mean, mean) or not np.isclose(new_std, std)) and iterations < max_iterations:
        new_mean = mean
        new_std = std
        new_gray = gray.copy()
        new_gray[np.where(threshold == 1)] = 0
        mean = new_gray[np.where(new_gray != 0)].mean() if np.any(new_gray != 0) else 0
        std = new_gray[np.where(new_gray != 0)].std() if np.any(new_gray != 0) else 1
        threshold_val = mean + cutoff_std * std
        _, threshold = cv.threshold(
            gray, threshold_val, 1, cv.THRESH_BINARY
        )
        iterations += 1
    
    final_detections = np.sum(threshold > 0)
    logger.debug(f"Line detection converged after {iterations} iterations: {initial_detections} → {final_detections} pixels")
    
    # Log per-spectrum statistics
    mask = threshold.astype(bool)
    line_channels = np.any(mask, axis=0)
    if np.any(line_channels):
        line_indices = np.where(line_channels)[0]
        logger.debug(f"Detected lines in {np.sum(line_channels)} channels: {line_indices[0]}-{line_indices[-1]}")
    
    return mask


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
                 smoothing_kernel_size=None):
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
        """
        self.config = config or {}
        self.line_kernel_size = line_kernel_size
        self.line_cutoff_std = line_cutoff_std
        self.smoothing_kernel_size = smoothing_kernel_size
        
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
        
        for decomp_file in decomp_files:
            try:
                # Extract mission_id and telescope from filename
                # New format: decomposition_{mission_id}_{telescop}_{date}_components.pkl
                # Old format: decomposition_{mission_id}_{date}_components.pkl
                filename = os.path.basename(decomp_file)
                parts = filename.replace('decomposition_', '').replace('_components.pkl', '').split('_')
                
                # Try to identify mission_id, telescope, and date
                # Date can be 7-8 digits (YYYYMMD or YYYYMMDD)
                flight_date = None
                remaining = parts
                if len(parts) >= 2:
                    # Check last part for date (7-8 consecutive digits)
                    if len(parts[-1]) >= 7 and parts[-1][:7].isdigit():
                        flight_date = parts[-1]
                        remaining = parts[:-1]
                
                if flight_date is None:
                    # No valid date found, skip
                    logger.warning(f"  Could not parse date from {filename}")
                    continue
                
                # Now determine if we have mission_id/telescope or just mission_id
                # New files have: mission_id, telescope, date
                # Old files have: mission_id, date (where mission_id may have underscores)
                # FULL format: decomposition_2017-02-01_GR_F367_LFAH_PX00_S_2017021_components.pkl
                # When split by '_': ['2017', '02', '01', 'GR', 'F367', 'LFAH', 'PX00', 'S', '2017021']
                # Need to find LFAH/LFAV which come AFTER F367
                
                # Try to detect if we have mission/telescope pairing
                # Look for telescope patterns - they usually have format like LFAH_PX##_S or LFAV_PX##_S
                mission_id = None
                telescope = None
                
                # Search from the end backwards to find a telescope pattern
                # Telescopes are typically: LFAH_PXxx_S, LFAV_PXxx_S, LFBH_PXxx_S, etc.
                for i in range(len(remaining) - 1, -1, -1):
                    part = remaining[i]
                    # Check if this looks like a 4-letter telescope code at the start
                    if len(part) >= 4 and part[:4] in ['LFAH', 'LFAI', 'LFAV', 'LFBI', 'LFBH', 'PRISM', 'HR', 'R100Q']:
                        # Check if the following parts look like PX##_S pattern
                        # e.g., remaining = [..., 'LFAH', 'PX00', 'S'] or [..., 'LFAH', 'PX00', 'S']
                        # We need: LFAH + PX## + S as telescope
                        if i + 2 <= len(remaining) - 1 and remaining[i+1].startswith('PX') and remaining[i+2] == 'S':
                            # Found LFAH_PX00_S pattern
                            mission_id = '_'.join(remaining[:i])
                            telescope = '_'.join(remaining[i:i+3])  # LFAH_PX00_S
                            break
                        elif i + 1 == len(remaining) - 1:
                            # Maybe it's just LFAH with no PX## suffix
                            mission_id = '_'.join(remaining[:i])
                            telescope = part
                            break
                
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
                else:
                    # Fall back to old format (mission_id only)
                    mission_id = '_'.join(remaining) if not mission_id else mission_id
                    decompositions[mission_id] = corrector
                    logger.info(f"  Loaded {mission_id}: {len(corrector.components)} components")
                
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
            Mission identifier
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
    
    def get_noise_ratio(self, spectrum, component, scaled_spectrum=None):
        """
        Calculate noise ratio: noise in spectrum / noise in component.
        
        Parameters
        ----------
        spectrum : ndarray
            Original spectrum
        component : ndarray
            PCA component
        scaled_spectrum : ndarray, optional
            Scaled spectrum (if different from original)
        
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
        
        # Noise in spectrum and component
        spectrum_std = np.std(scaled_spectrum)
        component_std = np.std(scaled_comp)
        
        if component_std < 1e-10:
            return np.inf
        
        return spectrum_std / component_std
    
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
                         cutoff_variance=None, cutoff_noise_ratio=None, verbose=False):
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
            
            # Calculate noise ratio
            scaled_comp = coeff[i] * comp
            noise_ratio = self.get_noise_ratio(spectrum, comp)
            
            if verbose:
                logger.info(f"  Comp {i}: coeff={coeff[i]:.6f}, var_ratio={var_ratio:.6f}, "
                           f"noise_ratio={noise_ratio:.6f}, scaled_comp_range=[{scaled_comp.min():.6f}, {scaled_comp.max():.6f}]")
            
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
                         telescope_filter=None, mission_decompositions=None):
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
            
            # Get velocity axis if available
            velocity_axis = None
            if 'VELOCITY_AXIS' in hdu.columns.names:
                velocity_axis = data['VELOCITY_AXIS'][0]  # Same for all spectra
                velocity_axis_kms = velocity_axis / 1000.0  # Convert m/s to km/s
                logger.info(f"  Using VELOCITY_AXIS: {len(velocity_axis)} channels")
            else:
                velocity_axis_kms = None
            
            # Detect science lines if requested
            # Note: Line detection is done on a representative sample of spectra
            # per mission/telescope/scan/subscan combination (following pca_correct.py)
            science_line_mask = None
            if detect_science_lines:
                # NOTE: Science line detection happens HERE in pca_correct, not during decomposition.
                # During decomposition, only the mission-specific telluric mask (from mission_id_parameters.yml)
                # is used to exclude telluric regions from the components. Here, we additionally detect
                # and mask the actual emission line being observed to prevent it from biasing the correction.
                logger.info("Detecting science lines per mission/telescope/scan/subscan...")
                try:
                    # Convert string columns to consistent format once
                    def to_string(x):
                        if isinstance(x, bytes):
                            return x.decode('utf-8').strip()
                        return str(x).strip()
                    
                    mission_ids = np.array([to_string(data['MISSION_ID'][i]) for i in indices])
                    telescopes = np.array([to_string(data['TELESCOP'][i]) for i in indices])
                    scans = np.array([to_string(x) for x in data['SCAN'][indices]])
                    subscans = np.array([to_string(x) for x in data['SUBSCAN'][indices]])
                    
                    # Create mask for all spectra
                    science_line_mask = np.zeros((len(indices), len(velocity_axis_kms)), dtype=bool)
                    
                    # Get unique mission/telescope/scan/subscan combinations
                    unique_groups = np.unique(
                        np.column_stack((mission_ids, telescopes, scans, subscans)),
                        axis=0
                    )
                    
                    logger.info(f"  Found {len(unique_groups)} unique mission/telescope/scan/subscan groups")
                    logger.info(f"  Total spectra to process: {len(indices)}")
                    
                    # DEBUG: Show first few groups and their sizes
                    group_sizes = {}
                    for mission_id_group, telescope_group, scan_id, subscan_id in unique_groups:
                        group_mask = (
                            (mission_ids == mission_id_group) &
                            (telescopes == telescope_group) &
                            (scans == scan_id) &
                            (subscans == subscan_id)
                        )
                        group_sizes[(mission_id_group, telescope_group, scan_id, subscan_id)] = np.sum(group_mask)
                    
                    # Show distribution
                    sizes = list(group_sizes.values())
                    logger.info(f"  Group sizes: min={min(sizes)}, max={max(sizes)}, mean={np.mean(sizes):.1f}, median={np.median(sizes):.0f}")
                    
                    # Show a few example groups
                    sample_groups = list(group_sizes.items())[:5]
                    for (mid, tel, scan, subscan), size in sample_groups:
                        logger.info(f"    Example group: {mid}/{tel}/scan{scan}/subscan{subscan} ({size} spectra)")
                    logger.info(f"  velocity_axis_kms is None: {velocity_axis_kms is None}, len: {len(velocity_axis_kms) if velocity_axis_kms is not None else 'N/A'}")
                    
                    group_detections = {}
                    # Store 2D threshold masks for each group for visualization
                    group_thresholds = {}
                    logger.info(f"Starting line detection for {len(unique_groups)} mission/telescope/scan/subscan groups...")
                    # Set seed for reproducible sampling across runs
                    np.random.seed(42)
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
                        
                        # Get sample of spectra from this group (up to 20 to keep computation reasonable)
                        # Use deterministic sampling with seed for reproducibility
                        sample_size = min(20, len(group_indices))
                        sample_indices_in_group = np.random.choice(np.arange(len(group_indices)), sample_size, replace=False)
                        sample_indices = group_indices[sample_indices_in_group]
                        group_sample = data[spectrum_col][indices[sample_indices]]
                        
                        # Detect lines in this group's sample (as 2D image to leverage smoothing/contours)
                        try:
                            group_threshold = find_science_lines(
                                group_sample,
                                kernel_size=self.line_kernel_size,
                                cutoff_std=self.line_cutoff_std,
                                smoothing_kernel=self.smoothing_kernel_size
                            )
                            
                            if group_threshold is not None:
                                # SAVE THE 2D THRESHOLD FOR VISUALIZATION
                                group_key = (mission_id_group, telescope_group, scan_id, subscan_id)
                                group_thresholds[group_key] = group_threshold
                                
                                # Get the detected line region from the 2D threshold
                                line_channels_window = np.any(group_threshold, axis=0)
                                if np.any(line_channels_window):
                                    # Get the continuous channel range
                                    line_indices = np.where(line_channels_window)[0]
                                    ch_min = line_indices[0]
                                    ch_max = line_indices[-1]
                                    
                                    logger.debug(f"  Sample detection for {mission_id_group}/{telescope_group}/scan{scan_id}/subscan{subscan_id}: channels {ch_min}-{ch_max}")
                                    
                                    # Apply the detected line mask to each spectrum in the group
                                    # IMPORTANT: Only apply masks to spectra that were actually sampled and detected
                                    # Do NOT apply fallback masks to unsampled spectra - they should have no mask
                                    n_spectra_with_lines = 0
                                    spectra_with_mask_list = []
                                    
                                    for group_spectrum_idx, global_spectrum_idx in enumerate(group_indices):
                                        # Check if this spectrum was in the sample
                                        # by seeing if its position matches any of the sample indices
                                        sample_position_in_group = np.where(sample_indices_in_group == group_spectrum_idx)[0]
                                        
                                        if len(sample_position_in_group) > 0:
                                            # This spectrum was in the sample, use its detection
                                            threshold_row = sample_position_in_group[0]
                                            line_mask = group_threshold[threshold_row]
                                            
                                            if np.any(line_mask):
                                                # Apply this mask to the spectrum
                                                # Convert global index to local index within filtered indices array
                                                local_idx = np.where(indices == global_spectrum_idx)[0][0]
                                                science_line_mask[local_idx] = line_mask
                                                n_spectra_with_lines += 1
                                                spectra_with_mask_list.append(global_spectrum_idx)
                                        # else: spectrum wasn't sampled, so don't apply any mask
                                    
                                    if n_spectra_with_lines > 0:
                                        # Record detection for logging
                                        group_key = (mission_id_group, telescope_group, scan_id, subscan_id)
                                        group_detections[group_key] = {
                                            'ch_min': ch_min,
                                            'ch_max': ch_max,
                                            'n_channels': ch_max - ch_min + 1,
                                            'n_spectra_with_mask': n_spectra_with_lines
                                        }
                                        logger.debug(f"    Applied to {n_spectra_with_lines} spectra")
                                        logger.debug(f"    Applied to spectra indices: {spectra_with_mask_list}")
                        except Exception as e:
                            logger.debug(f"Line detection failed for {mission_id_group}/{telescope_group}/scan{scan_id}/subscan{subscan_id}: {e}")
                    
                    if group_detections:
                        # Log summary of detections across all groups
                        line_channels_any = np.any(science_line_mask, axis=0)
                        n_line_channels = np.sum(line_channels_any)
                        line_indices = np.where(line_channels_any)[0]
                        ch_min = line_indices[0]
                        ch_max = line_indices[-1]
                        
                        if velocity_axis_kms is not None:
                            v_min = velocity_axis_kms[ch_min]
                            v_max = velocity_axis_kms[ch_max]
                            logger.info(f"  Detected science lines in {n_line_channels} channels "
                                      f"({ch_min}-{ch_max}, {v_min:.1f}-{v_max:.1f} km/s)")
                        else:
                            logger.info(f"  Detected science lines in {n_line_channels} channels ({ch_min}-{ch_max})")
                        
                        # Report per-group statistics
                        logger.info(f"  {len(group_detections)}/{len(unique_groups)} mission/telescope/scan/subscan groups with detected science lines")
                        
                        # DEBUG: Show how many unique mask patterns exist
                        unique_masks = set()
                        for row in science_line_mask:
                            unique_masks.add(tuple(row))
                        logger.info(f"  {len(unique_masks)} unique mask patterns across {len(science_line_mask)} spectra")
                    else:
                        logger.warning("  Could not detect science lines in any mission/telescope/scan/subscan group")
                        science_line_mask = None
                
                except Exception as e:
                    logger.warning(f"Line detection failed: {e}")
                    science_line_mask = None
            
            # NOTE: If no lines detected, don't apply fallback masking
            # The config_window is only used by the old reduction pipeline as a protection window.
            # For PCA correction, we want actual detected lines only, not speculative masking.
            # The config_window is still passed through for reference but not used as a mask.
            
            # Load telluric line mask from mission parameters
            telluric_line_mask = None
            if velocity_axis_kms is not None and mission_id != "UNKNOWN":
                try:
                    telluric_line_mask = get_telluric_line_mask(mission_id, velocity_axis_kms, len(indices))
                    if telluric_line_mask is not None:
                        telluric_channels = np.any(telluric_line_mask, axis=0)
                        n_telluric_channels = np.sum(telluric_channels)
                        logger.info(f"  Loaded telluric line mask: {n_telluric_channels} channels masked")
                        logger.debug(f"    telluric_line_mask shape: {telluric_line_mask.shape}, dtype: {telluric_line_mask.dtype}")
                        # Show which channels are masked
                        masked_ch_indices = np.where(telluric_channels)[0]
                        if len(masked_ch_indices) > 0:
                            logger.debug(f"    Masked channels: {masked_ch_indices[0]}-{masked_ch_indices[-1]} (indices)")
                        
                        # IMPORTANT: Remove telluric lines from science line detection
                        # Following the same pattern as pca_correct.py
                        # This ensures no line detection happens in telluric regions
                        if len(group_thresholds) > 0:
                            for group_key, threshold_2d in group_thresholds.items():
                                # Zero out telluric channels in the 2D threshold
                                threshold_2d[:, telluric_channels] = 0
                            logger.debug(f"  Removed telluric regions from {len(group_thresholds)} group thresholds")
                        
                        # Also remove from science_line_mask
                        for spec_idx in range(len(science_line_mask)):
                            science_line_mask[spec_idx][telluric_channels] = False
                        logger.debug(f"  Removed telluric regions from science_line_mask")
                except Exception as e:
                    logger.debug(f"  Could not load telluric mask: {e}")
            
            # Correct each spectrum
            # Initialize with original data so uncorrected spectra (outside indices) stay unchanged
            corrected_spectra = data[spectrum_col].copy()
            # IMPORTANT: Keep original spectra for plotting BEFORE any modifications
            original_spectra_for_plotting = data[spectrum_col].copy()
            correction_details = {}
            stats = {
                'total': len(indices),
                'corrected': 0,
                'failed': 0,
                'components_used': [],
                'lines_detected': science_line_mask is not None
            }
            
            for spec_idx, idx in enumerate(indices):
                try:
                    spectrum = data[spectrum_col][idx]
                    
                    # Select correct decomposition for this spectrum if using per-mission decompositions
                    if mission_decompositions:
                        spec_mission_id = data['MISSION_ID'][idx].strip()
                        spec_telescope = data['TELESCOP'][idx].strip()
                        decomp = self.get_decomposition_for_mission(spec_mission_id, spec_telescope, mission_decompositions)
                        if decomp is None:
                            logger.warning(f"Spectrum {idx}: no decomposition for {spec_mission_id}/{spec_telescope}, skipping")
                            corrected_spectra[idx] = spectrum
                            stats['failed'] += 1
                            continue
                        # Temporarily swap components for this spectrum
                        saved_components = self.components
                        saved_variance_ratio = self.explained_variance_ratio
                        self.components = decomp.components
                        self.explained_variance_ratio = decomp.explained_variance_ratio
                    
                    # Check for bad channels (NaN, Inf, 0)
                    bad_channels = np.isnan(spectrum) | np.isinf(spectrum) | (spectrum == 0)
                    good_channels_for_fitting = ~bad_channels  # Channels for fitting (no NaN/Inf/0)
                    good_channels_for_subtraction = ~bad_channels.copy()  # Channels where we actually subtract components
                    
                    n_good = np.sum(good_channels_for_fitting)
                    n_line = 0
                    
                    # Apply science line mask - exclude from SUBTRACTION, but still use for FITTING
                    if science_line_mask is not None:
                        # Science lines are protected - don't subtract components there
                        line_regions = science_line_mask[min(spec_idx, len(science_line_mask)-1)]
                        good_channels_for_subtraction = good_channels_for_subtraction & ~line_regions
                        n_line = np.sum(line_regions)
                    
                    # Apply telluric line mask - exclude from SUBTRACTION, but still use for FITTING
                    if telluric_line_mask is not None:
                        telluric_regions = telluric_line_mask[min(spec_idx, len(telluric_line_mask)-1)]
                        good_channels_for_subtraction = good_channels_for_subtraction & ~telluric_regions

                    
                    if not np.any(good_channels_for_fitting):
                        logger.warning(f"Spectrum {idx}: no good channels for fitting (had {n_good} good, {n_line if science_line_mask is not None else 0} lines), skipping")
                        corrected_spectra[idx] = spectrum
                        # Record as skipped but still track it
                        correction_details[idx] = {'n_components_used': 0, 'skipped': True}
                        stats['failed'] += 1
                        # Restore original decomposition if using per-mission
                        if mission_decompositions:
                            self.components = saved_components
                            self.explained_variance_ratio = saved_variance_ratio
                        continue
                    
                    # Apply correction
                    corrected, details = self.apply_correction(
                        spectrum, 
                        good_channels=good_channels_for_fitting,
                        good_channels_for_subtraction=good_channels_for_subtraction,
                        cutoff_variance=cutoff_variance,
                        cutoff_noise_ratio=cutoff_noise_ratio,
                        verbose=(spec_idx < 3)  # Log first 3 spectra
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
            
            # Write output FITS
            os.makedirs(os.path.dirname(output_fits) or '.', exist_ok=True)
            hdul_out = fits.HDUList([hdul[0], fits.BinTableHDU(data, header=header)])
            hdul_out.writeto(output_fits, overwrite=overwrite)
            
            logger.info(f"✓ Saved corrected spectra to {output_fits}")
            logger.info(f"  Corrected: {stats['corrected']}/{stats['total']}")
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
                
                # Note: science_line_mask and telluric_line_mask are already sized for the filtered spectra
                # (they were initialized as np.zeros((len(indices), n_channels))), so we pass them directly
                # without further filtering
                
                # Generate plots with all filtered spectra
                self._generate_plots(filtered_data, filtered_original, filtered_corrected, filtered_details,
                                   velocity_axis, mission_id, output_dir,
                                   science_line_mask=science_line_mask,
                                   telluric_line_mask=telluric_line_mask,
                                   input_fits=input_fits,
                                   group_thresholds=group_thresholds,
                                   global_indices=indices)
            
            return stats
    
    def _generate_plots(self, data, original_spectra, corrected_spectra, correction_details,
                       velocity_axis, mission_id, output_dir, science_line_mask=None, telluric_line_mask=None,
                       input_fits=None, group_thresholds=None, global_indices=None):
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
        
        # Debug logging for masks
        logger.debug(f"_generate_plots called with science_line_mask: {science_line_mask is not None}")
        logger.debug(f"_generate_plots called with telluric_line_mask: {telluric_line_mask is not None}")
        if science_line_mask is not None:
            logger.debug(f"  science_line_mask shape: {science_line_mask.shape}")
        if telluric_line_mask is not None:
            logger.debug(f"  telluric_line_mask shape: {telluric_line_mask.shape}")
        
        # Get unique scan/subscan/telescope combinations
        scans = data['SCAN']
        subscans = data['SUBSCAN']
        telescopes = data['TELESCOP']
        unique_combos = list(set(zip(scans, subscans, telescopes)))
        unique_combos = sorted(unique_combos)
        
        # Limit to ~20 plots by sampling
        max_plots = 20
        if len(unique_combos) > max_plots:
            step = max(1, len(unique_combos) // max_plots)
            unique_combos = unique_combos[::step]
        
        logger.info(f"Creating {len(unique_combos)} plots...")
        
        plot_count = 0
        for scan, subscan, telescope in unique_combos:
            plot_count += 1
            mask = (scans == scan) & (subscans == subscan) & (telescopes == telescope)
            local_indices = np.where(mask)[0]
            
            if len(local_indices) == 0:
                continue
            
            # Convert local indices to global indices if we have the mapping
            if global_indices is not None:
                global_idx_for_group = global_indices[local_indices]
            else:
                global_idx_for_group = local_indices
            
            # Filter to only spectra that were actually processed (exclude skipped ones with all bad channels)
            # Skipped spectra have correction_details[idx]['skipped'] = True
            usable_global_indices = [idx for idx in global_idx_for_group 
                                if idx in correction_details 
                                and not correction_details[idx].get('skipped', False)]
            
            if len(usable_global_indices) == 0:
                logger.info(f"scan={scan}, subscan={subscan}, telescope={telescope}: Found {len(global_idx_for_group)} total, {len([i for i in global_idx_for_group if i in correction_details])} in correction_details, {len(usable_global_indices)} usable (skipping plot)")
                continue
            
            # Convert usable global indices back to local indices for slicing the filtered arrays
            usable_local_indices = [np.where(global_indices == g_idx)[0][0] for g_idx in usable_global_indices] if global_indices is not None else usable_global_indices
            
            logger.info(f"Creating plot for scan={scan}, subscan={subscan}, telescope={telescope}: {len(usable_local_indices)} usable spectra (skipped {len(global_idx_for_group) - len(usable_local_indices)} bad ones)")
            
            # Get spectra for this scan/telescope (only usable ones)
            original_spectra_subset = original_spectra[usable_local_indices]
            corrected_spectra_subset = corrected_spectra[usable_local_indices]
            logger.debug(f"  original_spectra_subset.shape={original_spectra_subset.shape}")
            
            # Setup x_axis
            if velocity_axis is not None:
                x_axis = velocity_axis / 1000.0  # Convert m/s to km/s
                x_label = 'Velocity (km/s)'
            else:
                x_axis = np.arange(len(original_spectra_subset[0]))
                x_label = 'Channel'
            
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
            
            min_val = min(np.min(original_mean), np.min(corrected_mean)) * 1.1
            max_val = max(np.max(original_mean), np.max(corrected_mean)) * 1.1
            
            n_components_show = len(self.components)
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
                axi.plot(x_axis, self.components[i], 'k-', lw=1)
                label = f"Component {i+1}\n(Var: {self.explained_variance_ratio[i]*100:.1f}%)"
                axi.text(0.03, 0.9, label, transform=axi.transAxes, fontsize=8, va='top')
                axi.grid(True, alpha=0.3)
                axi.tick_params(labelsize=8)
                if i < n_components_show - 1:
                    axi.set_xticklabels([])
            
            # ==================================================================
            # COLUMN 2: WATERFALL + BEFORE HEATMAP (TOP-MID LEFT)
            # ==================================================================
            # Heatmap of spectra: prefer SKYCHOPDIFF reference spectra if available
            ax_waterfall = fig.add_axes([2*padding + plot_width, 0.55, plot_width, 0.35])
            
            # Try to load reference SKYCHOPDIFF spectra from FITS using stored indices
            waterfall_spectra = None
            if input_fits is not None:
                # Get indices from decomposition (handle both dict and object formats)
                indices_ref = None
                if isinstance(self.decomposition, dict):
                    indices_ref = self.decomposition.get('reference_spectrum_indices')
                    logger.info(f"Decomposition is dict, got reference_spectrum_indices: {type(indices_ref)}")
                elif hasattr(self.decomposition, 'reference_spectrum_indices'):
                    indices_ref = self.decomposition.reference_spectrum_indices
                    logger.info(f"Decomposition is object, got reference_spectrum_indices: {type(indices_ref)}")
                else:
                    logger.info(f"Decomposition type: {type(self.decomposition)}, has no reference_spectrum_indices")
                
                if indices_ref is not None and len(indices_ref) > 0:
                    logger.info(f"Found {len(indices_ref)} reference spectrum indices, loading...")
                    waterfall_spectra = self.load_reference_spectra_from_fits(input_fits, indices_ref)
                elif indices_ref is not None:
                    logger.warning(f"reference_spectrum_indices exists but is empty")
                else:
                    logger.info(f"reference_spectrum_indices is None")
            else:
                logger.info(f"input_fits is None")
            
            # Use reference SKYCHOPDIFF spectra if loaded, otherwise use science spectra
            if waterfall_spectra is not None:
                to_display = waterfall_spectra
                waterfall_title = f'Reference SKYCHOPDIFF Spectra\n(used for decomposition, {len(waterfall_spectra)} total)'
            else:
                to_display = original_spectra_subset
                waterfall_title = f'Science Spectra\n({len(original_spectra_subset)} spectra)'
                if input_fits is not None:
                    logger.warning(f"Reference spectra not available, displaying science spectra instead")
            
            n_display_waterfall = len(to_display)  # Display all spectra (not just first 200)
            im_waterfall = ax_waterfall.imshow(to_display[:n_display_waterfall], aspect='auto',
                                              extent=[x_axis[0], x_axis[-1], n_display_waterfall, 0],
                                              interpolation='nearest', cmap='viridis')
            ax_waterfall.set_title(waterfall_title, fontsize=10)
            ax_waterfall.set_xlabel(x_label, fontsize=9)
            ax_waterfall.set_ylabel('Spectrum #', fontsize=9)
            ax_waterfall.set_ylim(n_display_waterfall, 0)  # Top to bottom
            fig.colorbar(im_waterfall, ax=ax_waterfall, pad=0.02)
            ax_waterfall.tick_params(labelsize=8)
            
            # Calculate n_display early - all plots will use this for consistent sizing
            # Display all available processed spectra for this scan/subscan/telescope combination
            n_display = len(usable_local_indices)
            
            # Before correction heatmap with line detection contours overlay
            # IMPORTANT: Line detection was performed on sampled spectra (up to 20).
            # We show all available spectra for this scan/subscan/telescope.
            # The threshold overlay will only be drawn for the first N rows where detection was done.
            threshold_2d = None
            n_sampled_for_threshold = None
            
            if group_thresholds is not None and len(group_thresholds) > 0:
                # Try to find the threshold for THIS specific scan/subscan/telescope
                # group_thresholds keys are: (mission_id, telescope, scan, subscan)
                for group_key, thresh in group_thresholds.items():
                    group_mission_id, group_telescope, group_scan, group_subscan = group_key
                    # Match on scan, subscan, and telescope (mission_id should be same for all in this function)
                    if (group_scan == str(scan) and 
                        group_subscan == str(subscan) and 
                        group_telescope == str(telescope)):
                        threshold_2d = thresh
                        n_sampled_for_threshold = threshold_2d.shape[0]  # Number of spectra that were sampled
                        break
                
                # If no exact match found, use the first available as fallback
                if threshold_2d is None and len(group_thresholds) > 0:
                    threshold_2d = next(iter(group_thresholds.values()), None)
                    if threshold_2d is not None:
                        n_sampled_for_threshold = threshold_2d.shape[0]
            
            ax_before = fig.add_axes([2*padding + plot_width, 0.05, plot_width, 0.45])
            
            # Display original spectra as heatmap
            before_display = original_spectra_subset[:n_display].copy()
            # CRITICAL: extent MUST match the actual data shape being displayed
            # Shape of before_display: (n_display, n_channels)
            # extent format: [left, right, bottom, top] where bottom and top are in data coordinates
            im_before = ax_before.imshow(before_display, aspect='auto',
                                        extent=[x_axis[0], x_axis[-1], n_display, 0],
                                        interpolation='nearest', cmap='viridis')
            
            # Draw line detection contours as overlay if available
            if science_line_mask is not None and len(science_line_mask) > 0:
                try:
                    import cv2 as cv
                    
                    # Use actual 2D thresholds from group_thresholds if available
                    # This shows the actual line detections from OpenCV
                    if threshold_2d is not None and n_sampled_for_threshold is not None:
                        # Use threshold directly for the rows that were sampled
                        # Don't interpolate - just overlay it on the first n_sampled rows
                        threshold_2d_resized = threshold_2d[:min(n_sampled_for_threshold, n_display)].astype(bool)
                        
                        # Draw the 2D threshold mask as colored contours
                        # For each row that has any True values, draw the boundaries
                        for spectrum_idx in range(threshold_2d_resized.shape[0]):
                            row = threshold_2d_resized[spectrum_idx]
                            if np.any(row):
                                # Find contiguous regions in this row
                                diff = np.diff(row.astype(int))
                                # -1 means transition from True to False, +1 means False to True
                                starts = np.where(diff == 1)[0] + 1
                                ends = np.where(diff == -1)[0] + 1
                                
                                # Handle edge cases
                                if row[0]:
                                    starts = np.concatenate([[0], starts])
                                if row[-1]:
                                    ends = np.concatenate([ends, [len(row)]])
                                
                                # Draw vertical lines at boundaries
                                for start, end in zip(starts, ends):
                                    # Convert channel indices to frequency/velocity coordinates
                                    x_start = np.interp(start, [0, len(row)], [x_axis[0], x_axis[-1]])
                                    x_end = np.interp(end, [0, len(row)], [x_axis[0], x_axis[-1]])
                                    
                                    # Use normalized coordinates for consistency
                                    # spectrum_idx ranges from 0 to len(threshold_2d_resized)-1
                                    # We need to normalize to [0, n_display] range
                                    y_norm_start = spectrum_idx / n_display
                                    y_norm_end = (spectrum_idx + 0.9) / n_display
                                    
                                    # Draw vertical bars marking the detected line regions
                                    ax_before.axvline(x_start, ymin=y_norm_start, 
                                                    ymax=y_norm_end,
                                                    color='yellow', linewidth=1.5, alpha=0.8)
                                    ax_before.axvline(x_end, ymin=y_norm_start,
                                                    ymax=y_norm_end,
                                                    color='yellow', linewidth=1.5, alpha=0.8)
                                    # Fill the region between using NORMALIZED coordinates
                                    # Convert to 0-1 range for fill_betweenx with transform
                                    ax_before.fill_betweenx([spectrum_idx, spectrum_idx+0.9],
                                                           x_start, x_end,
                                                           color='yellow', alpha=0.15,
                                                           transform=ax_before.get_xaxis_transform())
                        has_line_contours = True
                    else:
                        has_line_contours = False
                except Exception as e:
                    logger.debug(f"Error drawing line contours: {e}")
                    has_line_contours = False
            else:
                has_line_contours = False
            
            # Add annotation showing line detection status
            title_text = f'Before Correction\n(first {n_display} specs)'
            if has_line_contours:
                line_channels = np.any(science_line_mask, axis=0)
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
            im_diff = ax_diff.imshow(difference, aspect='auto',
                                    extent=[x_axis[0], x_axis[-1], n_display, 0],
                                    interpolation='nearest', cmap='RdBu_r')
            ax_diff.set_title(f'Correction Applied\n(Original - Corrected)', fontsize=10)
            ax_diff.set_xlabel(x_label, fontsize=9)
            ax_diff.set_ylabel('Spectrum #', fontsize=9)
            ax_diff.set_ylim(n_display, 0)  # Top to bottom
            fig.colorbar(im_diff, ax=ax_diff, pad=0.02)
            ax_diff.tick_params(labelsize=8)
            
            # After correction heatmap
            ax_after = fig.add_axes([3*padding + 2*plot_width, 0.05, plot_width, 0.45])
            im_after = ax_after.imshow(corrected_spectra_subset[:n_display], aspect='auto',
                                      extent=[x_axis[0], x_axis[-1], n_display, 0],
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
            # Masked array visualization: show the spectra used for fitting with masks overlaid
            # This matches the original pca_correct.py which shows actual spectrum values
            # with masking information encoded in the data
            ax_masked = fig.add_axes([4*padding + 3*plot_width, 0.55, plot_width, 0.35])
            
            # Create a version of the spectra with masked channels replaced by NaN or zero
            # to visualize which regions are masked
            masked_spec_display = original_spectra_subset[:n_display].copy()
            
            # Mark masked regions (science lines + telluric) as NaN so they appear as white/empty
            for spec_idx in range(n_display):
                # Mark science line regions as NaN if present
                if science_line_mask is not None and spec_idx < len(science_line_mask):
                    science_regions = science_line_mask[spec_idx]
                    masked_spec_display[spec_idx, science_regions] = np.nan
                
                # Mark telluric line regions as NaN if present
                if telluric_line_mask is not None and spec_idx < len(telluric_line_mask):
                    telluric_regions = telluric_line_mask[spec_idx]
                    masked_spec_display[spec_idx, telluric_regions] = np.nan
            
            im_masked = ax_masked.imshow(masked_spec_display, aspect='auto',
                                        extent=[x_axis[0], x_axis[-1], n_display, 0],
                                        interpolation='nearest', cmap='viridis')
            ax_masked.set_title(f'Spectra Used for Fitting\n(white=masked regions)', fontsize=10)
            ax_masked.set_xlabel(x_label, fontsize=9)
            ax_masked.set_ylabel('Spectrum #', fontsize=9)
            ax_masked.set_ylim(n_display, 0)  # Top to bottom
            fig.colorbar(im_masked, ax=ax_masked, pad=0.02)
            ax_masked.tick_params(labelsize=8)
            
            # Coefficients heatmap
            ax_coeff = fig.add_axes([4*padding + 3*plot_width, 0.05, plot_width, 0.45])
            all_coeffs = []
            for idx in usable_global_indices:
                if idx in correction_details and 'coefficients' in correction_details[idx]:
                    all_coeffs.append(correction_details[idx]['coefficients'])
                else:
                    all_coeffs.append(np.zeros(len(self.components)))
            
            all_coeffs = np.array(all_coeffs)[:n_display]
            # Display with component number on X-axis, spectrum number on Y-axis
            im_coeff = ax_coeff.imshow(all_coeffs, aspect='auto',
                                      extent=[0, len(self.components), n_display, 0],
                                      interpolation='nearest', cmap='RdBu_r')
            ax_coeff.set_title(f'Coefficients\nCutoff: False', fontsize=10)
            ax_coeff.set_xlabel('Component #', fontsize=9)
            ax_coeff.set_ylabel('Spectrum #', fontsize=9)
            fig.colorbar(im_coeff, ax=ax_coeff, pad=0.02)
            ax_coeff.tick_params(labelsize=8)
            
            # ==================================================================
            # COLUMN 5: EXAMPLE SPECTRUM BREAKDOWN (RIGHT SIDE)
            # ==================================================================
            # Find example spectrum (median correction)
            correction_amounts = np.array([
                np.mean(np.abs(original_spectra_subset[i] - corrected_spectra_subset[i]))
                for i in range(len(usable_local_indices))
            ])
            
            if np.any(correction_amounts > 0):
                example_idx = np.argsort(correction_amounts)[len(correction_amounts)//2]
            else:
                example_idx = len(correction_amounts) // 2
            
            original_ex = original_spectra_subset[example_idx]
            corrected_ex = corrected_spectra_subset[example_idx]
            example_global_idx = usable_global_indices[example_idx]
            
            # Plot original spectrum
            ax_ex_orig = fig.add_axes([5*padding + 4*plot_width, 0.65, plot_width, 0.25])
            ax_ex_orig.plot(x_axis, original_ex, 'k-', lw=1.5)
            ax_ex_orig.set_title(f'Example Spectrum #{example_global_idx}\n(Median Correction)', fontsize=10)
            ax_ex_orig.grid(True, alpha=0.3)
            ax_ex_orig.tick_params(labelsize=8)
            
            # Plot correction for each component
            n_comp_show = min(3, len(self.components))
            spec_height = 0.55 / n_comp_show
            
            if example_global_idx in correction_details:
                details = correction_details[example_global_idx]
                comp_info = details.get('component_info', {})
                coeffs = details.get('coefficients', np.zeros(len(self.components)))
            else:
                comp_info = {}
                coeffs = np.zeros(len(self.components))
            
            for comp_idx in range(n_comp_show):
                ax_ex_comp = fig.add_axes([5*padding + 4*plot_width, 
                                          0.55 - (comp_idx + 1)*spec_height, 
                                          plot_width, spec_height])
                
                # Plot original and component contribution
                ax_ex_comp.plot(x_axis, original_ex, 'k-', lw=0.5, alpha=0.5, label='Orig')
                
                scaled_comp = coeffs[comp_idx] * self.components[comp_idx]
                ax_ex_comp.plot(x_axis, scaled_comp, 'b-', lw=1, label='Fitted Comp')
                
                info = comp_info.get(comp_idx, {})
                used = info.get('used', False)
                coeff_val = coeffs[comp_idx]
                noise_ratio = info.get('noise_ratio', 0)
                
                # Create label
                label_text = f"comp: {comp_idx+1}; coeff: {coeff_val:.2f}; noise_ratio: {noise_ratio:.1f}"
                if not used:
                    label_text += " (NOT USED)"
                    ax_ex_comp.text(0.03, 0.9, label_text, transform=ax_ex_comp.transAxes,
                                  fontsize=7, va='top', color='red')
                else:
                    ax_ex_comp.text(0.03, 0.9, label_text, transform=ax_ex_comp.transAxes,
                                  fontsize=7, va='top', color='green')
                
                ax_ex_comp.grid(True, alpha=0.3)
                ax_ex_comp.tick_params(labelsize=7)
                ax_ex_comp.set_ylabel(f'C{comp_idx+1}', fontsize=8)
                
                if comp_idx == n_comp_show - 1:
                    ax_ex_comp.set_xlabel(x_label, fontsize=9)
            
            fig.suptitle(f'{mission_id} | Scan {scan} | Subscan {subscan} | {telescope} | {len(usable_local_indices)} Spectra', 
                        fontsize=12, fontweight='bold')
            
            # Save plot
            plot_file = os.path.join(output_dir, 
                                    f'pca_correction_{mission_id}_scan{scan}_subscan{subscan}_{telescope}.png')
            fig.savefig(plot_file, dpi=100, bbox_inches='tight')
            plt.close(fig)
        
        logger.info(f"✓ Generated {plot_count} diagnostic plots")


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
        help='Input FITS file with spectra to correct (default: from config [output][reduced_fits])'
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
        '--overwrite',
        action='store_true',
        help='Overwrite output file if exists'
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
            config_input_file = output_config.get('reduced_fits')
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
    
    # Determine input file: explicit CLI arg > config [output][reduced_fits]
    if args.input is not None:
        input_file = args.input
        logger.info(f"✓ Input file = {input_file} (from command line)")
    elif config_input_file:
        input_file = config_input_file
        logger.info(f"✓ Input file = {input_file} (from config [output][reduced_fits])")
    else:
        logger.error("No input file specified in configuration or command line")
        logger.error("Specify via:")
        logger.error("  1. Command line: pca_correct --input /path/to/file.fits ...")
        logger.error("  2. Config [output][reduced_fits]")
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
                decomp_files = sorted(glob.glob(str(decomp_dir / "decomposition_*.pkl")))
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
        if decomposition_file:
            # Use single decomposition file
            corrector = PCACorrector(
                str(decomp_path),
                line_kernel_size=args.line_kernel_size,
                line_cutoff_std=args.line_cutoff_std,
                smoothing_kernel_size=args.smoothing_kernel
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
                    line_kernel_size=args.line_kernel_size,
                    line_cutoff_std=args.line_cutoff_std,
                    smoothing_kernel_size=args.smoothing_kernel
                )
            finally:
                # Clean up temp file
                os.unlink(temp_decomp_path)
            
            mission_decompositions_to_use = mission_decompositions
        
        logger.info("Applying PCA correction...")
        
        # Use CLI args if provided, otherwise use config values
        cutoff_variance = args.variance_cutoff if args.variance_cutoff is not None else config_variance_cutoff
        cutoff_noise_ratio = args.noise_ratio_cutoff if args.noise_ratio_cutoff is not None else config_noise_ratio_cutoff
        line_kernel_size = args.line_kernel_size if args.line_kernel_size != 51 else config_line_kernel_size
        line_cutoff_std = args.line_cutoff_std if args.line_cutoff_std != 2.0 else config_line_cutoff_std
        smoothing_kernel = args.smoothing_kernel if args.smoothing_kernel is not None else config_smoothing_kernel
        
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
            overwrite=args.overwrite,
            generate_plots=args.plot,
            output_dir=args.plot_dir,
            config_window=config_line_window if config_line_window else (tuple(args.config_window) if args.config_window else None),
            detect_science_lines=args.detect_science_lines,
            scan_filter=args.scan,
            subscan_filter=args.subscan,
            telescope_filter=args.telescope,
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
