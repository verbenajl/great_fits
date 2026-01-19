"""
Utility functions for PCA analysis.

Provides helper functions for baseline fitting, masking, and common operations
used throughout the PCA decomposition and correction workflows.
"""

import numpy as np
from typing import Tuple, Optional, Union
import logging
from pathlib import Path

from .errors import BaselineError

logger = logging.getLogger(__name__)


def fit_baseline_poly(spectrum: np.ndarray, order: int = 3,
                      mask: Optional[np.ndarray] = None,
                      max_iterations: int = 10,
                      threshold: float = 3.0) -> np.ndarray:
    """
    Fit polynomial baseline to spectrum using iterative sigma-clipping.
    
    Fits a polynomial to the spectrum, excluding outliers (bright lines, spikes).
    Uses iterative refinement: fit → identify outliers → refit.
    
    Parameters
    ----------
    spectrum : np.ndarray
        Input spectrum (1D array of intensity values)
    order : int, optional
        Polynomial order (default: 3)
    mask : np.ndarray, optional
        Boolean mask indicating valid channels (True = use, False = ignore)
        If None, all channels are used
    max_iterations : int, optional
        Maximum iterations for sigma-clipping refinement
    threshold : float, optional
        Sigma threshold for identifying outliers
        
    Returns
    -------
    np.ndarray
        Fitted baseline (same shape as input spectrum)
        
    Raises
    ------
    BaselineError
        If baseline fitting fails
    """
    try:
        spectrum = np.asarray(spectrum, dtype=float)
        n_channels = len(spectrum)
        
        # Create channel axis
        x = np.arange(n_channels)
        
        # Initialize mask if not provided
        if mask is None:
            mask = np.ones(n_channels, dtype=bool)
        else:
            mask = np.asarray(mask, dtype=bool)
        
        # Fit baseline iteratively
        y = spectrum.copy()
        valid = mask.copy()
        
        for iteration in range(max_iterations):
            if np.sum(valid) < order + 1:
                logger.warning(f"Not enough valid points for polynomial fit (iteration {iteration})")
                break
            
            # Fit polynomial to valid points
            coeffs = np.polyfit(x[valid], y[valid], order)
            baseline = np.polyval(coeffs, x)
            
            # Calculate residuals
            residuals = y - baseline
            sigma = np.std(residuals[valid])
            
            # Identify outliers
            outliers = np.abs(residuals) > threshold * sigma
            new_valid = valid & ~outliers
            
            # Check for convergence
            if np.sum(new_valid) == np.sum(valid):
                break
            
            valid = new_valid
        
        # Final fit with all valid points
        if np.sum(valid) >= order + 1:
            coeffs = np.polyfit(x[valid], y[valid], order)
            baseline = np.polyval(coeffs, x)
        else:
            raise BaselineError("Insufficient valid points for baseline fitting")
        
        return baseline
    
    except Exception as e:
        raise BaselineError(f"Baseline fitting failed: {e}")


def create_channel_mask(center: int, width: int, n_channels: int) -> np.ndarray:
    """
    Create a boolean mask for a channel range.
    
    Creates a mask where True indicates channels to exclude (mask region).
    
    Parameters
    ----------
    center : int
        Center channel of the mask region
    width : int
        Full width of the mask region
    n_channels : int
        Total number of channels
        
    Returns
    -------
    np.ndarray
        Boolean array (True = masked, False = valid)
    """
    mask = np.zeros(n_channels, dtype=bool)
    
    ch_start = max(0, center - width // 2)
    ch_end = min(n_channels, center + width // 2 + 1)
    
    mask[ch_start:ch_end] = True
    
    return mask


def combine_masks(*masks: np.ndarray) -> np.ndarray:
    """
    Combine multiple boolean masks with logical OR.
    
    Parameters
    ----------
    *masks : np.ndarray
        Variable number of boolean mask arrays
        
    Returns
    -------
    np.ndarray
        Combined mask (True if any input mask is True)
    """
    if not masks:
        return np.array([], dtype=bool)
    
    combined = np.zeros_like(masks[0], dtype=bool)
    for mask in masks:
        combined |= np.asarray(mask, dtype=bool)
    
    return combined


def apply_mask_to_spectrum(spectrum: np.ndarray, mask: np.ndarray,
                          value: float = 0.0) -> np.ndarray:
    """
    Apply a mask to spectrum (set masked channels to a value).
    
    Parameters
    ----------
    spectrum : np.ndarray
        Input spectrum
    mask : np.ndarray
        Boolean mask (True = mask region, False = keep)
    value : float, optional
        Value to set in masked regions (default: 0.0)
        
    Returns
    -------
    np.ndarray
        Masked spectrum (copy)
    """
    masked_spec = spectrum.copy()
    masked_spec[mask] = value
    return masked_spec


def normalize_spectrum(spectrum: np.ndarray, mean: Optional[np.ndarray] = None,
                      std: Optional[np.ndarray] = None) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Normalize spectrum by subtracting mean and dividing by std.
    
    Parameters
    ----------
    spectrum : np.ndarray
        Input spectrum
    mean : np.ndarray, optional
        Pre-calculated mean. If None, calculated from spectrum
    std : np.ndarray, optional
        Pre-calculated std. If None, calculated from spectrum
        
    Returns
    -------
    normalized : np.ndarray
        Normalized spectrum
    mean : np.ndarray
        Mean used for normalization
    std : np.ndarray
        Std used for normalization
    """
    if mean is None:
        mean = np.mean(spectrum)
    if std is None:
        std = np.std(spectrum)
    
    if std == 0:
        std = 1.0  # Avoid division by zero
    
    normalized = (spectrum - mean) / std
    
    return normalized, mean, std


def denormalize_spectrum(normalized: np.ndarray, mean: np.ndarray,
                        std: np.ndarray) -> np.ndarray:
    """
    Reverse normalization: recover original scale.
    
    Parameters
    ----------
    normalized : np.ndarray
        Normalized spectrum
    mean : np.ndarray
        Mean used for normalization
    std : np.ndarray
        Std used for normalization
        
    Returns
    -------
    np.ndarray
        Denormalized spectrum
    """
    return normalized * std + mean


def check_spectrum_validity(spectrum: np.ndarray, min_valid: int = 10) -> bool:
    """
    Check if spectrum is valid for analysis.
    
    Parameters
    ----------
    spectrum : np.ndarray
        Input spectrum
    min_valid : int, optional
        Minimum number of valid (non-NaN, non-Inf) values required
        
    Returns
    -------
    bool
        True if spectrum is valid
    """
    spectrum = np.asarray(spectrum)
    valid_count = np.sum(np.isfinite(spectrum))
    return valid_count >= min_valid


def get_spectrum_stats(spectrum: np.ndarray, mask: Optional[np.ndarray] = None) -> dict:
    """
    Calculate statistics for a spectrum.
    
    Parameters
    ----------
    spectrum : np.ndarray
        Input spectrum
    mask : np.ndarray, optional
        Boolean mask for excluded regions (True = exclude)
        
    Returns
    -------
    dict
        Dictionary with statistics: mean, std, min, max, rms
    """
    if mask is not None:
        valid_data = spectrum[~np.asarray(mask, dtype=bool)]
    else:
        valid_data = spectrum
    
    valid_data = valid_data[np.isfinite(valid_data)]
    
    return {
        'mean': np.mean(valid_data),
        'std': np.std(valid_data),
        'min': np.min(valid_data),
        'max': np.max(valid_data),
        'rms': np.sqrt(np.mean(valid_data**2)),
    }


def rolling_window(arr: np.ndarray, window_size: int) -> np.ndarray:
    """
    Create a rolling window view of an array.
    
    Parameters
    ----------
    arr : np.ndarray
        Input array (1D)
    window_size : int
        Size of the rolling window
        
    Returns
    -------
    np.ndarray
        2D array where each row is a window
    """
    arr = np.asarray(arr)
    shape = (arr.shape[0] - window_size + 1, window_size)
    strides = (arr.strides[0], arr.strides[0])
    return np.lib.stride_tricks.as_strided(arr, shape=shape, strides=strides)


def ensure_directory(path: Union[str, Path]) -> Path:
    """
    Ensure directory exists, create if necessary.
    
    Parameters
    ----------
    path : str or Path
        Directory path
        
    Returns
    -------
    Path
        Directory path (as Path object)
    """
    path = Path(path)
    path.mkdir(parents=True, exist_ok=True)
    return path


def safe_pickle_save(obj: object, filepath: Union[str, Path]) -> None:
    """
    Safely save object to pickle file.
    
    Parameters
    ----------
    obj : object
        Object to save
    filepath : str or Path
        Output pickle file path
    """
    try:
        import pickle
        filepath = Path(filepath)
        ensure_directory(filepath.parent)
        
        with open(filepath, 'wb') as f:
            pickle.dump(obj, f)
        
        logger.info(f"Saved pickle file: {filepath}")
    except Exception as e:
        logger.error(f"Failed to save pickle file: {e}")
        raise


def safe_pickle_load(filepath: Union[str, Path]) -> object:
    """
    Safely load object from pickle file.
    
    Parameters
    ----------
    filepath : str or Path
        Pickle file path
        
    Returns
    -------
    object
        Loaded object
    """
    try:
        import pickle
        filepath = Path(filepath)
        
        if not filepath.exists():
            raise FileNotFoundError(f"Pickle file not found: {filepath}")
        
        with open(filepath, 'rb') as f:
            obj = pickle.load(f)
        
        logger.info(f"Loaded pickle file: {filepath}")
        return obj
    except Exception as e:
        logger.error(f"Failed to load pickle file: {e}")
        raise
