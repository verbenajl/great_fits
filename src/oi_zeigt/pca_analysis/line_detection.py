"""
Science line detection and artifact masking.

Provides utilities for detecting the CII science line from spectral data
and creating masks for instrumental artifacts and science lines.
"""

import numpy as np
import logging
from typing import Tuple, Optional, Union
from pathlib import Path

from .errors import LineDetectionError

logger = logging.getLogger(__name__)


def detect_science_line_waterfall(spectrum_2d: np.ndarray, kernel_size: int = 5,
                                   prominence_threshold: float = 0.1) -> Tuple[int, float]:
    """
    Detect science line position from 2D spectral data using waterfall method.
    
    The waterfall method: average spectrum across spatial dimension, smooth with kernel,
    find peak. This is simple and robust for detecting bright emission lines like CII.
    
    Parameters
    ----------
    spectrum_2d : np.ndarray
        2D spectral array (shape: [n_spatial, n_channels] or [n_readouts, n_channels])
    kernel_size : int, optional
        Smoothing kernel size (default: 5)
    prominence_threshold : float, optional
        Prominence threshold relative to mean (default: 0.1)
        
    Returns
    -------
    center : int
        Detected center channel of science line
    confidence : float
        Confidence metric (0.0 to 1.0)
        
    Raises
    ------
    LineDetectionError
        If science line cannot be detected
    """
    try:
        spectrum_2d = np.asarray(spectrum_2d, dtype=float)
        
        # Average across first dimension (spatial/readouts)
        mean_spectrum = np.nanmean(spectrum_2d, axis=0)
        
        if mean_spectrum.size == 0:
            raise LineDetectionError("Empty spectrum")
        
        # Handle NaN/Inf
        valid_mask = np.isfinite(mean_spectrum)
        if np.sum(valid_mask) < kernel_size + 1:
            raise LineDetectionError("Insufficient valid channels")
        
        # Smooth spectrum
        smoothed = _smooth_spectrum(mean_spectrum, kernel_size)
        
        # Find peak
        peak_idx = np.nanargmax(smoothed)
        peak_value = smoothed[peak_idx]
        
        # Calculate confidence
        mean_val = np.nanmean(smoothed)
        std_val = np.nanstd(smoothed)
        
        if std_val == 0:
            raise LineDetectionError("Zero standard deviation in spectrum")
        
        # Normalize peak prominence
        prominence = (peak_value - mean_val) / std_val
        
        # Check if peak is significant
        if prominence < prominence_threshold:
            raise LineDetectionError(
                f"Detected peak not prominent enough (prominence={prominence:.2f})"
            )
        
        # Confidence is normalized prominence (0-1 scale)
        confidence = min(1.0, prominence / 5.0)  # 5-sigma is "perfect"
        
        return int(peak_idx), float(confidence)
    
    except LineDetectionError:
        raise
    except Exception as e:
        raise LineDetectionError(f"Science line detection failed: {e}")


def detect_science_line_max_intensity(spectrum_2d: np.ndarray) -> Tuple[int, float]:
    """
    Detect science line using maximum intensity across spatial dimension.
    
    Simpler alternative: just find channel with maximum intensity.
    
    Parameters
    ----------
    spectrum_2d : np.ndarray
        2D spectral array
        
    Returns
    -------
    center : int
        Channel of maximum intensity
    confidence : float
        Confidence metric (always 1.0 for this method)
    """
    try:
        spectrum_2d = np.asarray(spectrum_2d, dtype=float)
        
        # Maximum across all dimensions
        max_per_channel = np.nanmax(spectrum_2d, axis=0)
        
        if np.all(~np.isfinite(max_per_channel)):
            raise LineDetectionError("No valid data for line detection")
        
        center = int(np.nanargmax(max_per_channel))
        
        return center, 1.0
    
    except Exception as e:
        raise LineDetectionError(f"Max intensity detection failed: {e}")


def detect_science_line_centroid(spectrum_2d: np.ndarray, threshold: float = 0.5) -> Tuple[int, float]:
    """
    Detect science line using centroid method.
    
    Channels above threshold are used to calculate flux-weighted centroid.
    
    Parameters
    ----------
    spectrum_2d : np.ndarray
        2D spectral array
    threshold : float, optional
        Threshold as fraction of maximum (default: 0.5)
        
    Returns
    -------
    center : int
        Centroid channel position
    confidence : float
        Confidence metric
    """
    try:
        spectrum_2d = np.asarray(spectrum_2d, dtype=float)
        
        # Mean spectrum
        mean_spectrum = np.nanmean(spectrum_2d, axis=0)
        
        # Find peak
        peak_val = np.nanmax(mean_spectrum)
        peak_idx = np.nanargmax(mean_spectrum)
        
        if peak_val <= 0:
            raise LineDetectionError("No positive values in spectrum")
        
        # Create threshold mask
        threshold_val = threshold * peak_val
        above_threshold = mean_spectrum > threshold_val
        
        if np.sum(above_threshold) == 0:
            raise LineDetectionError("No channels above threshold")
        
        # Calculate flux-weighted centroid
        channels = np.arange(len(mean_spectrum))
        valid_channels = channels[above_threshold]
        valid_values = mean_spectrum[above_threshold]
        
        centroid = np.sum(valid_channels * valid_values) / np.sum(valid_values)
        
        # Confidence based on how concentrated the line is
        full_width = np.sum(above_threshold)
        confidence = min(1.0, 10.0 / full_width)  # Narrower = higher confidence
        
        return int(np.round(centroid)), float(confidence)
    
    except Exception as e:
        raise LineDetectionError(f"Centroid detection failed: {e}")


def refine_line_width(spectrum: np.ndarray, center: int,
                     fwhm_fraction: float = 0.5) -> int:
    """
    Estimate line width from spectrum using FWHM.
    
    Parameters
    ----------
    spectrum : np.ndarray
        1D spectrum
    center : int
        Approximate center channel
    fwhm_fraction : float, optional
        Fraction of peak for width estimate (default: 0.5 for FWHM)
        
    Returns
    -------
    int
        Estimated line width (full width)
    """
    try:
        spectrum = np.asarray(spectrum, dtype=float)
        
        # Get peak value
        peak_val = np.nanmax(spectrum)
        
        if peak_val <= 0:
            return 1
        
        # Find channels above threshold
        threshold = fwhm_fraction * peak_val
        above_threshold = spectrum > threshold
        
        if np.sum(above_threshold) == 0:
            return 1
        
        # Get leftmost and rightmost channels above threshold
        indices = np.where(above_threshold)[0]
        width = indices[-1] - indices[0] + 1
        
        return max(1, width)
    
    except Exception as e:
        logger.warning(f"Line width estimation failed: {e}")
        return 1


def _smooth_spectrum(spectrum: np.ndarray, kernel_size: int) -> np.ndarray:
    """
    Smooth spectrum using boxcar kernel.
    
    Parameters
    ----------
    spectrum : np.ndarray
        Input spectrum
    kernel_size : int
        Kernel size (must be odd)
        
    Returns
    -------
    np.ndarray
        Smoothed spectrum
    """
    if kernel_size < 1:
        return spectrum.copy()
    
    # Ensure kernel size is odd
    kernel_size = 2 * (kernel_size // 2) + 1
    
    kernel = np.ones(kernel_size) / kernel_size
    
    # Pad spectrum to handle edges
    pad_width = kernel_size // 2
    padded = np.pad(spectrum, pad_width, mode='edge')
    
    # Convolve
    smoothed = np.convolve(padded, kernel, mode='valid')
    
    return smoothed


class LineDetectionConfig:
    """Configuration for line detection."""
    
    def __init__(self, method: str = 'waterfall', kernel_size: int = 5,
                 prominence_threshold: float = 0.1):
        """
        Initialize line detection config.
        
        Parameters
        ----------
        method : str
            Detection method: 'waterfall', 'max_intensity', or 'centroid'
        kernel_size : int
            Smoothing kernel size
        prominence_threshold : float
            Prominence threshold for waterfall method
        """
        self.method = method
        self.kernel_size = kernel_size
        self.prominence_threshold = prominence_threshold
    
    def detect(self, spectrum_2d: np.ndarray) -> Tuple[int, float]:
        """
        Detect science line using configured method.
        
        Parameters
        ----------
        spectrum_2d : np.ndarray
            2D spectral data
            
        Returns
        -------
        center : int
            Detected center channel
        confidence : float
            Confidence metric
        """
        if self.method == 'waterfall':
            return detect_science_line_waterfall(
                spectrum_2d,
                kernel_size=self.kernel_size,
                prominence_threshold=self.prominence_threshold
            )
        elif self.method == 'max_intensity':
            return detect_science_line_max_intensity(spectrum_2d)
        elif self.method == 'centroid':
            return detect_science_line_centroid(spectrum_2d)
        else:
            raise ValueError(f"Unknown detection method: {self.method}")


def create_science_line_mask(spectrum: np.ndarray, center: int,
                            width: int) -> np.ndarray:
    """
    Create boolean mask for science line region.
    
    Parameters
    ----------
    spectrum : np.ndarray
        1D spectrum (used for shape only)
    center : int
        Center channel of science line
    width : int
        Full width of science line region
        
    Returns
    -------
    np.ndarray
        Boolean mask (True = science line, False = off-line)
    """
    n_channels = len(spectrum)
    mask = np.zeros(n_channels, dtype=bool)
    
    ch_start = max(0, center - width // 2)
    ch_end = min(n_channels, center + width // 2 + 1)
    
    mask[ch_start:ch_end] = True
    
    return mask


def create_artifact_mask(n_channels: int, artifact_center: int,
                        artifact_width: int) -> np.ndarray:
    """
    Create boolean mask for instrumental artifact region.
    
    Parameters
    ----------
    n_channels : int
        Total number of channels
    artifact_center : int
        Center channel of artifact
    artifact_width : int
        Full width of artifact region
        
    Returns
    -------
    np.ndarray
        Boolean mask (True = artifact, False = valid)
    """
    mask = np.zeros(n_channels, dtype=bool)
    
    ch_start = max(0, artifact_center - artifact_width // 2)
    ch_end = min(n_channels, artifact_center + artifact_width // 2 + 1)
    
    mask[ch_start:ch_end] = True
    
    return mask


def validate_line_detection(center: int, n_channels: int,
                           width: int, safety_margin: int = 10) -> bool:
    """
    Validate that detected line is reasonable.
    
    Parameters
    ----------
    center : int
        Detected center channel
    n_channels : int
        Total number of channels
    width : int
        Line width
    safety_margin : int
        Minimum distance from spectrum edges
        
    Returns
    -------
    bool
        True if detection is valid
    """
    # Check center is in valid range
    if center < safety_margin or center >= n_channels - safety_margin:
        return False
    
    # Check width is reasonable
    if width < 1 or width > n_channels // 2:
        return False
    
    # Check line doesn't extend beyond spectrum
    if center - width // 2 < 0 or center + width // 2 >= n_channels:
        return False
    
    return True
