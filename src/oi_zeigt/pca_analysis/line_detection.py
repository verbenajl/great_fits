"""
Science line detection and artifact masking.

Provides utilities for detecting the CII science line from spectral data
and creating masks for instrumental artifacts and science lines.
"""

import numpy as np
import logging
from typing import Tuple, Optional, Union
from pathlib import Path
from astropy.io import fits

from .errors import LineDetectionError

logger = logging.getLogger(__name__)


def velocity_to_channel(velocity_km_s: float, crval1: float, cdelt1: float, 
                       crpix1: float = 1.0) -> float:
    """
    Convert velocity (km/s) to channel number using WCS parameters.
    
    Uses the standard FITS WCS formula:
    channel = crpix1 - (velocity - crval1) / cdelt1
    
    Parameters
    ----------
    velocity_km_s : float
        Velocity value in km/s
    crval1 : float
        Reference velocity (CRVAL1) in km/s
    cdelt1 : float
        Velocity step per channel (CDELT1) in km/s/channel
        (usually negative for increasing velocity with decreasing channel)
    crpix1 : float, optional
        Reference pixel (CRPIX1), default=1.0 (FITS convention)
        
    Returns
    -------
    float
        Channel number (0-indexed)
        
    Raises
    ------
    ValueError
        If cdelt1 is zero or if inputs are invalid
    """
    if cdelt1 == 0:
        raise ValueError("CDELT1 (velocity step) cannot be zero")
    
    # WCS formula: channel = crpix1 - (velocity - crval1) / cdelt1
    # Subtract 1 to convert from FITS (1-indexed) to Python (0-indexed)
    channel = (crpix1 - 1.0) - (velocity_km_s - crval1) / cdelt1
    
    return channel


def velocity_range_to_channels(velocity_center_km_s: float, velocity_width_km_s: float,
                               crval1: float, cdelt1: float, crpix1: float = 1.0,
                               n_channels: Optional[int] = None) -> Tuple[int, int]:
    """
    Convert velocity range (center ± width/2) to channel range.
    
    Parameters
    ----------
    velocity_center_km_s : float
        Center velocity in km/s
    velocity_width_km_s : float
        Full width of range in km/s (will use ±width/2)
    crval1 : float
        Reference velocity (CRVAL1) in km/s
    cdelt1 : float
        Velocity step per channel (CDELT1) in km/s/channel
    crpix1 : float, optional
        Reference pixel (CRPIX1), default=1.0
    n_channels : int, optional
        Total number of channels (for bounds checking)
        
    Returns
    -------
    ch_start : int
        Start channel (clipped to valid range)
    ch_end : int
        End channel (clipped to valid range)
        
    Raises
    ------
    ValueError
        If cdelt1 is zero or if inputs are invalid
    """
    # Calculate velocity range
    v_min = velocity_center_km_s - velocity_width_km_s / 2.0
    v_max = velocity_center_km_s + velocity_width_km_s / 2.0
    
    # Convert to channels
    ch_min = velocity_to_channel(v_min, crval1, cdelt1, crpix1)
    ch_max = velocity_to_channel(v_max, crval1, cdelt1, crpix1)
    
    # Handle negative cdelt (velocity increases with decreasing channel)
    ch_start = int(np.floor(min(ch_min, ch_max)))
    ch_end = int(np.ceil(max(ch_min, ch_max))) + 1
    
    # Clip to valid range
    if n_channels is not None:
        ch_start = max(0, ch_start)
        ch_end = min(n_channels, ch_end)
    else:
        ch_start = max(0, ch_start)
        ch_end = max(ch_end, 0)
    
    return ch_start, ch_end


def extract_wcs_from_header(header: Union[fits.Header, dict]) -> Tuple[float, float, float]:
    """
    Extract WCS parameters from FITS header.
    
    Parameters
    ----------
    header : fits.Header or dict
        FITS header or dictionary-like object
        
    Returns
    -------
    crval1 : float
        Reference velocity (CRVAL1) in km/s
    cdelt1 : float
        Velocity step (CDELT1) in km/s/channel
    crpix1 : float
        Reference pixel (CRPIX1)
        
    Raises
    ------
    ValueError
        If required WCS keywords are missing
    KeyError
        If header doesn't support dictionary-like access
    """
    try:
        crval1 = float(header['CRVAL1'])
    except (KeyError, TypeError, ValueError):
        raise ValueError("CRVAL1 (reference velocity) not found in header")
    
    try:
        cdelt1 = float(header['CDELT1'])
    except (KeyError, TypeError, ValueError):
        raise ValueError("CDELT1 (velocity step) not found in header")
    
    try:
        crpix1 = float(header.get('CRPIX1', 1.0))
    except (TypeError, ValueError):
        crpix1 = 1.0
    
    if cdelt1 == 0:
        raise ValueError("CDELT1 cannot be zero")
    
    return crval1, cdelt1, crpix1


def find_channels_from_velocity_axis(velocity_center_km_s: float, velocity_width_km_s: float,
                                      velocity_axis: np.ndarray) -> Tuple[int, int]:
    """
    Find channel range for velocity region using actual velocity axis from FITS.
    
    This is more accurate than WCS conversion because it uses the actual velocity
    values stored in the FITS file, which may have been resampled or shifted.
    
    Parameters
    ----------
    velocity_center_km_s : float
        Center velocity in km/s
    velocity_width_km_s : float
        Full width of range in km/s (will use ±width/2)
    velocity_axis : np.ndarray
        1D array of velocity values for each channel
        Can be in either m/s or km/s (auto-detected based on magnitude)
        
    Returns
    -------
    ch_start : int
        Start channel (inclusive)
    ch_end : int
        End channel (exclusive)
        
    Raises
    ------
    ValueError
        If velocity_axis is invalid or telluric range doesn't overlap with data
    """
    if velocity_axis is None or len(velocity_axis) == 0:
        raise ValueError("velocity_axis is empty or None")
    
    velocity_axis = np.asarray(velocity_axis).flatten()
    
    # Auto-detect if velocity_axis is in m/s or km/s
    # If max value > 10000, assume m/s and convert to km/s
    v_max = np.nanmax(np.abs(velocity_axis))
    if v_max > 10000:  # Likely in m/s
        velocity_axis_km_s = velocity_axis / 1000.0
        logger.debug(f"VELOCITY_AXIS detected as m/s, converting to km/s (max: {v_max/1000:.1f} km/s)")
    else:
        velocity_axis_km_s = velocity_axis
        logger.debug(f"VELOCITY_AXIS detected as km/s (max: {v_max:.1f} km/s)")
    
    # Calculate velocity range
    v_min = velocity_center_km_s - velocity_width_km_s / 2.0
    v_max_range = velocity_center_km_s + velocity_width_km_s / 2.0
    
    # Find indices where velocity falls within range
    # Use searchsorted for efficient lookup
    
    # Handle both increasing and decreasing velocity arrays
    if velocity_axis_km_s[0] < velocity_axis_km_s[-1]:
        # Increasing velocity (typical)
        ch_start = np.searchsorted(velocity_axis_km_s, v_min, side='left')
        ch_end = np.searchsorted(velocity_axis_km_s, v_max_range, side='right')
    else:
        # Decreasing velocity
        ch_end_temp = np.searchsorted(velocity_axis_km_s[::-1], v_min, side='right')
        ch_start_temp = np.searchsorted(velocity_axis_km_s[::-1], v_max_range, side='left')
        ch_end = len(velocity_axis_km_s) - ch_end_temp
        ch_start = len(velocity_axis_km_s) - ch_start_temp
        ch_start, ch_end = min(ch_start, ch_end), max(ch_start, ch_end)
    
    # Ensure valid range
    ch_start = max(0, ch_start)
    ch_end = min(len(velocity_axis_km_s), ch_end)
    
    if ch_start >= ch_end:
        logger.warning(
            f"Telluric velocity range [{v_min:.1f}, {v_max_range:.1f}] km/s "
            f"does not overlap with available velocity range "
            f"[{velocity_axis_km_s.min():.1f}, {velocity_axis_km_s.max():.1f}] km/s"
        )
    
    return ch_start, ch_end


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
                        artifact_width: int, 
                        header: Optional[Union[fits.Header, dict]] = None,
                        velocity_axis: Optional[np.ndarray] = None,
                        use_velocity: bool = True) -> np.ndarray:
    """
    Create boolean mask for instrumental artifact region (e.g., telluric lines).
    
    Can work in multiple modes, with preference given to:
    1. VELOCITY_AXIS column (if provided) - most accurate, uses actual velocity data
    2. WCS parameters (CRVAL1, CDELT1, CRPIX1 from header) - standard conversion
    3. Channel mode (default) - artifact_center and artifact_width in channels
    
    Parameters
    ----------
    n_channels : int
        Total number of channels
    artifact_center : float
        Center of artifact region (in channels or km/s depending on use_velocity)
    artifact_width : float
        Full width of artifact region (in channels or km/s depending on use_velocity)
    header : fits.Header or dict, optional
        FITS header with WCS information (CRVAL1, CDELT1, CRPIX1)
        If provided with use_velocity=True, input values are treated as km/s
    velocity_axis : np.ndarray, optional
        1D array of velocity values (from VELOCITY_AXIS column in FITS file)
        If provided, this is used preferentially over WCS conversion
    use_velocity : bool, optional
        If True and header/velocity_axis is provided, treat artifact_center/width as km/s
        Default: True
        
    Returns
    -------
    np.ndarray
        Boolean mask (True = artifact, False = valid)
    """
    mask = np.zeros(n_channels, dtype=bool)
    ch_start = None
    ch_end = None
    
    # Method 1: Use VELOCITY_AXIS from FITS if available (most accurate)
    if velocity_axis is not None and use_velocity:
        try:
            velocity_axis = np.asarray(velocity_axis).flatten()
            
            # Find channels corresponding to telluric velocity range
            ch_start, ch_end = find_channels_from_velocity_axis(
                velocity_center_km_s=artifact_center,
                velocity_width_km_s=artifact_width,
                velocity_axis=velocity_axis
            )
            
            logger.info(
                f"✓ Telluric mask: {artifact_center:.1f}±{artifact_width/2:.1f} km/s "
                f"→ channels {ch_start}-{ch_end} "
                f"(from VELOCITY_AXIS column, {len(velocity_axis)} channels)"
            )
        
        except (ValueError, IndexError, TypeError) as e:
            logger.warning(
                f"Could not use VELOCITY_AXIS: {e}. "
                f"Falling back to WCS conversion..."
            )
            velocity_axis = None  # Fall through to WCS method
    
    # Method 2: Use WCS parameters from header if VELOCITY_AXIS not available
    if ch_start is None and velocity_axis is None and header is not None and use_velocity:
        try:
            # Extract WCS parameters from header
            crval1, cdelt1, crpix1 = extract_wcs_from_header(header)
            
            # Convert velocity range to channel range
            ch_start, ch_end = velocity_range_to_channels(
                velocity_center_km_s=artifact_center,
                velocity_width_km_s=artifact_width,
                crval1=crval1,
                cdelt1=cdelt1,
                crpix1=crpix1,
                n_channels=n_channels
            )
            
            logger.info(
                f"Telluric mask: velocity {artifact_center:.1f}±{artifact_width/2:.1f} km/s "
                f"→ channels {ch_start}-{ch_end} "
                f"(using WCS: CRVAL1={crval1:.1f}, CDELT1={cdelt1:.3f}, CRPIX1={crpix1:.1f})"
            )
        
        except (ValueError, KeyError, TypeError) as e:
            logger.warning(
                f"Could not convert velocity to channels: {e}. "
                f"Falling back to channel-based interpretation."
            )
            ch_start = None  # Fall through to channel-based mode
    
    # Method 3: Channel-based mode (fallback if both velocity methods failed)
    if ch_start is None:
        ch_start = max(0, int(artifact_center - artifact_width / 2))
        ch_end = min(n_channels, int(artifact_center + artifact_width / 2) + 1)
    
    if ch_start is not None and ch_end is not None:
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
