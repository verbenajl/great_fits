"""
Quality analysis functions for OI ZEIGT spectral data.

Provides tools for analyzing spectral quality metrics such as RMS ratios,
signal-to-noise ratios, and other quality indicators.
"""

import warnings
from typing import Optional, Union, Tuple, List, Dict
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits


def _get_column_data(hdul: fits.HDUList,
                     column_name: str,
                     object_filter: Optional[str] = None) -> np.ndarray:
    """
    Extract column data from FITS file with optional object filtering.
    
    Parameters
    ----------
    hdul : astropy.io.fits.HDUList
        FITS HDUList containing the data.
    column_name : str
        Name of the column to extract.
    object_filter : str, optional
        Filter by object name (substring match, case-insensitive).
    
    Returns
    -------
    np.ndarray
        Extracted column data (with NaN filtering applied).
        For 1D columns, returns 1D array with NaN values removed.
        For 2D columns (e.g., ROLL_RMS_N), returns 2D array as-is.
    
    Raises
    ------
    ValueError
        If column not found or no data matches filter.
    """
    # Find the binary table HDU
    matrix_hdu = None
    for idx, hdu in enumerate(hdul):
        if hasattr(hdu, 'data') and hdu.data is not None:
            if column_name in hdu.data.dtype.names:
                matrix_hdu = hdu
                break

    if matrix_hdu is None:
        raise ValueError(f"Column '{column_name}' not found in FITS file")

    data = matrix_hdu.data

    # Filter by object if specified
    if object_filter:
        if 'OBJECT' not in data.dtype.names:
            raise ValueError("OBJECT column not found in FITS file")
        
        objects = data['OBJECT']
        object_mask = np.array([
            object_filter.upper() in (obj.decode().strip() if isinstance(obj, bytes) else str(obj).strip()).upper()
            for obj in objects
        ])
        
        if not np.any(object_mask):
            raise ValueError(f"No objects matching '{object_filter}' found in FITS file")
        
        col_data = data[column_name][object_mask]
    else:
        col_data = data[column_name]

    # Remove NaN values (only for 1D arrays - 2D arrays are handled specially elsewhere)
    if col_data.ndim == 1:
        col_data = col_data[np.isfinite(col_data)]
    
    return col_data


def get_spechistogram(hdul: fits.HDUList,
                      metrics: List[str],
                      object_filter: Optional[str] = None,
                      bins: int = 30,
                      figsize: Optional[Tuple[int, int]] = None) -> Tuple[plt.Figure, Dict]:
    """
    Generate histograms of multiple spectral quality metrics in a single figure.

    Parameters
    ----------
    hdul : astropy.io.fits.HDUList
        FITS HDUList containing spectral data.
    metrics : list of str
        List of metric names to plot. Available metrics:
        - 'rmsratio': RMS ratio quality metric
        - 'squality': Signal quality flag
        - 'roll_rms_n': RMS rolloff (uses mean across channels)
        - 'mh2o': Water vapor column density
        - 'tsys': System temperature
        - 'tau_atm': Atmospheric optical depth
        - 'chi_sqr': Chi-square fit value
        - 'err_pwv': PWV error
        - 'rms_baseline': RMS baseline metric
    object_filter : str, optional
        Filter spectra by object name (substring match, case-insensitive).
        If None, all spectra are included.
    bins : int, optional
        Number of histogram bins (default: 30).
    figsize : tuple, optional
        Figure size as (width, height). If None, auto-sized based on number of metrics.

    Returns
    -------
    fig : matplotlib.figure.Figure
        The generated figure.
    stats_dict : dict
        Dictionary with statistics for each metric.

    Raises
    ------
    ValueError
        If required columns are not found or invalid metric names provided.
    """
    # Map metric names to FITS column names
    metric_mapping = {
        'rmsratio': ('RMSRATIO', 'RMS Ratio'),
        'rmsratiob': ('RMSRATIOB', 'RMS Ratio (radiometer)'),
        'squality': ('SQUALITY', 'Signal Quality'),
        'roll_rms_n': ('ROLL_RMS_N', 'Roll RMS'),
        'mh2o': ('MH2O', 'Water Vapor (MH₂O)'),
        'tsys': ('TSYS', 'System Temperature (Tsys)'),
        'tau_atm': ('TAU-ATM', 'Atmospheric Optical Depth (τ)'),
        'chi_sqr': ('CHI_SQR', 'Chi-square'),
        'err_pwv': ('ERR_PWV', 'PWV Error'),
        'rms_baseline': ('RMS_BASELINE', 'RMS Baseline'),
    }

    # Validate metrics
    invalid_metrics = [m for m in metrics if m.lower() not in metric_mapping]
    if invalid_metrics:
        raise ValueError(f"Invalid metrics: {invalid_metrics}. Valid options: {list(metric_mapping.keys())}")

    # Auto-size figure if not specified
    if figsize is None:
        n_metrics = len(metrics)
        cols = min(3, n_metrics)
        rows = (n_metrics + cols - 1) // cols
        figsize = (6 * cols, 5 * rows)

    # Create subplots
    n_metrics = len(metrics)
    cols = min(3, n_metrics)
    rows = (n_metrics + cols - 1) // cols
    
    fig, axes = plt.subplots(rows, cols, figsize=figsize)
    if n_metrics == 1:
        axes = np.array([axes])
    axes = axes.flatten()

    stats_dict = {}

    # Generate histogram for each metric
    for idx, metric in enumerate(metrics):
        metric_lower = metric.lower()
        if metric_lower not in metric_mapping:
            continue
        
        col_name, label = metric_mapping[metric_lower]
        
        try:
            # Get data for this metric
            col_data = _get_column_data(hdul, col_name, object_filter)
            
            # Handle 2D arrays (like ROLL_RMS_N which has one value per channel)
            if col_data.ndim == 2:
                # For 2D arrays, compute the mean across channels (axis=1)
                col_data = np.nanmean(col_data, axis=1)
            
            # Remove any remaining NaN values
            col_data = col_data[np.isfinite(col_data)]
            
            # Compute statistics
            stats = {
                'count': len(col_data),
                'mean': np.mean(col_data),
                'median': np.median(col_data),
                'std': np.std(col_data),
                'min': np.min(col_data),
                'max': np.max(col_data),
                'percentile_5': np.percentile(col_data, 5),
                'percentile_95': np.percentile(col_data, 95),
            }
            stats_dict[metric_lower] = stats
            
            ax = axes[idx]
            
            # Plot histogram
            ax.hist(col_data, bins=bins, edgecolor='black', alpha=0.7, color='steelblue')
            
            # Add mean and median lines
            mean_val = stats['mean']
            median_val = stats['median']
            ax.axvline(mean_val, color='red', linestyle='--', linewidth=2, label=f'μ={mean_val:.3f}')
            ax.axvline(median_val, color='green', linestyle='--', linewidth=2, label=f'M={median_val:.3f}')
            
            ax.set_xlabel(label, fontsize=11)
            ax.set_ylabel('Count', fontsize=11)
            ax.set_title(f'{label}\n(n={stats["count"]}, σ={stats["std"]:.3f})', fontsize=12)
            ax.legend(loc='best', fontsize=9)
            ax.grid(True, alpha=0.3)
        
        except ValueError as e:
            ax = axes[idx]
            ax.text(0.5, 0.5, f'Error:\n{str(e)}', ha='center', va='center', 
                   transform=ax.transAxes, fontsize=10, color='red')
            ax.set_xticks([])
            ax.set_yticks([])

    # Hide unused subplots
    for idx in range(n_metrics, len(axes)):
        axes[idx].set_visible(False)

    # Add overall title
    obj_label = f" ({object_filter})" if object_filter else ""
    fig.suptitle(f'Spectral Quality Metrics{obj_label}', fontsize=14, fontweight='bold')

    plt.tight_layout(pad=2.0, h_pad=8.0, w_pad=2.5, rect=[0, 0, 1, 0.96])

    return fig, stats_dict


def get_rmsratio_histogram(hdul: fits.HDUList,
                           object_filter: Optional[str] = None,
                           spectrum_column: str = 'SPECTRUM',
                           rmsratio_column: str = 'RMSRATIO',
                           bins: int = 30,
                           figsize: Tuple[int, int] = (12, 6)) -> Tuple[plt.Figure, np.ndarray, np.ndarray]:
    """
    Generate a histogram of RMSRATIO values for spectra.

    Parameters
    ----------
    hdul : astropy.io.fits.HDUList
        FITS HDUList containing spectral data.
    object_filter : str, optional
        Filter spectra by object name (substring match, case-insensitive).
        If None, all spectra are included.
    spectrum_column : str, optional
        Name of the spectrum column (default: 'SPECTRUM').
    rmsratio_column : str, optional
        Name of the RMS ratio column (default: 'RMSRATIO').
    bins : int, optional
        Number of histogram bins (default: 30).
    figsize : tuple, optional
        Figure size as (width, height) in inches (default: (12, 6)).

    Returns
    -------
    fig : matplotlib.figure.Figure
        The generated figure.
    hist : np.ndarray
        Histogram values.
    bin_edges : np.ndarray
        Bin edges.

    Raises
    ------
    ValueError
        If required columns are not found in the FITS file.
    """
    # Find the binary table HDU with spectra
    matrix_hdu = None
    for idx, hdu in enumerate(hdul):
        if hasattr(hdu, 'data') and hdu.data is not None:
            if spectrum_column in hdu.data.dtype.names:
                matrix_hdu = hdu
                break

    if matrix_hdu is None:
        raise ValueError(f"No HDU with '{spectrum_column}' column found")

    if rmsratio_column not in matrix_hdu.data.dtype.names:
        raise ValueError(f"Column '{rmsratio_column}' not found in FITS file")

    data = matrix_hdu.data

    # Filter by object if specified
    if object_filter:
        if 'OBJECT' not in data.dtype.names:
            raise ValueError("OBJECT column not found in FITS file")
        
        objects = data['OBJECT']
        object_mask = np.array([
            object_filter.upper() in (obj.decode().strip() if isinstance(obj, bytes) else str(obj).strip()).upper()
            for obj in objects
        ])
        
        if not np.any(object_mask):
            raise ValueError(f"No objects matching '{object_filter}' found in FITS file")
        
        rmsratio_values = data[rmsratio_column][object_mask]
        filtered_count = np.sum(object_mask)
    else:
        rmsratio_values = data[rmsratio_column]
        filtered_count = len(data)

    # Remove NaN values
    rmsratio_values = rmsratio_values[np.isfinite(rmsratio_values)]

    # Create histogram
    fig, ax = plt.subplots(figsize=figsize)
    hist, bin_edges, patches = ax.hist(rmsratio_values, bins=bins, edgecolor='black', alpha=0.7)

    # Add statistics
    mean_val = np.mean(rmsratio_values)
    median_val = np.median(rmsratio_values)
    std_val = np.std(rmsratio_values)

    # Plot mean and median lines
    ax.axvline(mean_val, color='red', linestyle='--', linewidth=2, label=f'Mean: {mean_val:.3f}')
    ax.axvline(median_val, color='green', linestyle='--', linewidth=2, label=f'Median: {median_val:.3f}')

    # Labels and title
    obj_label = f" ({object_filter})" if object_filter else ""
    ax.set_xlabel('RMSRATIO', fontsize=12)
    ax.set_ylabel('Count', fontsize=12)
    ax.set_title(f'RMSRATIO Distribution{obj_label}\n(n={len(rmsratio_values)}, σ={std_val:.3f})', fontsize=14)
    ax.legend(loc='best', fontsize=11)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()

    return fig, hist, bin_edges


def rmsratio_statistics(hdul: fits.HDUList,
                       object_filter: Optional[str] = None,
                       rmsratio_column: str = 'RMSRATIO') -> dict:
    """
    Compute statistics of RMSRATIO values for spectra.

    Parameters
    ----------
    hdul : astropy.io.fits.HDUList
        FITS HDUList containing spectral data.
    object_filter : str, optional
        Filter spectra by object name (substring match, case-insensitive).
        If None, all spectra are included.
    rmsratio_column : str, optional
        Name of the RMS ratio column (default: 'RMSRATIO').

    Returns
    -------
    dict
        Dictionary containing statistics:
        - count: Number of spectra
        - mean: Mean RMS ratio
        - median: Median RMS ratio
        - std: Standard deviation
        - min: Minimum value
        - max: Maximum value
        - percentile_5: 5th percentile
        - percentile_95: 95th percentile
    """
    # Find the binary table HDU with spectra
    matrix_hdu = None
    for idx, hdu in enumerate(hdul):
        if hasattr(hdu, 'data') and hdu.data is not None:
            if rmsratio_column in hdu.data.dtype.names:
                matrix_hdu = hdu
                break

    if matrix_hdu is None:
        raise ValueError(f"Column '{rmsratio_column}' not found in FITS file")

    data = matrix_hdu.data

    # Filter by object if specified
    if object_filter:
        if 'OBJECT' not in data.dtype.names:
            raise ValueError("OBJECT column not found in FITS file")
        
        objects = data['OBJECT']
        object_mask = np.array([
            object_filter.upper() in (obj.decode().strip() if isinstance(obj, bytes) else str(obj).strip()).upper()
            for obj in objects
        ])
        
        if not np.any(object_mask):
            raise ValueError(f"No objects matching '{object_filter}' found in FITS file")
        
        rmsratio_values = data[rmsratio_column][object_mask]
    else:
        rmsratio_values = data[rmsratio_column]

    # Remove NaN values
    rmsratio_values = rmsratio_values[np.isfinite(rmsratio_values)]

    return {
        'count': len(rmsratio_values),
        'mean': np.mean(rmsratio_values),
        'median': np.median(rmsratio_values),
        'std': np.std(rmsratio_values),
        'min': np.min(rmsratio_values),
        'max': np.max(rmsratio_values),
        'percentile_5': np.percentile(rmsratio_values, 5),
        'percentile_95': np.percentile(rmsratio_values, 95),
    }


def ripple_ratio(spectra: np.ndarray, linefree_mask: np.ndarray) -> np.ndarray:
    """
    Baseline-ripple statistic per spectrum: total RMS / white-noise RMS.

    The white-noise RMS comes from channel-to-channel differences,
    std(diff(x)) / sqrt(2), which slow structure (ripples, standing waves,
    curvature, band-edge droop) barely changes. The total RMS includes that
    structure. For pure white noise the ratio is ~1; correlated baseline
    structure pushes it above 1. Unlike the Allan WHITENESS test, the scale is
    absolute: it does not depend on the rest of the file.

    Only line-free channels are used, and differences are never taken across a
    gap in `linefree_mask` (e.g. across the line window). Each contiguous
    segment is mean-subtracted separately, so an offset between the two sides
    of the line window does not count as ripple. NaN channels are dropped.

    Smoothed spectra are NOT white (neighbouring channels are correlated), so
    the ratio is only meaningful on unsmoothed spectra.

    Parameters
    ----------
    spectra : np.ndarray
        2D array (n_spectra, n_channels).
    linefree_mask : np.ndarray of bool
        1D (n_channels,) mask, True for the channels to use.

    Returns
    -------
    np.ndarray
        Ratio per spectrum (float64); NaN where fewer than 8 usable channels.
    """
    spectra = np.asarray(spectra, dtype=np.float64)
    linefree_mask = np.asarray(linefree_mask, dtype=bool)
    # Contiguous runs of line-free channels
    edges = np.diff(np.r_[0, linefree_mask.astype(np.int8), 0])
    starts, stops = np.where(edges == 1)[0], np.where(edges == -1)[0]

    sum_sq = np.zeros(len(spectra))
    n_tot = np.zeros(len(spectra))
    diffs = []
    for a, b in zip(starts, stops):
        if b - a < 2:
            continue
        seg = spectra[:, a:b]
        with np.errstate(invalid='ignore'), warnings.catch_warnings():
            warnings.simplefilter('ignore', RuntimeWarning)
            dev = seg - np.nanmean(seg, axis=1, keepdims=True)
        sum_sq += np.nansum(dev ** 2, axis=1)
        n_tot += np.sum(np.isfinite(seg), axis=1)
        diffs.append(np.diff(seg, axis=1))   # NaN wherever either channel is NaN
    if not diffs:
        return np.full(len(spectra), np.nan)

    d = np.concatenate(diffs, axis=1)
    n_d = np.sum(np.isfinite(d), axis=1)
    with np.errstate(invalid='ignore', divide='ignore'), warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        sigma_total = np.sqrt(sum_sq / n_tot)
        sigma_white = np.nanstd(d, axis=1) / np.sqrt(2.0)
        ratio = sigma_total / sigma_white
    ratio[(n_tot < 8) | (n_d < 8) | ~(sigma_white > 0)] = np.nan
    return ratio
