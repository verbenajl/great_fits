"""
Core data reduction functions for OI GREAT data.
"""

from typing import Optional, Union, Tuple
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.io import fits


def detect_blank_channels(spectrum: np.ndarray, blank_value: float = 0.0, 
                         threshold: float = 1e-10) -> Tuple[np.ndarray, float]:
    """
    Detect blank/missing channels in a spectrum.
    
    Parameters
    ----------
    spectrum : np.ndarray
        Input spectrum array.
    blank_value : float, optional
        Value considered as blank (default: 0.0).
    threshold : float, optional
        Tolerance for considering a value as blank (default: 1e-10).
    
    Returns
    -------
    tuple
        - blank_mask : np.ndarray (bool)
            Boolean array indicating blank channels (True = blank).
        - fraction_blank : float
            Fraction of blank channels (0.0 to 1.0).
    
    Examples
    --------
    >>> mask, frac = detect_blank_channels(spectrum)
    >>> print(f"Blank fraction: {frac:.2%}")
    """
    blank_mask = np.abs(spectrum - blank_value) < threshold
    fraction_blank = np.sum(blank_mask) / len(spectrum)
    
    return blank_mask, fraction_blank


def detect_nan_channels(spectrum: np.ndarray) -> Tuple[np.ndarray, float]:
    """
    Detect NaN (missing) channels in a spectrum.
    
    Parameters
    ----------
    spectrum : np.ndarray
        Input spectrum array.
    
    Returns
    -------
    tuple
        - nan_mask : np.ndarray (bool)
            Boolean array indicating NaN channels (True = NaN).
        - fraction_nan : float
            Fraction of NaN channels (0.0 to 1.0).
    
    Examples
    --------
    >>> mask, frac = detect_nan_channels(spectrum)
    >>> print(f"NaN fraction: {frac:.2%}")
    """
    nan_mask = np.isnan(spectrum)
    fraction_nan = np.sum(nan_mask) / len(spectrum)
    
    return nan_mask, fraction_nan


def split_spectra_by_blanks(hdul: fits.HDUList, 
                            blank_threshold: float = 0.20,
                            blank_value: float = 0.0,
                            blank_tolerance: float = 1e-10) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Split spectra into two dataframes based on blank channel fraction.
    
    Parameters
    ----------
    hdul : astropy.io.fits.HDUList
        FITS HDU list containing the data.
    blank_threshold : float, optional
        Threshold fraction for blank channels (default: 0.20 = 20%).
        Spectra with blank fraction > threshold go to "with_blanks" dataframe.
    blank_value : float, optional
        Value considered as blank (default: 0.0).
    blank_tolerance : float, optional
        Tolerance for considering a value as blank (default: 1e-10).
    
    Returns
    -------
    tuple
        - df_clean : pd.DataFrame
            Dataframe with spectra having <= blank_threshold fraction of blanks.
        - df_with_blanks : pd.DataFrame
            Dataframe with spectra having > blank_threshold fraction of blanks.
    
    Examples
    --------
    >>> hdul = fits.open('file.fits')
    >>> df_clean, df_blanks = split_spectra_by_blanks(hdul, blank_threshold=0.20)
    >>> print(f"Clean: {len(df_clean)}, With blanks: {len(df_with_blanks)}")
    """
    # Find the binary table HDU with spectra
    matrix_hdu = None
    for hdu in hdul:
        if hasattr(hdu, 'data') and hdu.data is not None:
            if 'SPECTRUM' in hdu.data.dtype.names:
                matrix_hdu = hdu
                break
    
    if matrix_hdu is None:
        raise ValueError("No HDU with SPECTRUM column found")
    
    data = matrix_hdu.data
    
    # Check for blank channels in each spectrum
    blank_fractions = []
    blank_masks_list = []
    
    for spectrum in data['SPECTRUM']:
        mask, frac = detect_blank_channels(spectrum, blank_value, blank_tolerance)
        blank_fractions.append(frac)
        blank_masks_list.append(mask)
    
    blank_fractions = np.array(blank_fractions)
    
    # Create mask for clean vs blanks
    clean_mask = blank_fractions <= blank_threshold
    blanks_mask = ~clean_mask
    
    # Convert to dataframe for easier handling
    df_full = pd.DataFrame(data)
    
    # Add blank fraction as a new column
    df_full['blank_fraction'] = blank_fractions
    
    # Split into two dataframes
    df_clean = df_full[clean_mask].reset_index(drop=True)
    df_with_blanks = df_full[blanks_mask].reset_index(drop=True)
    
    return df_clean, df_with_blanks


def split_spectra_by_nans(hdul: fits.HDUList, 
                         nan_threshold: float = 0.80) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Split spectra into two dataframes based on NaN fraction.
    
    Parameters
    ----------
    hdul : astropy.io.fits.HDUList
        FITS HDU list containing the data.
    nan_threshold : float, optional
        Threshold fraction for NaN channels (default: 0.80 = 80%).
        Spectra with NaN fraction >= threshold go to "with_many_nans" dataframe.
    
    Returns
    -------
    tuple
        - df_clean : pd.DataFrame
            Dataframe with spectra having < nan_threshold fraction of NaNs.
        - df_many_nans : pd.DataFrame
            Dataframe with spectra having >= nan_threshold fraction of NaNs.
    
    Examples
    --------
    >>> hdul = fits.open('file.fits')
    >>> df_clean, df_many_nans = split_spectra_by_nans(hdul, nan_threshold=0.80)
    >>> print(f"Clean: {len(df_clean)}, Many NaNs: {len(df_many_nans)}")
    """
    # Find the binary table HDU with spectra
    matrix_hdu = None
    for hdu in hdul:
        if hasattr(hdu, 'data') and hdu.data is not None:
            if 'SPECTRUM' in hdu.data.dtype.names:
                matrix_hdu = hdu
                break
    
    if matrix_hdu is None:
        raise ValueError("No HDU with SPECTRUM column found")
    
    data = matrix_hdu.data
    
    # Check for NaN channels in each spectrum
    nan_fractions = []
    
    for spectrum in data['SPECTRUM']:
        _, frac = detect_nan_channels(spectrum)
        nan_fractions.append(frac)
    
    nan_fractions = np.array(nan_fractions)
    
    # Create mask for clean vs many nans
    clean_mask = nan_fractions < nan_threshold
    many_nans_mask = ~clean_mask
    
    # Convert to dataframe for easier handling
    df_full = pd.DataFrame(data)
    
    # Add NaN fraction as a new column
    df_full['nan_fraction'] = nan_fractions
    
    # Split into two dataframes
    df_clean = df_full[clean_mask].reset_index(drop=True)
    df_many_nans = df_full[many_nans_mask].reset_index(drop=True)
    
    return df_clean, df_many_nans


def get_blank_statistics(hdul: fits.HDUList, 
                         blank_value: float = 0.0,
                         blank_tolerance: float = 1e-10) -> dict:
    """
    Get statistics about blank channels in the FITS file.
    
    Parameters
    ----------
    hdul : astropy.io.fits.HDUList
        FITS HDU list containing the data.
    blank_value : float, optional
        Value considered as blank (default: 0.0).
    blank_tolerance : float, optional
        Tolerance for considering a value as blank (default: 1e-10).
    
    Returns
    -------
    dict
        Statistics including:
        - total_spectra: Total number of spectra
        - mean_blank_fraction: Mean fraction of blanks across all spectra
        - min_blank_fraction: Minimum blank fraction
        - max_blank_fraction: Maximum blank fraction
        - spectra_with_no_blanks: Number of spectra with no blanks
        - spectra_with_some_blanks: Number of spectra with some blanks
    
    Examples
    --------
    >>> stats = get_blank_statistics(hdul)
    >>> print(f"Mean blank fraction: {stats['mean_blank_fraction']:.2%}")
    """
    # Find the binary table HDU with spectra
    matrix_hdu = None
    for hdu in hdul:
        if hasattr(hdu, 'data') and hdu.data is not None:
            if 'SPECTRUM' in hdu.data.dtype.names:
                matrix_hdu = hdu
                break
    
    if matrix_hdu is None:
        raise ValueError("No HDU with SPECTRUM column found")
    
    data = matrix_hdu.data
    
    # Calculate blank fractions
    blank_fractions = []
    for spectrum in data['SPECTRUM']:
        _, frac = detect_blank_channels(spectrum, blank_value, blank_tolerance)
        blank_fractions.append(frac)
    
    blank_fractions = np.array(blank_fractions)
    
    stats = {
        'total_spectra': len(data),
        'mean_blank_fraction': np.mean(blank_fractions),
        'median_blank_fraction': np.median(blank_fractions),
        'std_blank_fraction': np.std(blank_fractions),
        'min_blank_fraction': np.min(blank_fractions),
        'max_blank_fraction': np.max(blank_fractions),
        'spectra_with_no_blanks': np.sum(blank_fractions == 0.0),
        'spectra_with_some_blanks': np.sum(blank_fractions > 0.0),
    }
    
    return stats


def analyze_spectrum_values(hdul: fits.HDUList, sample_size: int = 100) -> dict:
    """
    Analyze the distribution of values in spectra to identify blank values.
    
    This function samples spectra and shows you the most common values,
    which helps identify what might be used as a blank/missing value marker.
    
    Parameters
    ----------
    hdul : astropy.io.fits.HDUList
        FITS HDU list containing the data.
    sample_size : int, optional
        Number of spectra to sample for analysis (default: 100).
    
    Returns
    -------
    dict
        Analysis results including:
        - common_values: List of (value, count) tuples for most common values
        - min_value: Minimum value across all sampled spectra
        - max_value: Maximum value across all sampled spectra
        - mean_value: Mean value across all sampled spectra
        - std_value: Standard deviation across all sampled spectra
        - zero_count: Number of channels with value 0.0
        - negative_count: Number of negative values
        - very_small_count: Number of values < 1e-10
    
    Examples
    --------
    >>> analysis = analyze_spectrum_values(hdul, sample_size=50)
    >>> for value, count in analysis['common_values'][:10]:
    ...     print(f"Value {value}: appears {count} times")
    """
    # Find the binary table HDU with spectra
    matrix_hdu = None
    for hdu in hdul:
        if hasattr(hdu, 'data') and hdu.data is not None:
            if 'SPECTRUM' in hdu.data.dtype.names:
                matrix_hdu = hdu
                break
    
    if matrix_hdu is None:
        raise ValueError("No HDU with SPECTRUM column found")
    
    data = matrix_hdu.data
    
    # Sample spectra
    n_spectra = len(data)
    sample_indices = np.random.choice(n_spectra, size=min(sample_size, n_spectra), replace=False)
    
    all_values = []
    for idx in sample_indices:
        spectrum = data['SPECTRUM'][idx]
        all_values.extend(spectrum)
    
    all_values = np.array(all_values)
    
    # Count occurrences of each value
    from collections import Counter
    value_counts = Counter(all_values)
    
    # Sort by frequency
    common_values = sorted(value_counts.items(), key=lambda x: x[1], reverse=True)
    
    analysis = {
        'common_values': common_values[:20],  # Top 20 most common
        'min_value': float(np.min(all_values)),
        'max_value': float(np.max(all_values)),
        'mean_value': float(np.mean(all_values)),
        'median_value': float(np.median(all_values)),
        'std_value': float(np.std(all_values)),
        'zero_count': np.sum(all_values == 0.0),
        'negative_count': np.sum(all_values < 0.0),
        'very_small_count': np.sum(np.abs(all_values) < 1e-10),
        'total_values_analyzed': len(all_values),
    }
    
    return analysis


def filter_and_save_fits(hdul: fits.HDUList, 
                        object_name: str,
                        nan_threshold: float = 0.20,
                        output_clean: Optional[Union[str, Path]] = None,
                        output_rejected: Optional[Union[str, Path]] = None,
                        remove_column: Optional[str] = None,
                        remove_values: Optional[list] = None) -> Tuple[Path, Path]:
    """
    Filter FITS data by object and NaN content, with optional column value removal.
    
    This function:
    1. Processes ALL HDUs with SPECTRUM columns (important for combined files)
    2. Separates data into target object and other objects
    3. Filters target object spectra to keep only those with < nan_threshold NaN channels
    4. Optionally removes rows where a specified column matches given values
    5. Saves two FITS files:
       - output_clean: All non-target objects + filtered target object spectra
       - output_rejected: Target object spectra that were filtered out or removed
    
    Parameters
    ----------
    hdul : astropy.io.fits.HDUList
        Input FITS HDU list containing the data.
    object_name : str
        Name of the object to filter (e.g., "M51").
        Spectra with this string in OBJECT column will be filtered by NaN content.
    nan_threshold : float, optional
        Maximum acceptable NaN fraction (default: 0.20 = 20%).
        Keep spectra with NaN fraction < threshold.
        Reject spectra with NaN fraction >= threshold.
    output_clean : str or Path, optional
        Output path for clean FITS file (non-target + filtered target objects).
        If None, uses "clean_data.fits".
    output_rejected : str or Path, optional
        Output path for rejected FITS file (filtered out target objects).
        If None, uses "rejected_data.fits".
    remove_column : str, optional
        Column name to check for removal (e.g., "AOR_ID").
        If specified, rows matching remove_values will be removed.
    remove_values : list, optional
        List of values to remove from remove_column.
        Rows with these values will be moved to rejected file.
    
    Returns
    -------
    tuple
        - clean_path : Path
            Path to the clean FITS file
        - rejected_path : Path
            Path to the rejected FITS file
    
    Examples
    --------
    >>> hdul = fits.open('original.fits')
    >>> clean_path, rejected_path = filter_and_save_fits(
    ...     hdul, 
    ...     object_name="M51",
    ...     nan_threshold=0.20,
    ...     output_clean="m51_clean.fits",
    ...     output_rejected="m51_rejected.fits",
    ...     remove_column="AOR_ID",
    ...     remove_values=["04_0116_0020609", "04_0116_0020506"]
    ... )
    >>> print(f"Clean: {clean_path}, Rejected: {rejected_path}")
    """
    # Find ALL binary table HDUs with spectra
    matrix_hdus = []
    for idx, hdu in enumerate(hdul):
        if hasattr(hdu, 'data') and hdu.data is not None:
            if 'SPECTRUM' in hdu.data.dtype.names:
                matrix_hdus.append((idx, hdu))
    
    if not matrix_hdus:
        raise ValueError("No HDU with SPECTRUM column found")
    
    # Combine data from all HDUs
    all_data = []
    for idx, hdu in matrix_hdus:
        all_data.append(hdu.data)
    
    # Concatenate all data
    data = np.concatenate(all_data)
    
    # Apply removal filter FIRST if specified (before object separation)
    removed_data = np.array([])  # Track removed rows for rejected file
    if remove_column and remove_values:
        if remove_column not in data.dtype.names:
            raise ValueError(f"Column '{remove_column}' not found in FITS data")
        
        # Create a mask for values to KEEP (inverse of remove)
        keep_mask = np.ones(len(data), dtype=bool)
        for i, row in enumerate(data):
            value = row[remove_column]
            # Handle both bytes and string values
            if isinstance(value, bytes):
                value = value.decode('utf-8').strip()
            else:
                value = str(value).strip()
            
            if value in remove_values:
                keep_mask[i] = False
        
        # Save the removed data for the rejected file
        removed_data = data[~keep_mask]
        data = data[keep_mask]
    
    # Separate data by object
    target_mask = np.array([object_name.lower() in str(obj).lower() 
                           for obj in data['OBJECT']])
    other_mask = ~target_mask
    
    other_data = data[other_mask]
    target_data = data[target_mask]
    
    # Filter target object by NaN content
    nan_fractions = []
    for spectrum in target_data['SPECTRUM']:
        _, frac = detect_nan_channels(spectrum)
        nan_fractions.append(frac)
    
    nan_fractions = np.array(nan_fractions)
    clean_target_mask = nan_fractions < nan_threshold
    rejected_target_mask = ~clean_target_mask
    
    target_clean = target_data[clean_target_mask]
    target_rejected = target_data[rejected_target_mask]
    
    # Add the removed data to rejected
    if len(removed_data) > 0:
        target_rejected = np.concatenate([target_rejected, removed_data])
    
    # Combine other objects with clean target objects
    clean_combined = np.concatenate([other_data, target_clean])
    
    # Set default output paths
    if output_clean is None:
        output_clean = Path("clean_data.fits")
    else:
        output_clean = Path(output_clean)
    
    if output_rejected is None:
        output_rejected = Path("rejected_data.fits")
    else:
        output_rejected = Path(output_rejected)
    
    # Create FITS files with same structure as original
    # Copy primary HDU
    primary_hdu = hdul[0].copy()
    
    # Create binary table HDU for clean data using Table
    from astropy.table import Table
    
    # Convert to table, then back to FITS
    table_clean = Table(clean_combined)
    hdu_clean = fits.BinTableHDU(table_clean)
    hdu_clean.name = matrix_hdus[0][1].name  # Use name from first spectrum HDU
    
    # Copy header information from original
    for key in matrix_hdus[0][1].header:
        if key not in ['NAXIS1', 'NAXIS2', 'TFIELDS'] and key != '':
            try:
                hdu_clean.header[key] = matrix_hdus[0][1].header[key]
            except (ValueError, KeyError):
                pass
    
    # Create FITS file for clean data
    hdul_clean = fits.HDUList([primary_hdu, hdu_clean])
    hdul_clean.writeto(output_clean, overwrite=True)
    
    # Create binary table HDU for rejected data (only if there are rejected spectra)
    if len(target_rejected) > 0:
        table_rejected = Table(target_rejected)
        hdu_rejected = fits.BinTableHDU(table_rejected)
        hdu_rejected.name = matrix_hdus[0][1].name
        
        # Copy header information from original
        for key in matrix_hdus[0][1].header:
            if key not in ['NAXIS1', 'NAXIS2', 'TFIELDS'] and key != '':
                try:
                    hdu_rejected.header[key] = matrix_hdus[0][1].header[key]
                except (ValueError, KeyError):
                    pass
        
        # Create FITS file for rejected data
        hdul_rejected = fits.HDUList([primary_hdu, hdu_rejected])
        hdul_rejected.writeto(output_rejected, overwrite=True)
    else:
        # Create empty rejected file if no rejections
        primary_rejected = fits.PrimaryHDU()
        hdul_rejected = fits.HDUList([primary_rejected])
        hdul_rejected.writeto(output_rejected, overwrite=True)
    
    return output_clean, output_rejected


def baseline_subtract(spectrum: np.ndarray,
                      order: int = 1,
                      mask_nan: bool = True,
                      sigma_clip: bool = True,
                      niter: int = 3,
                      clip_sigma: float = 3.0,
                      window: Optional[Tuple[int, int]] = None) -> np.ndarray:
    """
    Subtract a polynomial baseline from a 1D spectrum.

    The function fits a polynomial of given `order` to the spectrum while
    ignoring NaNs, optionally an emission window, and (optionally) outlier 
    pixels using iterative sigma-clipping. The baseline is evaluated across 
    the full spectral axis and subtracted from the input. The input NaN 
    positions are preserved in the output.

    Parameters
    ----------
    spectrum : np.ndarray
        1D input spectrum array (floats). Can contain NaNs.
    order : int, optional
        Polynomial order to fit (default: 1 = linear).
    mask_nan : bool, optional
        If True, ignore NaN values when fitting (default: True).
    sigma_clip : bool, optional
        If True, perform iterative sigma-clipping to mask emission/line
        features before the final baseline fit (default: True).
    niter : int, optional
        Maximum number of sigma-clipping iterations (default: 3).
    clip_sigma : float, optional
        Sigma threshold for clipping (default: 3.0).
    window : tuple of (int, int), optional
        Channel/pixel range [start, end) to exclude from baseline fitting.
        Useful to ignore strong emission lines or absorption features.
        Default: None (no exclusion).

    Returns
    -------
    np.ndarray
        Baseline-subtracted spectrum (same shape as input). Input NaNs
        are preserved in the returned array.

    Notes
    -----
    - Uses numpy.polyfit for the polynomial fit. For very large orders or
      pathological data, consider using a more robust fitter.
    - The function is intentionally conservative: it only masks outliers
      relative to the current baseline estimate during sigma-clipping.
    - The window is always excluded from fitting, separate from
      sigma-clipping.

    Examples
    --------
    >>> corrected = baseline_subtract(spectrum, order=2)
    >>> corrected = baseline_subtract(spectrum, order=2, window=(90, 110))
    """
    # Ensure float copy
    y = np.array(spectrum, dtype=float)
    x = np.arange(y.size)

    # Initial good-pixel mask: finite values
    if mask_nan:
        good = np.isfinite(y)
    else:
        good = np.ones_like(y, dtype=bool)

    # Exclude pixels in the window
    if window is not None:
        start, end = window
        start = max(0, int(start))
        end = min(len(y), int(end))
        good[start:end] = False

    # If nothing to fit, return input
    if np.sum(good) <= order:
        return y

    # Iterative sigma-clipping to avoid lines influencing the baseline
    if sigma_clip:
        mask = good.copy()
        for _ in range(niter):
            try:
                coeffs = np.polyfit(x[mask], y[mask], order)
            except Exception:
                break
            baseline = np.polyval(coeffs, x)
            resid = y - baseline
            # Compute robust sigma only on currently masked good pixels
            sigma = np.nanstd(resid[mask])
            if sigma == 0 or not np.isfinite(sigma):
                break
            new_mask = mask & (np.abs(resid) <= clip_sigma * sigma)
            # Keep finite values forced
            new_mask &= np.isfinite(y)
            if new_mask.sum() == mask.sum():
                break
            mask = new_mask
        # Final fit on clipped mask
        if mask.sum() <= order:
            # fallback to mask=good
            mask = good
        coeffs = np.polyfit(x[mask], y[mask], order)
    else:
        coeffs = np.polyfit(x[good], y[good], order)

    baseline = np.polyval(coeffs, x)
    corrected = y - baseline

    # Preserve NaNs from input
    corrected[~np.isfinite(y)] = np.nan

    return corrected


def apply_baseline_to_hdul(hdul: fits.HDUList,
                           order: int = 1,
                           spectrum_column: str = 'SPECTRUM',
                           window: Optional[Tuple[int, int]] = None,
                           inplace: bool = False) -> fits.HDUList:
    """
    Apply baseline subtraction to every spectrum in the FITS table HDU.

    Parameters
    ----------
    hdul : astropy.io.fits.HDUList
        FITS HDUList containing a binary table HDU with a spectrum column.
    order : int, optional
        Polynomial order for baseline subtraction (default: 1).
    spectrum_column : str, optional
        Name of the column that contains the spectrum arrays (default: 'SPECTRUM').
    window : tuple of (int, int), optional
        Channel range [start, end) to exclude from baseline fitting across all spectra.
        Default: None (no exclusion).
    inplace : bool, optional
        If True, modify the provided HDUList in-place and return it. If False,
        a new HDUList is returned and the original is left unmodified.

    Returns
    -------
    astropy.io.fits.HDUList
        HDUList with baseline-subtracted spectra in the table HDU.
    """
    # Locate binary table HDU containing the spectrum column
    matrix_hdu = None
    matrix_hdu_index = None
    for idx, hdu in enumerate(hdul):
        if hasattr(hdu, 'data') and hdu.data is not None:
            if spectrum_column in hdu.data.dtype.names:
                matrix_hdu = hdu
                matrix_hdu_index = idx
                break

    if matrix_hdu is None:
        raise ValueError(f"No HDU with '{spectrum_column}' column found")

    data = matrix_hdu.data

    # Build new spectra array by applying baseline_subtract per row
    new_spectra = []
    for spec in data[spectrum_column]:
        corrected = baseline_subtract(np.array(spec, dtype=float), 
                                      order=order, 
                                      window=window)
        new_spectra.append(corrected)

    # Convert to astropy Table for easy column replacement
    from astropy.table import Table

    table = Table(data)
    # astropy Table accepts a 2D numpy array (nrows, nspec) for the column
    new_col = np.vstack(new_spectra)
    table[spectrum_column] = new_col

    # Build new HDUList preserving primary HDU
    primary = hdul[0].copy()
    new_hdu = fits.BinTableHDU(table)
    new_hdu.name = matrix_hdu.name

    # Copy header keywords where it makes sense
    for key in matrix_hdu.header:
        if key not in ['NAXIS1', 'NAXIS2', 'TFIELDS'] and key != '':
            try:
                new_hdu.header[key] = matrix_hdu.header[key]
            except (ValueError, KeyError):
                pass

    new_hdul = fits.HDUList([primary, new_hdu])

    if inplace:
        # Replace the data in the original HDUList
        hdul[matrix_hdu_index] = new_hdu
        hdul[0] = primary
        return hdul

    return new_hdul


def apply_baseline_from_config(config_path: Optional[Union[str, Path]] = None,
                               hdul: Optional[fits.HDUList] = None,
                               output_path: Optional[Union[str, Path]] = None,
                               overwrite: bool = True,
                               window: Optional[Tuple[int, int]] = None) -> Path:
    """
    Read baseline order from config.toml and apply baseline subtraction.

    Behavior:
    - If `hdul` is None the FITS file is read from the config [input].fits_file
    - The baseline order is read from [reduction].baseline or
      [reduction].baseline_order (both supported). If not present, defaults to 1.
    - The window is read from [reduction].window (format: [start, end])
      or can be passed as a parameter (parameter overrides config).
    - The resulting reduced FITS is written to `output_path` if provided,
      otherwise to config [output].reduced_fits or 'reduced_data.fits'.

    Parameters
    ----------
    config_path : str or Path, optional
        Path to config.toml. If None, searches in project hierarchy.
    hdul : astropy.io.fits.HDUList, optional
        FITS HDUList to process. If None, reads from config [input].fits_file.
    output_path : str or Path, optional
        Output FITS file path. If None, uses config [output].reduced_fits.
    overwrite : bool, optional
        If True, overwrite output file if it exists (default: True).
    window : tuple of (int, int), optional
        Channel range [start, end) to exclude from baseline fitting.
        If provided, overrides config [reduction].window.

    Returns
    -------
    Path
        Path to the written reduced FITS file.
    """
    # Local import to avoid circular imports at module import time
    try:
        from ..basic_io import get_config, read_fits_from_config
    except Exception:
        from oi_zeigt.basic_io import get_config, read_fits_from_config

    if hdul is None:
        # Read HDU from config
        cfg = get_config(config_path)
        fits_file = cfg.get('input', {}).get('fits_file')
        if fits_file is None:
            raise KeyError("'fits_file' not found in [input] section of config")
        hdul = read_fits_from_config(config_path)
    else:
        cfg = get_config(config_path) if config_path is not None else {}

    # Read baseline order from config (support both keys)
    reduction_cfg = cfg.get('reduction', {}) if isinstance(cfg, dict) else {}
    order = reduction_cfg.get('baseline_order', reduction_cfg.get('baseline', 1))
    try:
        order = int(order)
    except Exception:
        order = 1

    # Read window from config if not provided as parameter
    if window is None:
        window_cfg = reduction_cfg.get('window', None)
        if window_cfg is not None:
            # Try to parse as list/tuple [start, end]
            try:
                if isinstance(window_cfg, (list, tuple)) and len(window_cfg) == 2:
                    window = (int(window_cfg[0]), int(window_cfg[1]))
            except (ValueError, TypeError, IndexError):
                window = None

    # Apply baseline subtraction
    new_hdul = apply_baseline_to_hdul(hdul, order=order, 
                                      window=window, inplace=False)

    # Determine output path
    if output_path is None:
        output_path = cfg.get('output', {}).get('reduced_fits', 'reduced_data.fits')

    output_path = Path(output_path)
    new_hdul.writeto(output_path, overwrite=overwrite)

    return output_path


def smooth_spectrum(spectrum: np.ndarray, window_size: int = 5) -> np.ndarray:
    """
    Apply boxcar smoothing to a spectrum.
    
    Parameters
    ----------
    spectrum : np.ndarray
        Input spectrum array.
    window_size : int, optional
        Size of the smoothing window (default: 5).
    
    Returns
    -------
    np.ndarray
        Smoothed spectrum.
    
    Examples
    --------
    >>> smoothed = smooth_spectrum(spectrum, window_size=7)
    """
    kernel = np.ones(window_size) / window_size
    smoothed = np.convolve(spectrum, kernel, mode='same')
    
    return smoothed


def reduce_spectra(hdul: fits.HDUList,
                   spectrum_column: str = 'SPECTRUM',
                   methods: Optional[dict] = None) -> fits.HDUList:
    """
    Apply a series of reduction methods to spectra in a FITS HDU list.

    This is a general-purpose function that orchestrates multiple reduction
    steps (baseline subtraction, smoothing, etc.) in sequence.

    Parameters
    ----------
    hdul : astropy.io.fits.HDUList
        FITS HDUList containing a binary table HDU with spectra.
    spectrum_column : str, optional
        Name of the spectrum column (default: 'SPECTRUM').
    methods : dict, optional
        Dictionary specifying which reduction methods to apply and their parameters.
        Format: {'unblank': {}, 'baseline': {...}, 'smooth': {...}, ...}
        Methods are applied in order: unblank → baseline → smooth
        Supported methods:
        - 'unblank': Fill NaN values using linear interpolation.
          Parameters: {} (no parameters needed)
        - 'baseline': Apply polynomial baseline subtraction.
          Parameters: {'order': int, 'window': (int, int) or None}
        - 'smooth': Apply boxcar smoothing.
          Parameters: {'window_size': int}
        If None, returns a copy of hdul unchanged.

    Returns
    -------
    astropy.io.fits.HDUList
        HDUList with reduced spectra, preserving original data types.

    Examples
    --------
    >>> methods = {
    ...     'baseline': {'order': 2, 'window': (100, 120)},
    ...     'smooth': {'window_size': 5}
    ... }
    >>> reduced_hdul = reduce_spectra(hdul, methods=methods)
    """
    if methods is None:
        methods = {}

    # Locate binary table HDU
    matrix_hdu = None
    matrix_hdu_index = None
    for idx, hdu in enumerate(hdul):
        if hasattr(hdu, 'data') and hdu.data is not None:
            if spectrum_column in hdu.data.dtype.names:
                matrix_hdu = hdu
                matrix_hdu_index = idx
                break

    if matrix_hdu is None:
        raise ValueError(f"No HDU with '{spectrum_column}' column found")

    data = matrix_hdu.data

    # Save the original data type of the spectrum column
    original_spectrum_dtype = data[spectrum_column].dtype
    
    # Start with the original spectra (convert to float for processing)
    spectra = np.array([np.array(spec, dtype=float) for spec in data[spectrum_column]])
    
    # Track which rows to keep (for when we filter spectra)
    keep_mask = np.ones(len(spectra), dtype=bool)

    # Apply each reduction method in sequence
    if 'unblank' in methods:
        params = methods['unblank']
        nan_threshold = params.get('nan_threshold', 0.20)
        # Filter spectra and track which ones we keep
        filtered_spectra = _reduce_unblank(spectra, nan_threshold=nan_threshold)
        # Figure out which rows were kept
        nan_fractions = []
        for spec in spectra:
            nan_mask = np.isnan(spec)
            nan_frac = np.sum(nan_mask) / len(spec)
            nan_fractions.append(nan_frac)
        nan_fractions = np.array(nan_fractions)
        keep_mask = nan_fractions < nan_threshold
        spectra = filtered_spectra

    if 'baseline' in methods:
        params = methods['baseline']
        order = params.get('order', 1)
        window = params.get('window', None)
        spectra = _reduce_baseline(spectra, order=order, window=window)

    if 'smooth' in methods:
        params = methods['smooth']
        window_size = params.get('window_size', 5)
        spectra = _reduce_smooth(spectra, window_size=window_size)

    # Convert spectra back to original data type to preserve file size
    spectra = spectra.astype(original_spectrum_dtype)

    # Create new HDU with reduced spectra
    from astropy.table import Table

    # Filter the data table to match the filtered spectra (if unblank was used)
    if 'unblank' in methods:
        filtered_data = data[keep_mask]
    else:
        filtered_data = data
    
    table = Table(filtered_data)
    # Update the spectrum column with reduced spectra (now in original dtype)
    table[spectrum_column] = spectra

    # Build new HDUList
    primary = hdul[0].copy()
    new_hdu = fits.BinTableHDU(table)
    new_hdu.name = matrix_hdu.name

    # Copy header keywords from original
    for key in matrix_hdu.header:
        if key not in ['NAXIS1', 'NAXIS2', 'TFIELDS'] and key != '':
            try:
                new_hdu.header[key] = matrix_hdu.header[key]
            except (ValueError, KeyError):
                pass

    new_hdul = fits.HDUList([primary, new_hdu])

    return new_hdul


def _reduce_unblank(spectra: np.ndarray,
                    nan_threshold: float = 0.20) -> np.ndarray:
    """
    Filter out spectra with too many NaN (blank) channels.

    This function removes spectra that have a NaN fraction >= nan_threshold,
    keeping only spectra with sufficient valid data. This is useful for 
    preparing clean datasets before baseline fitting.

    Parameters
    ----------
    spectra : np.ndarray
        2D array of shape (nspectra, nchannel).
    nan_threshold : float, optional
        Maximum acceptable NaN fraction (default: 0.20 = 20%).
        Spectra with NaN fraction >= threshold are removed.

    Returns
    -------
    np.ndarray
        2D array with spectra containing >= nan_threshold NaNs removed.
    """
    # Calculate NaN fraction for each spectrum
    nan_fractions = []
    for spec in spectra:
        nan_mask = np.isnan(spec)
        nan_frac = np.sum(nan_mask) / len(spec)
        nan_fractions.append(nan_frac)
    
    nan_fractions = np.array(nan_fractions)
    
    # Keep only spectra with NaN fraction < threshold
    keep_mask = nan_fractions < nan_threshold
    
    # Return filtered spectra
    return spectra[keep_mask]


def _reduce_baseline(spectra: np.ndarray,
                     order: int = 1,
                     window: Optional[Tuple[int, int]] = None) -> np.ndarray:
    """
    Apply baseline subtraction to an array of spectra (2D array).

    Parameters
    ----------
    spectra : np.ndarray
        2D array of shape (nspectra, nchannel).
    order : int, optional
        Polynomial order (default: 1).
    window : tuple of (int, int), optional
        Channel range to exclude from fitting.

    Returns
    -------
    np.ndarray
        2D array of baseline-subtracted spectra.
    """
    result = np.zeros_like(spectra)
    for i, spec in enumerate(spectra):
        result[i] = baseline_subtract(spec, order=order, window=window)
    return result


def _reduce_smooth(spectra: np.ndarray,
                   window_size: int = 5) -> np.ndarray:
    """
    Apply smoothing to an array of spectra (2D array).

    Parameters
    ----------
    spectra : np.ndarray
        2D array of shape (nspectra, nchannel).
    window_size : int, optional
        Size of smoothing kernel (default: 5).

    Returns
    -------
    np.ndarray
        2D array of smoothed spectra.
    """
    result = np.zeros_like(spectra)
    for i, spec in enumerate(spectra):
        result[i] = smooth_spectrum(spec, window_size=window_size)
    return result


def reduce_spectra_from_config(config_path: Optional[Union[str, Path]] = None,
                               hdul: Optional[fits.HDUList] = None,
                               output_path: Optional[Union[str, Path]] = None,
                               overwrite: bool = True,
                               methods: Optional[dict] = None) -> Path:
    """
    Read reduction methods from config.toml and perform reduction.

    Behavior:
    - If `hdul` is None, read FITS from [input].fits_file
    - If `output_path` is None, try config [output].reduced_fits
    - If no config output path, use 'reduced_data.fits' and print a warning
    - Reduction methods can be specified via the `methods` parameter or config

    Parameters
    ----------
    config_path : str or Path, optional
        Path to config.toml.
    hdul : astropy.io.fits.HDUList, optional
        FITS HDUList to reduce. If None, read from config.
    output_path : str or Path, optional
        Output file path. If None, read from config or use default.
    overwrite : bool, optional
        If True, overwrite output file (default: True).
    methods : dict, optional
        Reduction methods dict. If provided, overrides config settings.

    Returns
    -------
    Path
        Path to the written reduced FITS file.
    """
    # Local import to avoid circular imports
    try:
        from ..basic_io import get_config, read_fits_from_config
    except Exception:
        from oi_zeigt.basic_io import get_config, read_fits_from_config

    if hdul is None:
        cfg = get_config(config_path)
        hdul = read_fits_from_config(config_path)
    else:
        cfg = get_config(config_path) if config_path is not None else {}

    # Determine output path with warning if using default
    if output_path is None:
        output_path = cfg.get('output', {}).get('reduced_fits', None)
        if output_path is None:
            output_path = 'reduced_data.fits'
            print(f"\n⚠ WARNING: No output path specified in config or command line.")
            print(f"  Using default: {output_path}\n")

    output_path = Path(output_path)

    # Apply reduction
    reduced_hdul = reduce_spectra(hdul, methods=methods)

    # Write output
    reduced_hdul.writeto(output_path, overwrite=overwrite)

    return output_path


def average_spectra(hdul: fits.HDUList,
                   spectrum_column: str = 'SPECTRUM',
                   group_by: Optional[str] = None) -> dict:
    """
    Compute average spectrum(s) from a FITS HDUList.

    Can compute a single average of all spectra, or group-averaged spectra
    if a grouping column is specified.

    Parameters
    ----------
    hdul : astropy.io.fits.HDUList
        FITS HDUList containing spectral data.
    spectrum_column : str, optional
        Name of the spectrum column (default: 'SPECTRUM').
    group_by : str, optional
        Column name to group by (e.g., 'OBJECT'). If None, compute
        single global average.

    Returns
    -------
    dict
        If group_by is None:
        {
            'avg_spectrum': np.ndarray,
            'std_spectrum': np.ndarray,
            'count': int,
            'rms': float
        }
        
        If group_by is specified:
        {
            'group_name_1': {...},
            'group_name_2': {...},
            ...
        }
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

    data = matrix_hdu.data
    spectra = np.array([np.array(spec, dtype=float) for spec in data[spectrum_column]])

    if group_by is None:
        # Single global average
        avg_spectrum = np.nanmean(spectra, axis=0)
        std_spectrum = np.nanstd(spectra, axis=0)
        count = len(spectra)
        rms = np.nanstd(spectra, axis=0).mean()  # Mean RMS across channels
        
        return {
            'avg_spectrum': avg_spectrum,
            'std_spectrum': std_spectrum,
            'count': count,
            'rms': rms
        }
    else:
        # Group-based averaging
        if group_by not in data.dtype.names:
            raise ValueError(f"Column '{group_by}' not found in data")
        
        groups = data[group_by]
        unique_groups = np.unique(groups)
        result = {}
        
        for group_val in unique_groups:
            # Decode if bytes
            group_name = group_val.decode().strip() if isinstance(group_val, bytes) else str(group_val).strip()
            mask = groups == group_val
            group_spectra = spectra[mask]
            
            avg_spectrum = np.nanmean(group_spectra, axis=0)
            std_spectrum = np.nanstd(group_spectra, axis=0)
            count = len(group_spectra)
            rms = np.nanstd(group_spectra, axis=0).mean()
            
            result[group_name] = {
                'avg_spectrum': avg_spectrum,
                'std_spectrum': std_spectrum,
                'count': count,
                'rms': rms
            }
        
        return result


def average_spectra_from_config(config_path: Optional[Union[str, Path]] = None,
                                hdul: Optional[fits.HDUList] = None,
                                output_fits: Optional[Union[str, Path]] = None,
                                group_by: Optional[str] = None,
                                overwrite: bool = True) -> Path:
    """
    Compute and save average spectrum(s) from config and/or command line options.

    Parameters
    ----------
    config_path : str or Path, optional
        Path to config.toml.
    hdul : astropy.io.fits.HDUList, optional
        FITS HDUList to average. If None, read from config.
    output_fits : str or Path, optional
        Output FITS file path. If None, read from config or use default.
    group_by : str, optional
        Column to group by (e.g., 'OBJECT') before averaging.
    overwrite : bool, optional
        If True, overwrite output file (default: True).

    Returns
    -------
    Path
        Path to the written averaged FITS file.
    """
    try:
        from ..basic_io import get_config, read_fits_from_config
    except Exception:
        from oi_zeigt.basic_io import get_config, read_fits_from_config

    if hdul is None:
        cfg = get_config(config_path)
        hdul = read_fits_from_config(config_path)
    else:
        cfg = get_config(config_path) if config_path is not None else {}

    # Determine output path
    if output_fits is None:
        output_fits = cfg.get('output', {}).get('averaged_fits', None)
        if output_fits is None:
            output_fits = 'averaged_data.fits'
            print(f"\n⚠ WARNING: No output path specified in config or command line.")
            print(f"  Using default: {output_fits}\n")

    output_fits = Path(output_fits)

    # Compute averages
    avg_results = average_spectra(hdul, group_by=group_by)

    # Create output FITS file
    primary = hdul[0].copy()
    
    if isinstance(avg_results, dict) and 'avg_spectrum' in avg_results:
        # Single average: create a simple 1-row table
        avg_data = avg_results
        from astropy.table import Table
        
        table = Table()
        table['SPECTRUM'] = [avg_data['avg_spectrum']]
        table['STD'] = [avg_data['std_spectrum']]
        table['COUNT'] = [avg_data['count']]
        table['RMS'] = [avg_data['rms']]
        
        new_hdu = fits.BinTableHDU(table)
        new_hdu.name = 'AVERAGE'
    else:
        # Multiple group averages: create multi-row table
        from astropy.table import Table
        
        group_names = []
        spectra_list = []
        stds_list = []
        counts = []
        rmss = []
        
        for group_name, data in sorted(avg_results.items()):
            group_names.append(group_name)
            spectra_list.append(data['avg_spectrum'])
            stds_list.append(data['std_spectrum'])
            counts.append(data['count'])
            rmss.append(data['rms'])
        
        table = Table()
        table['OBJECT'] = group_names
        table['SPECTRUM'] = spectra_list
        table['STD'] = stds_list
        table['COUNT'] = counts
        table['RMS'] = rmss
        
        new_hdu = fits.BinTableHDU(table)
        new_hdu.name = 'AVERAGE'

    new_hdul = fits.HDUList([primary, new_hdu])
    new_hdul.writeto(output_fits, overwrite=overwrite)

    return output_fits
