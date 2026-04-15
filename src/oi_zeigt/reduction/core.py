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


def detect_missing_channels(spectrum: np.ndarray, 
                           include_blanks: bool = True,
                           blank_value: float = 0.0,
                           blank_tolerance: float = 1e-10,
                           gildas_blank: bool = False) -> Tuple[np.ndarray, float]:
    """
    Detect missing/blank channels in a spectrum using multiple detection methods.
    
    This function detects channels that are marked as missing/blank using various
    astronomical software conventions:
    
    - **NaN (IEEE standard)**: Used by NumPy, FITS, and most modern software
    - **Zero values**: Sometimes used as blanks in older data or specific archives
    - **GILDAS/CLASS blanks**: Values near -9.99e30 (GILDAS/CLASS convention)
    - **Infinity**: Sometimes used for invalid/missing data
    
    Parameters
    ----------
    spectrum : np.ndarray
        Input spectrum array.
    include_blanks : bool, optional
        If True, also detect blank/zero values (default: True).
        If False, only detect NaNs.
    blank_value : float, optional
        Value to consider as blank (default: 0.0).
        Only used if include_blanks=True.
    blank_tolerance : float, optional
        Tolerance for blank value detection (default: 1e-10).
        Only used if include_blanks=True.
    gildas_blank : bool, optional
        If True, also detect GILDAS-style blanks near -9.99e30 (default: False).
    
    Returns
    -------
    tuple
        - missing_mask : np.ndarray (bool)
            Boolean array indicating missing/blank channels (True = missing).
        - fraction_missing : float
            Fraction of missing/blank channels (0.0 to 1.0).
    
    Notes
    -----
    This function combines multiple blank detection methods:
    
    1. **NaN detection** (always applied):
       - Uses np.isnan() - catches IEEE NaN values
    
    2. **Blank/zero detection** (if include_blanks=True):
       - Detects values near blank_value within blank_tolerance
       - Useful for data that marks blanks as 0.0 or similar
    
    3. **GILDAS blank detection** (if gildas_blank=True):
       - Detects values near -9.99e30 (GILDAS/CLASS convention)
       - Catches extreme negative values used as missing markers
    
    4. **Infinity detection** (always applied):
       - Uses np.isinf() - catches infinite values
    
    Examples
    --------
    >>> # Basic usage - only detect NaNs
    >>> mask, frac = detect_missing_channels(spectrum)
    >>> print(f"Missing fraction: {frac:.2%}")
    
    >>> # Include zero-value blanks
    >>> mask, frac = detect_missing_channels(spectrum, include_blanks=True)
    
    >>> # Include GILDAS blanks
    >>> mask, frac = detect_missing_channels(spectrum, gildas_blank=True)
    
    >>> # Full detection - all methods
    >>> mask, frac = detect_missing_channels(
    ...     spectrum, 
    ...     include_blanks=True, 
    ...     gildas_blank=True
    ... )
    """
    # Start with NaN detection (always applied)
    missing_mask = np.isnan(spectrum)
    
    # Add infinity detection
    missing_mask |= np.isinf(spectrum)
    
    # Add blank/zero value detection if requested
    if include_blanks:
        blank_mask = np.abs(spectrum - blank_value) < blank_tolerance
        missing_mask |= blank_mask
    
    # Add GILDAS blank detection if requested
    if gildas_blank:
        # GILDAS uses -9.99e30 as blank marker
        # Use relative tolerance for very large numbers
        gildas_blank_mask = np.abs(spectrum) > 1e30
        missing_mask |= gildas_blank_mask
    
    fraction_missing = np.sum(missing_mask) / len(spectrum)
    
    return missing_mask, fraction_missing


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


def split_fits_by_mission(fits_path: Union[str, Path],
                          output_dir: Union[str, Path],
                          manifest_path: Optional[Union[str, Path]] = None,
                          manifest_science_path: Optional[Union[str, Path]] = None,
                          split_by_scan: bool = True) -> Tuple[list, list]:
    """
    Split a large FITS file into per-(MISSION_ID, TELESCOP, SCAN) chunks.

    TELESCOP is the receiver/pixel identifier (e.g. ``LFAV_PX00_S``).
    Each scan contains its own complete set of calibration rows (TSYS,
    TAU_SIG, SKYCHOPDIFF, …) alongside the science rows, so scan-level
    splitting is safe.

    Chunks that contain only ``TREC (SSB)`` rows are written but kept in a
    separate list so the pipeline can skip them for science-processing steps.

    TSYS_INDEX and TAU_SIG_INDEX (if present) are remapped to be valid
    within each chunk.

    Parameters
    ----------
    fits_path : str or Path
        Input FITS file (opened with memmap — no full-file RAM copy).
    output_dir : str or Path
        Directory where chunk files are written.
    manifest_path : str or Path, optional
        Text file listing *all* chunk paths, one per line.
    manifest_science_path : str or Path, optional
        Text file listing only chunks that contain science rows
        (i.e. not TREC-only chunks).  Pass this to the pipeline steps.
    split_by_scan : bool, optional
        If True (default) split by (MISSION_ID, TELESCOP, SCAN).
        If False split by (MISSION_ID, TELESCOP) only.

    Returns
    -------
    all_paths : list of Path
        All written chunk files.
    science_paths : list of Path
        Subset that contain at least one non-``TREC (SSB)`` row.
    """
    fits_path  = Path(fits_path)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    with fits.open(fits_path, memmap=True) as hdul:
        spec_hdu = None
        for hdu in hdul:
            if (hasattr(hdu, 'data') and hdu.data is not None
                    and hasattr(hdu.data, 'dtype')
                    and 'SPECTRUM' in hdu.data.dtype.names):
                spec_hdu = hdu
                break
        if spec_hdu is None:
            raise ValueError("No HDU with SPECTRUM column found in input file")

        data        = spec_hdu.data
        orig_header = spec_hdu.header
        hdu_name    = spec_hdu.name
        primary_hdu = hdul[0].copy()
        n_rows      = len(data)

        def _to_str_array(col):
            raw = data[col]
            if raw.dtype.kind in ('S', 'U'):
                return np.char.strip(raw.astype('U'))
            return np.array([str(v).strip() for v in raw])

        if 'MISSION_ID' not in data.dtype.names:
            raise ValueError("FITS file has no MISSION_ID column")

        mission_ids = _to_str_array('MISSION_ID')
        telescopes  = _to_str_array('TELESCOP') if 'TELESCOP' in data.dtype.names \
                      else np.full(n_rows, 'UNKNOWN', dtype='U20')
        objects     = _to_str_array('OBJECT')

        if split_by_scan and 'SCAN' in data.dtype.names:
            scans  = data['SCAN'].astype(np.int64)
            groups = sorted(set(zip(mission_ids, telescopes, scans)))
        else:
            scans  = None
            groups = sorted(set(zip(mission_ids, telescopes)))

        has_tsys = 'TSYS_INDEX'    in data.dtype.names
        has_tau  = 'TAU_SIG_INDEX' in data.dtype.names

        safe = lambda s: str(s).replace('/', '_').replace(' ', '_').replace(':', '-')

        all_paths     = []
        science_paths = []

        for group in groups:
            if len(group) == 3:
                mission_id, telescop, scan = group
                mask = (mission_ids == mission_id) & (telescopes == telescop) & (scans == scan)
                chunk_name = f"{safe(mission_id)}__{safe(telescop)}__scan{scan:05d}.fits"
            else:
                mission_id, telescop = group
                mask = (mission_ids == mission_id) & (telescopes == telescop)
                chunk_name = f"{safe(mission_id)}__{safe(telescop)}.fits"

            orig_indices = np.where(mask)[0]
            chunk = data[mask].copy()

            # Remap row-index columns if present
            if has_tsys or has_tau:
                lut_keys = orig_indices.astype(np.int64)
                lut_vals = np.arange(len(orig_indices), dtype=np.int32)

                def _remap(col_name):
                    old      = chunk[col_name].astype(np.int64)
                    new      = np.full(len(old), -1, dtype=np.int32)
                    pos_mask = old >= 0
                    if np.any(pos_mask):
                        ins = np.searchsorted(lut_keys, old[pos_mask])
                        hit = (ins < len(lut_keys)) & (lut_keys[ins] == old[pos_mask])
                        tmp = np.full(int(np.sum(pos_mask)), -1, dtype=np.int32)
                        tmp[hit] = lut_vals[ins[hit]]
                        new[pos_mask] = tmp
                    chunk[col_name] = new

                if has_tsys:
                    _remap('TSYS_INDEX')
                if has_tau:
                    _remap('TAU_SIG_INDEX')

            out_path = output_dir / chunk_name
            hdu_out  = fits.BinTableHDU(chunk)
            hdu_out.name = hdu_name
            for key in orig_header:
                if key not in ('NAXIS1', 'NAXIS2', 'TFIELDS', ''):
                    try:
                        hdu_out.header[key] = orig_header[key]
                    except (ValueError, KeyError):
                        pass

            fits.HDUList([primary_hdu, hdu_out]).writeto(out_path, overwrite=True)
            all_paths.append(out_path)

            # Science chunk = has at least one row that is not TREC (SSB)
            chunk_objects = objects[mask]
            if not np.all(chunk_objects == 'TREC (SSB)'):
                science_paths.append(out_path)

            del chunk

    def _write_manifest(path, paths):
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        with open(path, 'w') as f:
            for p in paths:
                f.write(str(p.resolve()) + '\n')

    if manifest_path is not None:
        _write_manifest(manifest_path, all_paths)
    if manifest_science_path is not None:
        _write_manifest(manifest_science_path, science_paths)

    return all_paths, science_paths


def filter_and_save_fits(hdul: fits.HDUList,
                        object_name: str,
                        nan_threshold: float = 0.20,
                        output_clean: Optional[Union[str, Path]] = None,
                        output_rejected: Optional[Union[str, Path]] = None,
                        remove_column: Optional[str] = None,
                        remove_values: Optional[list] = None,
                        apply_to_all: bool = False,
                        filter_zero_spectra: bool = False,
                        param_filters: Optional[list] = None,
                        spectrum_peak_threshold: Optional[float] = None,
                        filter_tau: bool = False,
                        filter_flights: Optional[list] = None,
                        filter_object_exact: Optional[list] = None,
                        filter_out_object_exact: Optional[list] = None,
                        exclude_obsmode: Optional[list] = None) -> Tuple[Path, Path]:
    """
    Filter FITS data by object and NaN content, with optional column value removal and zero-spectrum filtering.
    
    This function:
    1. Processes ALL HDUs with SPECTRUM columns (important for combined files)
    2. Optionally separates data into target object and other objects (unless apply_to_all=True)
    3. Filters spectra to keep only those with < nan_threshold NaN channels
       - If apply_to_all=False: Only filters target object spectra (default behavior)
       - If apply_to_all=True: Filters ALL spectra regardless of object
    4. Optionally removes spectra where all channels are 0 (if filter_zero_spectra=True)
    5. Optionally removes rows where a specified column matches given values
    6. Saves two FITS files:
       - output_clean: Filtered spectra (clean data)
       - output_rejected: Spectra that were filtered out or removed
    
    Parameters
    ----------
    hdul : astropy.io.fits.HDUList
        Input FITS HDU list containing the data.
    object_name : str
        Name of the object to filter (e.g., "M51").
        Only used when apply_to_all=False.
    nan_threshold : float, optional
        Maximum acceptable NaN fraction (default: 0.20 = 20%).
        Keep spectra with NaN fraction < threshold.
        Reject spectra with NaN fraction >= threshold.
    output_clean : str or Path, optional
        Output path for clean FITS file.
        If None, uses "clean_data.fits".
    output_rejected : str or Path, optional
        Output path for rejected FITS file.
        If None, uses "rejected_data.fits".
    remove_column : str, optional
        Column name to check for removal (e.g., "AOR_ID").
        If specified, rows matching remove_values will be removed.
    remove_values : list, optional
        List of values to remove from remove_column.
        Rows with these values will be moved to rejected file.
    apply_to_all : bool, optional
        If True, apply NaN filtering to ALL spectra (regardless of object).
        If False (default), apply NaN filtering only to target object spectra.
    filter_zero_spectra : bool, optional
        If True, also filter out spectra where all channels are 0.
        These spectra are moved to the rejected file.
        Default: False (keep zero spectra).
    param_filters : list of tuples, optional
        List of ``(column, direction, value)`` filter conditions, where
        ``direction`` is ``'below'`` (keep rows where column >= value) or
        ``'above'`` (keep rows where column <= value).
        Rows with non-finite values in the filter column always pass.
        All conditions are ANDed together.
        Example: ``[('RMSRATIOB', 'below', 2.0), ('RMSRATIOB', 'above', 0.5)]``
    
    Returns
    -------
    tuple
        - clean_path : Path
            Path to the clean FITS file
        - rejected_path : Path
            Path to the rejected FITS file
        - stats : dict
            Dictionary containing filtering statistics:
            - 'clean_count': Number of spectra in clean file
            - 'rejected_total': Total number of rejected spectra
            - 'rejected_nan': Count of spectra rejected for NaN threshold violation
            - 'rejected_zero': Count of spectra rejected for being all-zero
            - 'rejected_removed': Count of spectra rejected by --remove criteria
    
    Examples
    --------
    >>> hdul = fits.open('original.fits')
    >>> # Filter only M51 spectra (default)
    >>> clean_path, rejected_path = filter_and_save_fits(
    ...     hdul, 
    ...     object_name="M51",
    ...     nan_threshold=0.20,
    ...     output_clean="m51_clean.fits",
    ...     output_rejected="m51_rejected.fits",
    ...     apply_to_all=False  # Default: only M51 filtered by NaN
    ... )
    
    >>> # Filter ALL spectra by NaN content
    >>> clean_path, rejected_path = filter_and_save_fits(
    ...     hdul, 
    ...     object_name="M51",  # Still needed for other logic
    ...     nan_threshold=0.20,
    ...     output_clean="all_clean.fits",
    ...     output_rejected="all_rejected.fits",
    ...     apply_to_all=True  # Filter ALL spectra
    ... )
    
    >>> # Also filter out all-zero spectra
    >>> clean_path, rejected_path = filter_and_save_fits(
    ...     hdul, 
    ...     object_name="M51",
    ...     nan_threshold=0.20,
    ...     output_clean="clean_no_zeros.fits",
    ...     output_rejected="rejected_with_zeros.fits",
    ...     filter_zero_spectra=True  # Remove all-zero spectra
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
    
    # Combine data from all HDUs.
    # When there is only one spectrum HDU keep the memory-mapped reference so
    # we never load the entire file into RAM up front.  Multiple HDUs still
    # require a concatenation (and therefore a RAM copy).
    if len(matrix_hdus) == 1:
        data = matrix_hdus[0][1].data
    else:
        data = np.concatenate([hdu.data for _, hdu in matrix_hdus])

    # Track the original row index of every row through all filtering steps so
    # that TSYS_INDEX / TAU_SIG_INDEX can be remapped to the new positions in
    # the output file.
    orig_idx = np.arange(len(data), dtype=np.int64)

    param_rejected_data = np.array([])
    param_rejected_count = 0
    param_filter_details = []  # list of (col, direction, value, n_removed) per condition

    def _apply_param_filters(rows):
        """Apply param_filters to a subset of rows.

        Returns (kept_rows, rejected_rows, details, keep_mask).
        keep_mask is a boolean array aligned to the input rows.
        """
        if not param_filters or len(rows) == 0:
            return rows, np.array([]), [], np.ones(len(rows), dtype=bool)
        keep_mask = np.ones(len(rows), dtype=bool)
        details = []
        for col, direction, value in param_filters:
            if col not in rows.dtype.names:
                raise ValueError(f"Filter column '{col}' not found in FITS data")
            col_vals = rows[col].astype(np.float64)
            finite = np.isfinite(col_vals)
            if direction == 'below':
                cond_mask = (~finite) | (col_vals <= value)
            else:
                cond_mask = (~finite) | (col_vals >= value)
            n_removed_by_cond = int(np.sum(keep_mask & ~cond_mask))
            details.append((col, direction, value, n_removed_by_cond))
            keep_mask &= cond_mask
        return rows[keep_mask], rows[~keep_mask], details, keep_mask

    def _remap_index_column(arr, col_name, row_mapping):
        """Remap an integer index column using old→new row mapping. -1 if not found."""
        if col_name not in arr.dtype.names:
            return
        old_vals = arr[col_name].astype(np.int64)
        new_vals = np.full(len(old_vals), -1, dtype=np.int32)
        if row_mapping:
            keys = np.array(sorted(row_mapping.keys()), dtype=np.int64)
            vals = np.array([row_mapping[k] for k in keys], dtype=np.int32)
            pos_mask = old_vals >= 0
            if np.any(pos_mask):
                ins = np.searchsorted(keys, old_vals[pos_mask])
                valid = (ins < len(keys)) & (keys[ins] == old_vals[pos_mask])
                tmp = np.full(int(np.sum(pos_mask)), -1, dtype=np.int32)
                tmp[valid] = vals[ins[valid]]
                new_vals[pos_mask] = tmp
        arr[col_name] = new_vals

    # Chunk size for per-spectrum operations.  Limits peak RAM to
    # ~chunk_size × n_channels × 4 bytes regardless of total file size.
    _CHUNK = 2000

    def _nan_fractions_chunked(spectra_col):
        """NaN fraction per spectrum — processes in chunks to avoid a full float64 copy."""
        n = len(spectra_col)
        result = np.empty(n, dtype=np.float32)
        for i in range(0, n, _CHUNK):
            chunk = np.asarray(spectra_col[i:i + _CHUNK], dtype=np.float32)
            result[i:i + _CHUNK] = np.isnan(chunk).mean(axis=1)
        return result

    def _peak_mask_chunked(spectra_col, threshold):
        """Boolean mask: True where any |channel| > threshold."""
        n = len(spectra_col)
        result = np.zeros(n, dtype=bool)
        for i in range(0, n, _CHUNK):
            chunk = np.asarray(spectra_col[i:i + _CHUNK], dtype=np.float32)
            with np.errstate(invalid='ignore'):
                result[i:i + _CHUNK] = np.any(np.abs(chunk) > threshold, axis=1)
        return result

    def _zero_mask_chunked(spectra_col):
        """Boolean mask: True where spectrum is all-NaN or all-zero (ignoring NaN)."""
        n = len(spectra_col)
        result = np.zeros(n, dtype=bool)
        for i in range(0, n, _CHUNK):
            chunk = np.asarray(spectra_col[i:i + _CHUNK], dtype=np.float32)
            nan_2d = np.isnan(chunk)
            all_nan = np.all(nan_2d, axis=1)
            all_zero = np.all((chunk == 0) | nan_2d, axis=1)
            result[i:i + _CHUNK] = all_nan | all_zero
        return result

    # Apply removal filter if specified (before object separation)
    removed_data = np.array([])  # Track removed rows for rejected file
    if remove_column and remove_values:
        if remove_column not in data.dtype.names:
            raise ValueError(f"Column '{remove_column}' not found in FITS data")

        raw = data[remove_column]
        if raw.dtype.kind in ('S', 'U'):
            str_vals = np.char.strip(raw.astype('U'))
        else:
            str_vals = np.array([str(v).strip() for v in raw])
        keep_mask = ~np.isin(str_vals, remove_values)

        removed_data = data[~keep_mask]
        data     = data[keep_mask]
        orig_idx = orig_idx[keep_mask]

    # Apply flight filter: remove all rows whose MISSION_ID contains any of the given strings
    flight_removed_data = np.array([])
    flight_removed_missions = []
    flight_removed_count = 0
    if filter_flights:
        if 'MISSION_ID' not in data.dtype.names:
            raise ValueError("MISSION_ID column not found in FITS data — cannot use --filter-flight")
        mid_col = np.array([
            s.decode().strip() if isinstance(s, bytes) else str(s).strip()
            for s in data['MISSION_ID']
        ])
        keep_mask = np.ones(len(data), dtype=bool)
        for flight in filter_flights:
            flight_mask = np.array([flight in mid for mid in mid_col])
            matched = sorted(set(mid_col[flight_mask].tolist()))
            flight_removed_missions.extend(matched)
            keep_mask &= ~flight_mask
        flight_removed_data = data[~keep_mask]
        flight_removed_count = int(np.sum(~keep_mask))
        data     = data[keep_mask]
        orig_idx = orig_idx[keep_mask]

    # Apply OBSMODE exclusion filter: remove rows whose OBSMODE exactly matches any of the given strings.
    # Calibration rows (TSYS, TAU_SIG) are always kept regardless.
    obsmode_removed_data = np.array([])
    obsmode_removed_count = 0
    if exclude_obsmode:
        if 'OBSMODE' not in data.dtype.names:
            import warnings
            warnings.warn("--exclude-obsmode: OBSMODE column not found in FITS data — skipping.")
        else:
            obs_col = np.array([
                s.decode().strip() if isinstance(s, bytes) else str(s).strip()
                for s in data['OBSMODE']
            ])
            obj_col_obs = np.array([
                s.decode().strip() if isinstance(s, bytes) else str(s).strip()
                for s in data['OBJECT']
            ]) if 'OBJECT' in data.dtype.names else np.full(len(data), '', dtype='U20')
            cal_mask = np.isin(obj_col_obs, ['TSYS', 'TAU_SIG'])
            remove_mask = np.isin(obs_col, list(exclude_obsmode)) & ~cal_mask
            obsmode_removed_data = data[remove_mask]
            obsmode_removed_count = int(np.sum(remove_mask))
            data     = data[~remove_mask]
            orig_idx = orig_idx[~remove_mask]
            if obsmode_removed_count:
                removed_modes = sorted(set(obs_col[remove_mask].tolist()))
                print(f"OBSMODE exclusion: removed {obsmode_removed_count} rows with OBSMODE in {removed_modes}")

    # Calibration OBJECT values that must never be removed by object filters
    _PROTECTED = {'SKYCHOPDIFF', 'SKYDIFF', 'SKY-DIFF', 'TAU_SIG', 'TSYS'}

    # Apply exact object keep filter: keep only rows whose OBJECT exactly matches one of the given strings.
    # Calibration rows are always kept regardless.
    if filter_object_exact:
        if 'OBJECT' not in data.dtype.names:
            raise ValueError("OBJECT column not found in FITS data — cannot use --filter-object-exact")
        obj_exact = np.array([
            s.decode().strip() if isinstance(s, bytes) else str(s).strip()
            for s in data['OBJECT']
        ])
        keep_mask = np.zeros(len(data), dtype=bool)
        for name in filter_object_exact:
            keep_mask |= (obj_exact == name)
        # Always keep calibration rows
        keep_mask |= np.isin(obj_exact, list(_PROTECTED))
        data     = data[keep_mask]
        orig_idx = orig_idx[keep_mask]

    # Apply exact object remove filter: remove rows whose OBJECT exactly matches one of the given strings.
    # Calibration rows are always kept regardless.
    if filter_out_object_exact:
        if 'OBJECT' not in data.dtype.names:
            raise ValueError("OBJECT column not found in FITS data — cannot use --filter-out-object-exact")
        obj_exact = np.array([
            s.decode().strip() if isinstance(s, bytes) else str(s).strip()
            for s in data['OBJECT']
        ])
        keep_mask = np.ones(len(data), dtype=bool)
        for name in filter_out_object_exact:
            # Only remove if not a protected calibration type
            if name in _PROTECTED:
                import warnings
                warnings.warn(f"--filter-out-object-exact: '{name}' is a protected calibration type and will not be removed.")
            else:
                keep_mask &= (obj_exact != name)
        data     = data[keep_mask]
        orig_idx = orig_idx[keep_mask]

    # Always separate target object from calibration rows first.
    # param_filters are ONLY applied to science (target) rows — never to
    # TSYS, TAU_SIG, SKYCHOPDIFF or any other calibration source.
    obj_col = data['OBJECT']
    if obj_col.dtype.kind in ('S', 'U'):
        obj_strs = np.char.lower(obj_col.astype('U'))
    else:
        obj_strs = np.array([str(o).lower() for o in obj_col])
    target_mask = np.char.find(obj_strs, object_name.lower()) >= 0
    other_data      = data[~target_mask]         # calibration rows — always kept
    other_orig_idx  = orig_idx[~target_mask]
    target_data     = data[target_mask]
    target_orig_idx = orig_idx[target_mask]
    del data, orig_idx  # free the base arrays — all data now lives in other_data / target_data

    # Apply param_filters exclusively to the target (science) rows.
    target_data, pf_rejected, param_filter_details, pf_keep_mask = _apply_param_filters(target_data)
    target_orig_idx = target_orig_idx[pf_keep_mask]
    if len(pf_rejected) > 0:
        param_rejected_data = pf_rejected
        param_rejected_count = len(pf_rejected)

    # Apply spectrum peak filter exclusively to the target (science) rows.
    peaks_rejected_count = 0
    peaks_rejected_data  = np.array([])
    if spectrum_peak_threshold is not None and len(target_data) > 0:
        peak_mask = _peak_mask_chunked(target_data['SPECTRUM'], spectrum_peak_threshold)
        peaks_rejected_data  = target_data[peak_mask]
        peaks_rejected_count = int(np.sum(peak_mask))
        target_data     = target_data[~peak_mask]
        target_orig_idx = target_orig_idx[~peak_mask]

    # Apply --filter-tau: reject science spectra linked to bad TAU_SIG rows
    # (any channel outside [0.001, 1.0]), and also move the bad TAU_SIG rows
    # themselves to the rejected dataset.
    tau_rejected_science_count = 0
    tau_rejected_tau_count = 0
    tau_rejected_data = np.array([])
    tau_filter_details = []  # list of (tau_original_index, n_science_removed)

    if filter_tau and len(other_data) > 0 and 'TAU_SIG_INDEX' in target_data.dtype.names:
        # Identify TAU_SIG rows within other_data
        def _to_str(v):
            return v.decode().strip() if isinstance(v, bytes) else str(v).strip()

        tau_obj_col = other_data['OBJECT']
        if tau_obj_col.dtype.kind in ('S', 'U'):
            tau_obj_strs = np.char.strip(tau_obj_col.astype('U'))
        else:
            tau_obj_strs = np.array([_to_str(v) for v in tau_obj_col])
        tau_sig_pos = tau_obj_strs == 'TAU_SIG'

        if np.any(tau_sig_pos):
            tau_spectra_arr = np.asarray(other_data['SPECTRUM'][tau_sig_pos], dtype=np.float32)
            tau_orig_subset = other_orig_idx[tau_sig_pos]

            # A TAU_SIG spectrum is "bad" if any non-NaN channel is > 1 or < 0.001
            with np.errstate(invalid='ignore'):
                nan_tau = np.isnan(tau_spectra_arr)
                has_valid = ~np.all(nan_tau, axis=1)
                # fill NaN with an in-range value so they don't trigger the condition
                filled = np.where(nan_tau, np.float32(0.5), tau_spectra_arr)
                bad_tau_mask = has_valid & (
                    np.any(filled > np.float32(1.0), axis=1) |
                    np.any(filled < np.float32(0.001), axis=1)
                )

            if np.any(bad_tau_mask):
                bad_tau_orig_set = set(int(v) for v in tau_orig_subset[bad_tau_mask])

                # Remove linked science spectra
                sci_tau_idx = np.array(target_data['TAU_SIG_INDEX'], dtype=np.int64)
                sci_bad_mask = np.isin(sci_tau_idx, list(bad_tau_orig_set))

                # Build per-tau reporting details
                for orig_idx_val in sorted(bad_tau_orig_set):
                    n_sci = int(np.sum(sci_tau_idx == orig_idx_val))
                    tau_filter_details.append((int(orig_idx_val), n_sci))

                tau_rejected_science_count = int(np.sum(sci_bad_mask))
                if tau_rejected_science_count > 0:
                    tau_rejected_data = target_data[sci_bad_mask]
                target_data     = target_data[~sci_bad_mask]
                target_orig_idx = target_orig_idx[~sci_bad_mask]

                # Remove bad TAU_SIG rows from other_data
                full_bad_other = np.zeros(len(other_data), dtype=bool)
                tau_positions_in_other = np.where(tau_sig_pos)[0]
                full_bad_other[tau_positions_in_other[bad_tau_mask]] = True
                tau_rejected_tau_count = int(np.sum(full_bad_other))
                bad_tau_rows = other_data[full_bad_other]
                other_data     = other_data[~full_bad_other]
                other_orig_idx = other_orig_idx[~full_bad_other]

                # Collect all tau-rejected rows (science + bad TAU_SIG spectra)
                if tau_rejected_science_count > 0 and tau_rejected_tau_count > 0:
                    tau_rejected_data = np.concatenate([tau_rejected_data, bad_tau_rows])
                elif tau_rejected_tau_count > 0:
                    tau_rejected_data = bad_tau_rows

    if apply_to_all:
        # NaN filtering applied to ALL spectra (target + calibration)
        combined          = np.concatenate([other_data, target_data]) if len(other_data) > 0 else target_data
        combined_orig_idx = np.concatenate([other_orig_idx, target_orig_idx]) if len(other_data) > 0 else target_orig_idx
        nan_fractions  = _nan_fractions_chunked(combined['SPECTRUM'])
        clean_mask     = nan_fractions < nan_threshold
        clean_combined      = combined[clean_mask]
        clean_combined_orig = combined_orig_idx[clean_mask]
        all_rejected        = combined[~clean_mask]
        nan_rejected_count  = int(np.sum(~clean_mask))
        removed_count = len(removed_data)
        if len(removed_data) > 0:
            all_rejected = np.concatenate([all_rejected, removed_data])
        if flight_removed_count > 0:
            all_rejected = np.concatenate([all_rejected, flight_removed_data]) if len(all_rejected) > 0 else flight_removed_data
        if obsmode_removed_count > 0:
            all_rejected = np.concatenate([all_rejected, obsmode_removed_data]) if len(all_rejected) > 0 else obsmode_removed_data
        del removed_data, combined, combined_orig_idx, other_data, other_orig_idx, target_data, target_orig_idx
    else:
        # NaN filtering applied only to target object (default behavior)
        nan_fractions = _nan_fractions_chunked(target_data['SPECTRUM']) if len(target_data) > 0 else np.array([], dtype=np.float32)
        clean_target_mask = nan_fractions < nan_threshold if len(nan_fractions) > 0 else np.array([], dtype=bool)
        target_clean         = target_data[clean_target_mask]
        target_clean_orig    = target_orig_idx[clean_target_mask]
        target_rejected      = target_data[~clean_target_mask]
        nan_rejected_count   = int(np.sum(~clean_target_mask))
        removed_count = len(removed_data)
        if len(removed_data) > 0:
            target_rejected = np.concatenate([target_rejected, removed_data])
        if flight_removed_count > 0:
            target_rejected = np.concatenate([target_rejected, flight_removed_data]) if len(target_rejected) > 0 else flight_removed_data
        if obsmode_removed_count > 0:
            target_rejected = np.concatenate([target_rejected, obsmode_removed_data]) if len(target_rejected) > 0 else obsmode_removed_data
        del removed_data
        clean_combined      = np.concatenate([other_data, target_clean]) if len(other_data) > 0 else target_clean
        clean_combined_orig = np.concatenate([other_orig_idx, target_clean_orig]) if len(other_data) > 0 else target_clean_orig
        all_rejected        = target_rejected
        del other_data, other_orig_idx, target_data, target_orig_idx, target_clean, target_clean_orig, target_rejected
    
    # Filter out all-zero spectra if requested
    zero_rejected_count = 0
    if filter_zero_spectra:
        zero_mask = _zero_mask_chunked(clean_combined['SPECTRUM'])

        if np.any(zero_mask):
            zero_rejected_count = int(np.sum(zero_mask))
            zero_spectra = clean_combined[zero_mask]
            all_rejected = np.concatenate([all_rejected, zero_spectra]) if len(all_rejected) > 0 else zero_spectra
            clean_combined      = clean_combined[~zero_mask]
            clean_combined_orig = clean_combined_orig[~zero_mask]

    # Merge param_rejected_data and peaks_rejected_data into all_rejected
    if param_rejected_count > 0:
        all_rejected = np.concatenate([all_rejected, param_rejected_data]) if len(all_rejected) > 0 else param_rejected_data
    if peaks_rejected_count > 0:
        all_rejected = np.concatenate([all_rejected, peaks_rejected_data]) if len(all_rejected) > 0 else peaks_rejected_data
    if tau_rejected_science_count + tau_rejected_tau_count > 0:
        all_rejected = np.concatenate([all_rejected, tau_rejected_data]) if len(all_rejected) > 0 else tau_rejected_data

    # Remap TSYS_INDEX / TAU_SIG_INDEX to match the new row positions in the
    # output file.  clean_combined_orig[i] is the original row number that
    # landed at position i; we invert that to get old→new.
    has_tsys = 'TSYS_INDEX'    in clean_combined.dtype.names
    has_tau  = 'TAU_SIG_INDEX' in clean_combined.dtype.names
    if has_tsys or has_tau:
        row_mapping = {int(old): new for new, old in enumerate(clean_combined_orig)}
        if has_tsys:
            _remap_index_column(clean_combined, 'TSYS_INDEX',    row_mapping)
        if has_tau:
            _remap_index_column(clean_combined, 'TAU_SIG_INDEX', row_mapping)

    # Set default output paths
    if output_clean is None:
        output_clean = Path("clean_data.fits")
    else:
        output_clean = Path(output_clean)
    
    write_rejected = output_rejected is not None
    if output_rejected is not None:
        output_rejected = Path(output_rejected)
    
    # Create FITS files with same structure as original
    # Copy primary HDU
    primary_hdu = hdul[0].copy()

    # Build binary table HDU directly from the numpy recarray (no Table copy)
    hdu_clean = fits.BinTableHDU(clean_combined)
    hdu_clean.name = matrix_hdus[0][1].name  # Use name from first spectrum HDU
    
    # Copy header information from original
    for key in matrix_hdus[0][1].header:
        if key not in ['NAXIS1', 'NAXIS2', 'TFIELDS'] and key != '':
            try:
                hdu_clean.header[key] = matrix_hdus[0][1].header[key]
            except (ValueError, KeyError):
                pass
    
    # Create FITS file for clean data
    n_clean = len(clean_combined)
    hdul_clean = fits.HDUList([primary_hdu, hdu_clean])
    hdul_clean.writeto(output_clean, overwrite=True)
    del clean_combined, hdu_clean, hdul_clean

    # Write rejected file only if a path was provided
    if write_rejected:
        if len(all_rejected) > 0:
            hdu_rejected = fits.BinTableHDU(all_rejected)
            hdu_rejected.name = matrix_hdus[0][1].name

            for key in matrix_hdus[0][1].header:
                if key not in ['NAXIS1', 'NAXIS2', 'TFIELDS'] and key != '':
                    try:
                        hdu_rejected.header[key] = matrix_hdus[0][1].header[key]
                    except (ValueError, KeyError):
                        pass

            hdul_rejected = fits.HDUList([primary_hdu, hdu_rejected])
            hdul_rejected.writeto(output_rejected, overwrite=True)
        else:
            primary_rejected = fits.PrimaryHDU()
            hdul_rejected = fits.HDUList([primary_rejected])
            hdul_rejected.writeto(output_rejected, overwrite=True)
    
    # Return paths and statistics
    stats = {
        'clean_count': n_clean,
        'rejected_total': len(all_rejected),
        'rejected_nan': nan_rejected_count,
        'rejected_zero': zero_rejected_count,
        'rejected_removed': removed_count,
        'rejected_flights': flight_removed_count,
        'flight_removed_missions': flight_removed_missions,
        'rejected_param': param_rejected_count,
        'param_filter_details': param_filter_details,
        'rejected_peaks': peaks_rejected_count,
        'rejected_tau_science': tau_rejected_science_count,
        'rejected_tau_tau': tau_rejected_tau_count,
        'tau_filter_details': tau_filter_details,
    }
    
    return output_clean, output_rejected, stats


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

    # Determine which rows are calibration spectra (must not be baselined)
    calibration_objects = {'TSYS', 'TAU_SIG'}
    if 'OBJECT' in data.dtype.names:
        row_objects = np.array([str(s).strip() for s in data['OBJECT']])
    else:
        row_objects = np.array([''] * len(data))

    # Build new spectra array by applying baseline_subtract per row
    import logging as _logging
    _logger = _logging.getLogger(__name__)
    n_skipped = {obj: 0 for obj in calibration_objects}
    new_spectra = []
    for i, spec in enumerate(data[spectrum_column]):
        if row_objects[i] in calibration_objects:
            # Leave calibration spectra untouched — their absolute values are needed
            new_spectra.append(np.array(spec, dtype=float))
            n_skipped[row_objects[i]] += 1
        else:
            corrected = baseline_subtract(np.array(spec, dtype=float),
                                          order=order,
                                          window=window)
            new_spectra.append(corrected)

    for obj, count in n_skipped.items():
        if count > 0:
            _logger.info(f"Baseline subtraction: skipped {count} {obj} spectra "
                         f"(calibration rows — absolute values preserved)")

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


def _extract_spectral_params(hdul: fits.HDUList) -> dict:
    """
    Extract spectral axis parameters from FITS file.
    
    Extracts velocity reference, velocity spacing per channel, and reference pixel
    information needed to reconstruct the velocity axis.
    
    Parameters
    ----------
    hdul : astropy.io.fits.HDUList
        FITS HDUList with spectral data.
    
    Returns
    -------
    dict
        Dictionary with keys:
        - 'velo_ref': Reference velocity in m/s (from VELOCITY column)
        - 'deltav': Velocity spacing per channel in m/s (from DELTAV column)
        - 'crpix1_spec': Reference pixel index, 1-indexed (from CRPIX1 header)
        - 'nchans': Number of spectral channels (from SPECTRUM column)
    
    Notes
    -----
    Parameters are extracted from the first spectrum, assuming they're constant
    across all observations in the FITS file.
    """
    # Locate binary table HDU
    table = None
    table_hdu = None
    for hdu in hdul[1:]:
        if hasattr(hdu, 'data') and hdu.data is not None:
            if 'SPECTRUM' in hdu.data.dtype.names:
                table = hdu.data
                table_hdu = hdu
                break
    
    if table is None:
        raise ValueError("No SPECTRUM column found in FITS table")
    
    # Extract velocity parameters from first spectrum
    velo_ref = float(table['VELOCITY'][0]) if 'VELOCITY' in table.dtype.names else 0.0
    deltav = float(table['DELTAV'][0]) if 'DELTAV' in table.dtype.names else 1.0
    
    # Extract reference pixel from header
    crpix1_spec = 1.0  # Default
    if table_hdu is not None and 'CRPIX1' in table_hdu.header:
        crpix1_spec = float(table_hdu.header['CRPIX1'])
    elif len(hdul) > 0 and 'CRPIX1' in hdul[0].header:
        crpix1_spec = float(hdul[0].header['CRPIX1'])
    
    # Get number of spectral channels
    nchans = table['SPECTRUM'].shape[1]
    
    return {
        'velo_ref': velo_ref,
        'deltav': deltav,
        'crpix1_spec': crpix1_spec,
        'nchans': nchans,
    }


def _create_velocity_axis(velo_ref: float, deltav: float, crpix1_spec: float,
                         nchans: int) -> np.ndarray:
    """
    Create velocity axis array for spectral data.
    
    Constructs a velocity axis based on FITS WCS spectral parameters.
    The velocity at each channel is calculated as:
        v[i] = velo_ref + (i - (crpix1_spec - 1)) * deltav
    
    Parameters
    ----------
    velo_ref : float
        Reference velocity in m/s
    deltav : float
        Velocity spacing per channel in m/s
    crpix1_spec : float
        Reference pixel (1-indexed, FITS convention)
    nchans : int
        Number of spectral channels
    
    Returns
    -------
    np.ndarray
        Velocity axis array of shape (nchans,) in m/s
    
    Examples
    --------
    >>> vel = _create_velocity_axis(470000.0, 500.0, 506.0, 1264)
    >>> print(f"Channel 0: {vel[0]:.1f} m/s")
    >>> print(f"Channel 505: {vel[505]:.1f} m/s")  # Reference pixel
    """
    channel_indices = np.arange(nchans, dtype=np.float64)
    velocity_axis = velo_ref + (channel_indices - (crpix1_spec - 1.0)) * deltav
    return velocity_axis


def _velocity_to_channel_index(velocity: float, velo_ref: float, deltav: float,
                               crpix1_spec: float) -> int:
    """
    Convert a velocity value to the nearest channel index.
    
    Inverse of the velocity axis formula:
        v = velo_ref + (i - (crpix1_spec - 1)) * deltav
    Solving for i:
        i = (v - velo_ref) / deltav + (crpix1_spec - 1)
    
    Parameters
    ----------
    velocity : float
        Velocity in m/s
    velo_ref : float
        Reference velocity in m/s
    deltav : float
        Velocity spacing per channel in m/s
    crpix1_spec : float
        Reference pixel (1-indexed, FITS convention)
    
    Returns
    -------
    int
        Channel index (0-indexed)
    """
    channel_float = (velocity - velo_ref) / deltav + (crpix1_spec - 1.0)
    channel_index = int(np.round(channel_float))
    return channel_index


def _extract_velocity_range(spectra: np.ndarray, velocity_axis: np.ndarray,
                           velo_min: float, velo_max: float) -> Tuple[np.ndarray, int, int]:
    """
    Extract spectra in a given velocity range.
    
    Parameters
    ----------
    spectra : np.ndarray
        Spectra array [n_spectra, n_channels]
    velocity_axis : np.ndarray
        Velocity axis [n_channels]
    velo_min : float
        Minimum velocity in m/s
    velo_max : float
        Maximum velocity in m/s
    
    Returns
    -------
    extracted_spectra : np.ndarray
        Spectra trimmed to velocity range
    ch_min : int
        Starting channel index
    ch_max : int
        Ending channel index (inclusive)
    
    Raises
    ------
    ValueError
        If velocity range is invalid or outside data range
    """
    # Find channels corresponding to velocity range
    mask = (velocity_axis >= velo_min) & (velocity_axis <= velo_max)
    indices = np.where(mask)[0]
    
    if len(indices) == 0:
        raise ValueError(
            f"No channels found in velocity range [{velo_min}, {velo_max}]. "
            f"Valid range: [{velocity_axis.min():.0f}, {velocity_axis.max():.0f}]"
        )
    
    ch_min = indices[0]
    ch_max = indices[-1]
    
    # Extract spectra
    extracted = spectra[:, ch_min:ch_max+1]
    
    return extracted, ch_min, ch_max


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

    # Load all spectra as float32 — half the memory of float64 while sufficient
    # for the reductions applied here.  Baseline fitting temporarily converts
    # individual chunks to float64 for numerical precision.
    spectra = np.asarray(data[spectrum_column], dtype=np.float32)

    # Track which rows to keep (for when we filter spectra)
    keep_mask = np.ones(len(spectra), dtype=bool)

    # Identify calibration rows (TSYS, TAU_SIG) — these must NOT be modified
    # by baseline subtraction, smoothing, or any other reduction step that
    # would corrupt their calibration values.
    NON_SCIENCE = {'TSYS', 'TAU_SIG'}
    if 'OBJECT' in data.dtype.names:
        raw_obj = data['OBJECT']
        if raw_obj.dtype.kind in ('S', 'U'):
            obj_strs = np.char.strip(raw_obj.astype('U'))
        else:
            obj_strs = np.array([str(v).strip() for v in raw_obj])
        science_mask = ~np.isin(obj_strs, list(NON_SCIENCE))
    else:
        science_mask = np.ones(len(spectra), dtype=bool)

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

    # Extract velocity range BEFORE baseline (so baseline window uses extracted channels)
    if 'extract' in methods:
        try:
            extract_range = methods['extract']
            extract_mode = methods.get('extract_mode', 'velocity')
            
            # Extract by velocity range (m/s)
            spectral_params = _extract_spectral_params(hdul)
            velocity_axis_full = _create_velocity_axis(
                spectral_params['velo_ref'],
                spectral_params['deltav'],
                spectral_params['crpix1_spec'],
                spectral_params['nchans']
            )
            
            velo_min, velo_max = extract_range[0], extract_range[1]
            spectra, ch_min, ch_max = _extract_velocity_range(
                spectra, velocity_axis_full, velo_min, velo_max
            )
            
            # Store channel range for later use with baseline window
            methods['_extract_ch_min'] = ch_min
            methods['_extract_ch_max'] = ch_max
            
            print(f"Extracted velocity range [{velo_min/1000:.0f}, {velo_max/1000:.0f}] km/s "
                  f"({velo_min:.0f}-{velo_max:.0f} m/s) "
                  f"→ channels [{ch_min}, {ch_max}] ({ch_max - ch_min + 1} channels)")
        except Exception as e:
            import warnings
            warnings.warn(f"Could not extract spectrum range: {e}", UserWarning)

    # Calculate RMS outside baseline window (after baseline subtraction)
    rms_baseline_values = None
    baseline_window_info = None
    
    if 'baseline' in methods:
        params = methods['baseline']
        order = params.get('order', 1)
        window = None
        
        # Handle window parameter - could be in velocity (m/s) or channel format
        if 'window_m_s' in params:
            # Window is in m/s (from config, absolute velocities)
            velo_min_window, velo_max_window = params['window_m_s']
            
            # Get full velocity axis to map velocities to channels
            spectral_params = _extract_spectral_params(hdul)
            velocity_axis_full = _create_velocity_axis(
                spectral_params['velo_ref'],
                spectral_params['deltav'],
                spectral_params['crpix1_spec'],
                spectral_params['nchans']
            )
            
            # Find channels corresponding to baseline window velocities
            ch_min_window_idx = _velocity_to_channel_index(
                velo_min_window, spectral_params['velo_ref'], 
                spectral_params['deltav'], spectral_params['crpix1_spec']
            )
            ch_max_window_idx = _velocity_to_channel_index(
                velo_max_window, spectral_params['velo_ref'],
                spectral_params['deltav'], spectral_params['crpix1_spec']
            )
            
            ch_min_window = int(np.clip(ch_min_window_idx, 0, spectral_params['nchans'] - 1))
            ch_max_window = int(np.clip(ch_max_window_idx, 0, spectral_params['nchans'] - 1))
            
            # If we extracted, adjust window to extracted coordinate system
            if '_extract_ch_min' in methods:
                extract_ch_min = methods['_extract_ch_min']
                ch_min_window = max(0, ch_min_window - extract_ch_min)
                ch_max_window = max(0, ch_max_window - extract_ch_min)
            
            window = (ch_min_window, ch_max_window)
            print(f"Baseline window: [{velo_min_window/1000:.0f}, {velo_max_window/1000:.0f}] km/s "
                  f"→ channels [{ch_min_window}, {ch_max_window}]")
            
            # Store window info for RMS calculation after baseline subtraction
            baseline_window_info = (ch_min_window, ch_max_window)
        elif 'window' in params:
            window = params['window']
        
        if window is not None:
            # Final validation and clamping
            window = (
                max(0, window[0]),
                min(spectra.shape[1] - 1, window[1])
            )
        
        # Process baseline in chunks: convert each chunk to float64 for numerical
        # precision, store result back as float32.  Avoids a full float64 copy.
        _CHUNK = 1000
        sci_indices = np.where(science_mask)[0]
        if baseline_window_info is not None:
            ch_min_window, ch_max_window = baseline_window_info
            rms_baseline_values = np.full(spectra.shape[0], np.nan, dtype=np.float32)
            outside_mask = np.ones(spectra.shape[1], dtype=bool)
            outside_mask[ch_min_window:ch_max_window + 1] = False
        for start in range(0, len(sci_indices), _CHUNK):
            idx = sci_indices[start:start + _CHUNK]
            chunk = spectra[idx].astype(np.float64)
            spectra[idx] = _reduce_baseline(chunk, order=order, window=window).astype(np.float32)
            if baseline_window_info is not None:
                for j, i in enumerate(idx):
                    outside_channels = spectra[i][outside_mask]
                    n_valid = int(np.sum(~np.isnan(outside_channels)))
                    if n_valid >= 2:
                        rms_baseline_values[i] = np.nanstd(outside_channels.astype(np.float64))

    if 'smooth' in methods:
        params = methods['smooth']
        window_size = params.get('window_size', 5)
        _CHUNK = 1000
        sci_indices = np.where(science_mask)[0]
        for start in range(0, len(sci_indices), _CHUNK):
            idx = sci_indices[start:start + _CHUNK]
            spectra[idx] = _reduce_smooth(spectra[idx], window_size=window_size)

    if 'decimate' in methods:
        factor = methods['decimate'].get('factor', 1)
        if factor > 1:
            spectra = spectra[:, ::factor]
            methods['_decimate_factor'] = factor
            print(f"Decimated by factor {factor}: {spectra.shape[1]} channels remaining")

    # Convert spectra back to original data type to preserve file size
    spectra = spectra.astype(original_spectrum_dtype)

    # Build output recarray — filter rows if unblank was used, then replace SPECTRUM.
    # Using a recarray copy avoids the extra memory of an astropy Table intermediate.
    filtered_data = data[keep_mask] if 'unblank' in methods else data[:]

    # If the number of channels changed (extraction / decimation) we cannot simply
    # copy the original recarray and overwrite the SPECTRUM column because the dtype
    # is fixed to the old channel count.  Rebuild the recarray with an updated dtype.
    original_nchans = filtered_data[spectrum_column].shape[1]
    new_nchans = spectra.shape[1]
    if new_nchans != original_nchans:
        new_dtype = []
        for name in filtered_data.dtype.names:
            if name == spectrum_column:
                new_dtype.append((name, original_spectrum_dtype, (new_nchans,)))
            else:
                new_dtype.append((name, filtered_data.dtype[name]))
        output_data = np.recarray(len(filtered_data), dtype=new_dtype)
        for name in filtered_data.dtype.names:
            if name != spectrum_column:
                output_data[name] = filtered_data[name]
        output_data[spectrum_column] = spectra
    else:
        output_data = filtered_data.copy()
        output_data[spectrum_column] = spectra
    del spectra  # free processing array before building the HDU

    # Update DELTAV and FREQRES if decimation was applied
    if '_decimate_factor' in methods:
        factor = methods['_decimate_factor']
        if 'DELTAV' in output_data.dtype.names:
            output_data['DELTAV'] = output_data['DELTAV'] * factor
        if 'FREQRES' in output_data.dtype.names:
            output_data['FREQRES'] = output_data['FREQRES'] * factor

    # Add RMS_BASELINE column if it was calculated during baseline subtraction.
    # Velocity axis is NOT stored as a column — it is always reconstructable from
    # VELOCITY, DELTAV, and CRPIX1 via basic_io.reconstruct_velocity_axis().
    if rms_baseline_values is not None:
        from numpy.lib.recfunctions import append_fields
        output_data = append_fields(output_data, 'RMS_BASELINE',
                                    rms_baseline_values, dtypes=np.float32,
                                    usemask=False, asrecarray=True)

    # Build new HDUList
    primary = hdul[0].copy()
    new_hdu = fits.BinTableHDU(output_data)
    del output_data
    new_hdu.name = matrix_hdu.name

    # Copy header keywords from original
    for key in matrix_hdu.header:
        if key not in ['NAXIS1', 'NAXIS2', 'TFIELDS'] and key != '':
            try:
                new_hdu.header[key] = matrix_hdu.header[key]
            except (ValueError, KeyError):
                pass

    # Recalculate CRPIX1 to account for channel extraction and/or decimation.
    # Original CRPIX1 is 1-indexed; convert to 0-based, apply shifts, convert back.
    if '_extract_ch_min' in methods or '_decimate_factor' in methods:
        old_crpix1 = float(new_hdu.header.get('CRPIX1', 1.0))
        # Convert to 0-based index
        ref_ch = old_crpix1 - 1.0
        # Shift for channel extraction: reference channel moves left by ch_min
        if '_extract_ch_min' in methods:
            ref_ch -= methods['_extract_ch_min']
        # Scale for decimation: reference channel index divides by factor
        if '_decimate_factor' in methods:
            ref_ch /= methods['_decimate_factor']
        # Convert back to 1-based FITS pixel
        new_hdu.header['CRPIX1'] = ref_ch + 1.0

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

    # Build methods dict from config if not provided
    if methods is None:
        methods = {}
        reduction_config = cfg.get('reduction', {})
        
        # Parse baseline parameters
        # window: [velo_min, velo_max] in km/s (absolute velocities)
        # Example: window=[450, 500] → 450,000 to 500,000 m/s
        if 'baseline' in reduction_config and reduction_config['baseline'] is not None:
            baseline_order = reduction_config['baseline']
            baseline_window = reduction_config.get('line_window', reduction_config.get('window', None))
            methods['baseline'] = {'order': baseline_order}
            if baseline_window is not None:
                # Store as-is; will be converted to channel indices in reduce_spectra
                # Store original values in km/s for reference
                methods['baseline']['window_km_s'] = tuple(baseline_window)
                # Also convert to m/s for internal use
                window_ms = [baseline_window[0] * 1000.0, baseline_window[1] * 1000.0]
                methods['baseline']['window_m_s'] = window_ms
        
        # Parse extraction parameters
        # extract: [velo_min, velo_max] in km/s (will be converted to m/s)
        # Example: extract=[350, 700] → 350,000 to 700,000 m/s
        if 'extract' in reduction_config and reduction_config['extract'] is not None:
            extract_range = reduction_config['extract']
            if len(extract_range) == 2:
                # Convert from km/s to m/s
                velo_min_ms = extract_range[0] * 1000.0
                velo_max_ms = extract_range[1] * 1000.0
                methods['extract'] = [velo_min_ms, velo_max_ms]
                methods['extract_mode'] = 'velocity'
        
        # Parse smooth parameters if present
        if 'smooth' in reduction_config and reduction_config['smooth'] is not None:
            smooth_window = int(reduction_config['smooth'])
            if smooth_window > 1:
                methods['smooth'] = {'window_size': smooth_window}
        
        # Parse unblank parameters if present
        if 'unblank' in reduction_config and reduction_config['unblank'] is not None:
            methods['unblank'] = reduction_config['unblank']

    # Apply reduction
    reduced_hdul = reduce_spectra(hdul, methods=methods)

    # Write output (always overwrite)
    reduced_hdul.writeto(output_path, overwrite=True)

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
