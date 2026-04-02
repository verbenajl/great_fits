"""
Input/output utilities for reading FITS files and configuration.
"""

import os
from pathlib import Path
from typing import Optional, Union

import numpy as np

try:
    import tomllib  # Built-in for Python 3.11+
except ImportError:
    import tomli as tomllib  # Fallback for Python < 3.11

from astropy.io import fits


def reconstruct_velocity_axis(hdu) -> np.ndarray:
    """
    Reconstruct the velocity axis (m/s) from FITS header/column parameters.

    Uses the VELOCITY column (reference velocity at CRPIX1), DELTAV column
    (channel spacing), and CRPIX1 header key (reference pixel, 1-indexed).
    These values are always present and already account for any channel
    extraction or decimation applied by ``reduce_spectra``.

    Parameters
    ----------
    hdu : astropy.io.fits.BinTableHDU
        Binary table HDU containing SPECTRUM, VELOCITY, and DELTAV columns.

    Returns
    -------
    np.ndarray
        Velocity axis in m/s, shape ``(n_channels,)``.
    """
    data     = hdu.data
    velo_ref = float(data['VELOCITY'][0])
    deltav   = float(data['DELTAV'][0])
    crpix1   = float(hdu.header.get('CRPIX1', 1.0))
    nchans   = data['SPECTRUM'].shape[1]
    ch       = np.arange(nchans, dtype=np.float64)
    return velo_ref + (ch - (crpix1 - 1.0)) * deltav


def get_config(config_path: Optional[Union[str, Path]] = None) -> dict:
    """
    Load configuration from a TOML file.
    
    Parameters
    ----------
    config_path : str or Path, optional
        Path to the config.toml file. If None, looks for config.toml 
        in the project root and parent directories.
    
    Returns
    -------
    dict
        Configuration dictionary.
    
    Raises
    ------
    FileNotFoundError
        If config file is not found.
    """
    if config_path is None:
        # Search for config.toml in project hierarchy
        current = Path.cwd()
        for candidate in [current / "config.toml"] + [p / "config.toml" for p in current.parents]:
            if candidate.exists():
                config_path = candidate
                break
        else:
            raise FileNotFoundError("config.toml not found in current or parent directories")
    else:
        config_path = Path(config_path)
        if not config_path.exists():
            raise FileNotFoundError(f"Config file not found: {config_path}")
    
    with open(config_path, "rb") as f:
        config = tomllib.load(f)
    
    return config


def read_fits_from_config(config_path: Optional[Union[str, Path]] = None) -> fits.HDUList:
    """
    Read a FITS file specified in the config.toml file.
    
    Looks for the 'fits_file' key in the [input] section of config.toml.
    
    Parameters
    ----------
    config_path : str or Path, optional
        Path to the config.toml file. If None, searches in project hierarchy.
    
    Returns
    -------
    astropy.io.fits.HDUList
        The opened FITS file.
    
    Raises
    ------
    FileNotFoundError
        If config file or FITS file is not found.
    KeyError
        If 'fits_file' key is not in the [input] section.
    
    Examples
    --------
    >>> hdul = read_fits_from_config()
    >>> hdul.info()
    
    >>> hdul = read_fits_from_config("./config.toml")
    >>> data = hdul[0].data
    """
    config = get_config(config_path)
    
    # Get the FITS file path from config
    try:
        fits_file = config["input"]["fits_file"]
    except KeyError as e:
        raise KeyError("'fits_file' not found in [input] section of config") from e
    
    fits_path = Path(fits_file)
    if not fits_path.exists():
        raise FileNotFoundError(f"FITS file not found: {fits_path}")
    
    # Open and return the FITS file
    hdul = fits.open(fits_path)
    return hdul


def read_fits(fits_file: Union[str, Path]) -> fits.HDUList:
    """
    Read a FITS file from a given path.
    
    Parameters
    ----------
    fits_file : str or Path
        Path to the FITS file.
    
    Returns
    -------
    astropy.io.fits.HDUList
        The opened FITS file.
    
    Raises
    ------
    FileNotFoundError
        If the FITS file is not found.
    
    Examples
    --------
    >>> hdul = read_fits("./data/cycle6.fits")
    >>> hdul.info()
    """
    fits_path = Path(fits_file)
    if not fits_path.exists():
        raise FileNotFoundError(f"FITS file not found: {fits_path}")
    
    return fits.open(fits_path)


def combine_fits_files(fits_file_list: list, output_hdul: fits.HDUList = None, 
                       single_hdu: bool = False) -> fits.HDUList:
    """
    Combine multiple FITS files into a single HDUList.
    
    This function takes a list of FITS file paths and combines them by
    either appending all HDUs into a single HDUList, or merging all data
    into a single binary table HDU.
    
    Parameters
    ----------
    fits_file_list : list
        List of paths to FITS files to combine.
    output_hdul : astropy.io.fits.HDUList, optional
        Existing HDUList to append to. If None, starts with the first FITS file.
    single_hdu : bool, optional
        If True, combines all spectral data into a single binary table HDU 
        (in addition to the primary HDU). If False (default), appends all HDUs 
        from all files into separate extension HDUs.
    
    Returns
    -------
    astropy.io.fits.HDUList
        Combined HDUList containing all HDUs from input files.
    
    Raises
    ------
    FileNotFoundError
        If any FITS file in the list is not found.
    ValueError
        If the input list is empty.
    
    Examples
    --------
    >>> files = ["/path/to/file1.fits", "/path/to/file2.fits"]
    >>> # Keep as separate HDUs
    >>> combined_hdul = combine_fits_files(files)
    >>> combined_hdul.writeto("/path/to/combined.fits", overwrite=True)
    
    >>> # Merge all data into single HDU
    >>> combined_hdul = combine_fits_files(files, single_hdu=True)
    >>> combined_hdul.writeto("/path/to/combined_single.fits", overwrite=True)
    """
    import numpy as np
    import tempfile
    import os

    if not fits_file_list:
        raise ValueError("fits_file_list cannot be empty")

    if single_hdu:
        # Streaming implementation: peak RAM ≈ size of largest single input file.
        #
        # Strategy:
        #   1. Scan all files with memmap=True → get dtype, total rows, headers.
        #      RAM: negligible (only headers and metadata).
        #   2. Allocate a disk-backed numpy.memmap of (total_rows,) with the
        #      correct dtype.  RAM: zero (backed by a temp file on disk).
        #   3. Fill the memmap one file at a time.
        #      RAM: one file's data at a time (≈ total / n_files).
        #   4. Build the output HDUList from the memmap and write to disk.
        #      RAM: one file's data at a time (astropy streams the memmap).
        #   5. Delete the temp file.

        # ── Phase 1: inventory ───────────────────────────────────────────────
        primary_hdu = None
        first_table_header = None
        dtype = None
        total_rows = 0

        for fits_file_path in fits_file_list:
            fits_path = Path(fits_file_path)
            if not fits_path.exists():
                raise FileNotFoundError(f"FITS file not found: {fits_path}")
            with fits.open(fits_path, memmap=True) as hdul:
                if primary_hdu is None:
                    primary_hdu = hdul[0].copy()
                    first_table_header = hdul[1].header.copy()
                    dtype = hdul[1].data.dtype
                for hdu in hdul[1:]:
                    if hdu.data is not None:
                        total_rows += len(hdu.data)

        if total_rows == 0:
            raise ValueError("No spectral data found in any FITS files")

        # ── Phase 2: disk-backed staging array ──────────────────────────────
        tmp_fd, tmp_path = tempfile.mkstemp(suffix='.memmap')
        os.close(tmp_fd)
        try:
            combined = np.memmap(tmp_path, dtype=dtype, mode='w+', shape=(total_rows,))

            # ── Phase 3: fill one file at a time ────────────────────────────
            offset = 0
            for fits_file_path in fits_file_list:
                with fits.open(fits_file_path, memmap=True) as hdul:
                    for hdu in hdul[1:]:
                        if hdu.data is not None:
                            n = len(hdu.data)
                            combined[offset:offset + n] = hdu.data[:]
                            offset += n
                combined.flush()

            # ── Phase 4: build output HDUList ────────────────────────────────
            combined_bintable_hdu = fits.BinTableHDU(combined)
            combined_bintable_hdu.name = "SPECTRA"

            # Copy WCS and spectral keywords from the first table header
            keys_to_copy = [
                'CRVAL2', 'CRVAL3', 'CRPIX2', 'CRPIX3',
                'CTYPE2', 'CTYPE3', 'CUNIT2', 'CUNIT3', 'CDELT2', 'CDELT3',
                'CRVAL1', 'CRPIX1', 'CTYPE1', 'CUNIT1', 'CDELT1',
                'VELOCITY', 'DELTAV', 'RESTFREQ', 'VELDEF',
            ]
            for key in keys_to_copy:
                if key in first_table_header:
                    primary_hdu.header[key] = first_table_header[key]
                    combined_bintable_hdu.header[key] = first_table_header[key]

            output_hdul = fits.HDUList([primary_hdu, combined_bintable_hdu])

        finally:
            # ── Phase 5: remove temp file ────────────────────────────────────
            try:
                os.unlink(tmp_path)
            except OSError:
                pass
    else:
        # Mode: Keep separate HDUs (original behavior)
        # Initialize with the first file's HDUList
        if output_hdul is None:
            first_file = Path(fits_file_list[0])
            if not first_file.exists():
                raise FileNotFoundError(f"FITS file not found: {first_file}")
            
            # Open, copy data to memory, and close immediately
            with fits.open(first_file) as hdul:
                # Create a new HDUList by copying the first file's HDUs into memory
                output_hdul = fits.HDUList([hdu.copy() for hdu in hdul])
            
            start_index = 1
        else:
            start_index = 0
        
        # Append HDUs from remaining files
        for fits_file_path in fits_file_list[start_index:]:
            fits_path = Path(fits_file_path)
            if not fits_path.exists():
                raise FileNotFoundError(f"FITS file not found: {fits_path}")
            
            with fits.open(fits_path) as hdul:
                # Append all HDUs except the primary HDU (which is usually redundant)
                # Copy each HDU to avoid reference issues with closed files
                for hdu in hdul[1:]:
                    output_hdul.append(hdu.copy())
    
    return output_hdul


def combine_fits_from_config(config_path: Optional[Union[str, Path]] = None) -> fits.HDUList:
    """
    Combine multiple FITS files from a configuration file into a single FITS file.
    
    Reads a list of FITS files from the [input] section 'fits_files' key in the
    config.toml file and combines them into a single HDUList. Writes the combined
    file to the output path specified in the [output] section 'combined_fits' key.
    
    Parameters
    ----------
    config_path : str or Path, optional
        Path to the config.toml file. If None, searches in project hierarchy.
    
    Returns
    -------
    astropy.io.fits.HDUList
        The combined HDUList.
    
    Raises
    ------
    FileNotFoundError
        If config file or any FITS file in the list is not found.
    KeyError
        If required keys are not in the config file.
    ValueError
        If the fits_files list is empty.
    
    Examples
    --------
    >>> combined = combine_fits_from_config()
    >>> combined.info()
    
    >>> combined = combine_fits_from_config("./combine_config.toml")
    """
    config = get_config(config_path)
    
    # Get the list of FITS files from config
    try:
        fits_files_list = config["input"]["fits_files"]
    except KeyError as e:
        raise KeyError("'fits_files' list not found in [input] section of config") from e
    
    if not isinstance(fits_files_list, list):
        raise ValueError("'fits_files' in config must be a list")
    
    if not fits_files_list:
        raise ValueError("'fits_files' list in config is empty")
    
    # Get the output path from config
    try:
        output_path = config["output"]["combined_fits"]
    except KeyError as e:
        raise KeyError("'combined_fits' not found in [output] section of config") from e
    
    # Combine the FITS files
    combined_hdul = combine_fits_files(fits_files_list)
    
    # Write to output file
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    combined_hdul.writeto(output_path, overwrite=True)
    
    print(f"Combined FITS file written to: {output_path}")
    
    return combined_hdul


def combine_fits_from_list(input_list: Union[str, Path], 
                           output_file: Union[str, Path],
                           single_hdu: bool = False) -> fits.HDUList:
    """
    Combine multiple FITS files from a text list file into a single FITS file.
    
    Reads a text file containing paths to FITS files (one per line, ignoring
    empty lines and comments starting with '#'), combines them into a single
    HDUList, and writes to the output file.
    
    Parameters
    ----------
    input_list : str or Path
        Path to a text file containing FITS file paths (one per line).
        Lines starting with '#' are treated as comments and ignored.
        Empty lines are ignored.
    output_file : str or Path
        Path where the combined FITS file will be written.
    single_hdu : bool, optional
        If True, combines all spectral data into a single binary table HDU.
        If False (default), keeps all HDUs as separate extensions.
    
    Returns
    -------
    astropy.io.fits.HDUList
        The combined HDUList.
    
    Raises
    ------
    FileNotFoundError
        If the input list file or any FITS file is not found.
    ValueError
        If the input list is empty or has no valid entries.
    
    Examples
    --------
    >>> # Combine with separate HDUs
    >>> combined = combine_fits_from_list("list.txt", "combined.fits")
    >>> combined.info()
    
    >>> # Combine into single HDU
    >>> combined = combine_fits_from_list("list.txt", "combined_single.fits", single_hdu=True)
    >>> combined.info()
    """
    input_path = Path(input_list)
    if not input_path.exists():
        raise FileNotFoundError(f"Input list file not found: {input_path}")
    
    # Read the list of FITS files from text file
    fits_files_list = []
    with open(input_path, 'r') as f:
        for line in f:
            # Strip whitespace and skip empty lines and comments
            line = line.strip()
            if line and not line.startswith('#'):
                fits_files_list.append(line)
    
    if not fits_files_list:
        raise ValueError(f"No valid FITS file paths found in {input_path}")
    
    # Combine the FITS files
    combined_hdul = combine_fits_files(fits_files_list, single_hdu=single_hdu)
    
    # Write to output file
    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    combined_hdul.writeto(output_path, overwrite=True)
    
    return combined_hdul
