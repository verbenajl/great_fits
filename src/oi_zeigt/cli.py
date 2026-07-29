"""
Command-line interface for OI ZEIGT.
"""

import sys
from pathlib import Path
from typing import Optional, Tuple

import click
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

from .basic_io import read_fits_from_config, read_fits, get_config, combine_fits_from_list, combine_fits_files
from .reduction.core import (analyze_spectrum_values, detect_blank_channels,
                            detect_nan_channels, filter_and_save_fits,
                            apply_baseline_from_config, reduce_spectra_from_config,
                            average_spectra_from_config, split_fits_by_mission)
from .statistics.quality import get_spechistogram, get_rmsratio_histogram, rmsratio_statistics


def _print_fits_details(hdul):
    """
    Helper function to print FITS file details.
    
    Parameters
    ----------
    hdul : astropy.io.fits.HDUList
        The FITS HDU list to print info about.
    """
    # Print basic info
    click.echo("\n" + "="*70)
    click.echo("FITS FILE INFORMATION")
    click.echo("="*70)
    
    # Try to print file info, but handle corrupted headers gracefully
    try:
        hdul.info()
    except (OSError, UnicodeDecodeError, ValueError) as e:
        click.echo(f"WARNING: Could not read full file info: {e}")
        click.echo("Attempting basic HDU listing...")
        # Fall back to manual HDU listing
        try:
            click.echo(f"Number of HDUs: {len(hdul)}")
            for i in range(len(hdul)):
                try:
                    hdu = hdul[i]
                    click.echo(f"  HDU {i}: {hdu.name if hasattr(hdu, 'name') else 'UNKNOWN'}")
                except (OSError, UnicodeDecodeError, ValueError):
                    click.echo(f"  HDU {i}: (Could not read header)")
        except Exception as e2:
            click.echo(f"Could not list HDUs: {e2}")
    
    click.echo("="*70 + "\n")
    
    # Print additional details
    click.echo(f"Number of HDUs: {len(hdul)}\n")
    
    for i in range(len(hdul)):
        try:
            hdu = hdul[i]
            click.echo(f"HDU {i}: {hdu.name} ({type(hdu).__name__})")
            
            # For BinTableHDU, get shape from header without loading data
            if isinstance(hdu, fits.BinTableHDU):
                nrows = hdu.header.get('NAXIS2', 0)
                ncols = hdu.header.get('TFIELDS', 0)
                click.echo(f"  Data shape: ({nrows}, {ncols})")
                click.echo(f"  Data type: BinTable (Large file - data not loaded)")
                
                # Show column names
                if hasattr(hdu, 'columns') and hdu.columns.names:
                    click.echo(f"  Columns: {', '.join(hdu.columns.names)}")
            elif hasattr(hdu, 'data'):
                try:
                    # Try to get shape safely
                    if hasattr(hdu.data, 'shape'):
                        click.echo(f"  Data shape: {hdu.data.shape}")
                        if hasattr(hdu.data, 'dtype'):
                            click.echo(f"  Data type: {hdu.data.dtype}")
                except (TypeError, MemoryError):
                    # File too large to load
                    if 'NAXIS2' in hdu.header:
                        click.echo(f"  Data shape: ({hdu.header['NAXIS2']}, {hdu.header.get('TFIELDS', '?')})")
                    click.echo(f"  Data type: (Large file - data not loaded)")
            
            if hdu.header:
                click.echo(f"  Header keywords: {len(hdu.header)}")
            click.echo()
        except (OSError, UnicodeDecodeError, ValueError) as e:
            click.echo(f"HDU {i}: (Could not read: {e})")
            click.echo()


def _print_velocity_resolution(hdul):
    """
    Print the velocity resolution (channel width) of every HDU that carries a
    spectral axis - both 3D cubes and CLASS-style spectral tables.

    Cubes are read from the WCS spectral axis (CTYPE3/CDELT3/CUNIT3): a velocity
    axis is reported directly, a frequency axis is converted via RESTFRQ. Spectral
    tables are read from the DELTAV column (m/s per channel); if DELTAV is not
    uniform across rows the min/max spread is reported as well.

    Parameters
    ----------
    hdul : astropy.io.fits.HDUList
        The FITS HDU list to analyze.
    """
    C_KMS = 299792.458  # speed of light, km/s

    def _cunit_to_kms(value, cunit):
        """Convert a velocity increment in the given CUNIT to km/s (or None)."""
        u = str(cunit).strip().lower()
        if u in ('m/s', 'm s-1', 'ms-1'):
            return value / 1000.0
        if u in ('km/s', 'km s-1', 'kms-1'):
            return value
        return None

    printed_header = False

    def _ensure_header():
        nonlocal printed_header
        if not printed_header:
            click.echo("="*70)
            click.echo("VELOCITY RESOLUTION")
            click.echo("="*70)
            printed_header = True

    for i, hdu in enumerate(hdul):
        header = getattr(hdu, 'header', None)
        if header is None:
            continue

        # --- Cubes / images with a spectral third axis ---
        if not isinstance(hdu, fits.BinTableHDU) and int(header.get('NAXIS', 0)) >= 3:
            ctype3 = str(header.get('CTYPE3', '')).strip().upper()
            cdelt3 = header.get('CDELT3', None)
            cunit3 = header.get('CUNIT3', '')
            if cdelt3 is not None:
                _ensure_header()
                label = f"HDU {i} ({hdu.name or type(hdu).__name__}) [cube, CTYPE3={ctype3 or '?'}]"
                if ctype3.startswith(('VELO', 'VRAD', 'VOPT', 'FELO')):
                    dv = _cunit_to_kms(abs(float(cdelt3)), cunit3)
                    if dv is not None:
                        click.echo(f"  {label}: {dv:.4f} km/s/channel")
                    else:
                        click.echo(f"  {label}: {abs(float(cdelt3)):.6g} {cunit3}/channel "
                                   f"(unrecognized CUNIT3, not converted)")
                elif ctype3.startswith('FREQ'):
                    restfrq = header.get('RESTFRQ', header.get('RESTFREQ', None))
                    df_hz = abs(float(cdelt3))  # CUNIT3 assumed Hz for FREQ axes
                    if restfrq:
                        dv = C_KMS * df_hz / float(restfrq)
                        click.echo(f"  {label}: {dv:.4f} km/s/channel "
                                   f"({df_hz/1e6:.4f} MHz, RESTFRQ={float(restfrq)/1e9:.4f} GHz)")
                    else:
                        click.echo(f"  {label}: {df_hz/1e6:.4f} MHz/channel "
                                   f"(no RESTFRQ - cannot convert to km/s)")
                else:
                    click.echo(f"  {label}: CDELT3={abs(float(cdelt3)):.6g} {cunit3} "
                               f"(unrecognized spectral CTYPE3)")
            continue

        # --- CLASS-style spectral tables (DELTAV column, m/s per channel) ---
        if isinstance(hdu, fits.BinTableHDU):
            colnames = hdu.columns.names if hasattr(hdu, 'columns') else []
            if 'DELTAV' not in colnames:
                continue
            try:
                deltav = np.asarray(hdu.data['DELTAV'], dtype=np.float64)
            except (OSError, ValueError) as e:
                _ensure_header()
                click.echo(f"  HDU {i} ({hdu.name or 'table'}): DELTAV present but unreadable ({e})")
                continue
            deltav = deltav[np.isfinite(deltav)]
            if deltav.size == 0:
                continue
            _ensure_header()
            dv_ref_kms = abs(float(deltav[0])) / 1000.0
            label = f"HDU {i} ({hdu.name or 'table'}) [spectral table]"
            spread = float(np.ptp(deltav))
            if spread > 1e-3 * abs(float(deltav[0])):
                click.echo(f"  {label}: {dv_ref_kms:.4f} km/s/channel (row 0); "
                           f"NON-UNIFORM: min {abs(deltav).min()/1000.0:.4f}, "
                           f"max {abs(deltav).max()/1000.0:.4f} km/s")
            else:
                click.echo(f"  {label}: {dv_ref_kms:.4f} km/s/channel")

    if printed_header:
        click.echo("="*70 + "\n")
    else:
        click.echo("="*70)
        click.echo("VELOCITY RESOLUTION")
        click.echo("="*70)
        click.echo("  No spectral axis (cube CDELT3 or table DELTAV) found.")
        click.echo("="*70 + "\n")


def _print_object_info(hdul, object_filter=None):
    """
    Print information about unique objects, AOR_IDs, and MISSION_IDs in the FITS file.
    
    Parameters
    ----------
    hdul : astropy.io.fits.HDUList
        The FITS HDU list to analyze.
    object_filter : str, optional
        If provided, only display entries matching this substring in the OBJECT column.
    """
    # Look for OBJECT column in binary tables
    for hdu in hdul:
        try:
            if isinstance(hdu, fits.BinTableHDU):
                # Check if 'OBJECT' column exists via header
                if 'OBJECT' in hdu.columns.names:
                    # Try to access the data with error handling
                    try:
                        objects = hdu.data['OBJECT']
                    except (OSError, ValueError):
                        # If normal access fails, try with memmap=False by re-opening
                        # For now, just warn the user
                        click.echo("="*70)
                        click.echo("OBJECT COLUMN INFORMATION")
                        click.echo("="*70)
                        click.echo(f"OBJECT column exists but could not be read (file may have encoding issues)")
                        click.echo("Consider regenerating the FITS file with proper encoding.\n")
                        return
                    
                    # Get unique objects and their counts
                    from collections import Counter
                    object_counts = Counter(objects)
                    
                    click.echo("="*70)
                    click.echo("UNIQUE OBJECTS IN FITS FILE")
                    click.echo("="*70)
                    click.echo(f"Total unique objects: {len(object_counts)}\n")
                    
                    # Sort by count (descending)
                    sorted_objects = sorted(object_counts.items(), key=lambda x: x[1], reverse=True)
                    
                    max_obj_len = max(len(obj.decode() if isinstance(obj, bytes) else obj) 
                                      for obj, _ in sorted_objects)
                    
                    for obj, count in sorted_objects:
                        obj_str = obj.decode().strip() if isinstance(obj, bytes) else str(obj).strip()
                        click.echo(f"  {obj_str:<{max_obj_len}}  : {count:6d} entries")
                    
                    click.echo("="*70 + "\n")
                
                # Filter data by object if object_filter is provided
                if object_filter and 'OBJECT' in hdu.columns.names:
                    try:
                        objects = hdu.data['OBJECT']
                        mask = np.array([
                            object_filter.upper() in (obj.decode().strip() if isinstance(obj, bytes) else str(obj).strip()).upper()
                            for obj in objects
                        ])
                        filtered_data = hdu.data[mask]
                        
                        click.echo("="*70)
                        click.echo(f"FILTERED DATA (OBJECT contains '{object_filter}')")
                        click.echo("="*70)
                        click.echo(f"Total entries matching filter: {np.sum(mask)}\n")
                        
                        # Show AOR_IDs for filtered data
                        if 'AOR_ID' in hdu.columns.names:
                            filtered_aor_ids = filtered_data['AOR_ID']
                            aor_id_counts = Counter(filtered_aor_ids)
                            
                            click.echo(f"AOR_IDs in filtered data ({len(aor_id_counts)} unique):")
                            sorted_aor_ids = sorted(aor_id_counts.items(), key=lambda x: x[1], reverse=True)
                            max_aor_len = max(len(aor_id.decode() if isinstance(aor_id, bytes) else aor_id) 
                                              for aor_id, _ in sorted_aor_ids)
                            
                            for aor_id, count in sorted_aor_ids:
                                aor_id_str = aor_id.decode().strip() if isinstance(aor_id, bytes) else str(aor_id).strip()
                                click.echo(f"  {aor_id_str:<{max_aor_len}}  : {count:6d} entries")
                        
                        # Show MISSION_IDs for filtered data
                        if 'MISSION_ID' in hdu.columns.names:
                            filtered_mission_ids = filtered_data['MISSION_ID']
                            mission_id_counts = Counter(filtered_mission_ids)

                            click.echo(f"\nMISSION_IDs in filtered data ({len(mission_id_counts)} unique):")
                            sorted_mission_ids = sorted(mission_id_counts.items(), key=lambda x: x[1], reverse=True)
                            max_mission_len = max(len(mission_id.decode() if isinstance(mission_id, bytes) else mission_id)
                                                  for mission_id, _ in sorted_mission_ids)

                            for mission_id, count in sorted_mission_ids:
                                mission_id_str = mission_id.decode().strip() if isinstance(mission_id, bytes) else str(mission_id).strip()
                                click.echo(f"  {mission_id_str:<{max_mission_len}}  : {count:6d} entries")

                        # Show TELESCOP for filtered data
                        if 'TELESCOP' in hdu.columns.names:
                            filtered_telescopes = filtered_data['TELESCOP']
                            telescop_counts = Counter(filtered_telescopes)

                            click.echo(f"\nBackends/TELESCOPs in filtered data ({len(telescop_counts)} unique):")
                            sorted_telescopes = sorted(telescop_counts.items(), key=lambda x: x[0])
                            max_tel_len = max(len(t.decode() if isinstance(t, bytes) else t)
                                              for t, _ in sorted_telescopes)

                            for tel, count in sorted_telescopes:
                                tel_str = tel.decode().strip() if isinstance(tel, bytes) else str(tel).strip()
                                click.echo(f"  {tel_str:<{max_tel_len}}  : {count:6d} entries")

                        click.echo("="*70 + "\n")
                        return  # Skip the general AOR/MISSION display if we're showing filtered data
                    except (OSError, ValueError, TypeError):
                        # If filtered data access fails, skip filtering info
                        pass
                
                # Check for AOR_ID column (only if not filtered)
                try:
                    if 'AOR_ID' in hdu.columns.names:
                        aor_ids = hdu.data['AOR_ID']
                        
                        # Get unique AOR_IDs and their counts
                        aor_id_counts = Counter(aor_ids)
                        
                        click.echo("="*70)
                        click.echo("UNIQUE AOR_IDS IN FITS FILE")
                        click.echo("="*70)
                        click.echo(f"Total unique AOR_IDs: {len(aor_id_counts)}\n")
                        
                        # Sort by count (descending)
                        sorted_aor_ids = sorted(aor_id_counts.items(), key=lambda x: x[1], reverse=True)
                        
                        max_aor_len = max(len(aor_id.decode() if isinstance(aor_id, bytes) else aor_id) 
                                          for aor_id, _ in sorted_aor_ids)
                        
                        for aor_id, count in sorted_aor_ids:
                            aor_id_str = aor_id.decode().strip() if isinstance(aor_id, bytes) else str(aor_id).strip()
                            click.echo(f"  {aor_id_str:<{max_aor_len}}  : {count:6d} entries")
                        
                        click.echo("="*70 + "\n")
                except (OSError, ValueError, TypeError):
                    pass
                
                # Check for MISSION_ID column (only if not filtered)
                try:
                    if 'MISSION_ID' in hdu.columns.names:
                        mission_ids = hdu.data['MISSION_ID']
                        scans = hdu.data['SCAN'] if 'SCAN' in hdu.columns.names else None

                        # Get unique MISSION_IDs and their counts
                        mission_id_counts = Counter(mission_ids)

                        click.echo("="*70)
                        click.echo("UNIQUE MISSION_IDS IN FITS FILE")
                        click.echo("="*70)
                        click.echo(f"Total unique MISSION_IDs: {len(mission_id_counts)}\n")

                        # Sort by count (descending)
                        sorted_mission_ids = sorted(mission_id_counts.items(), key=lambda x: x[1], reverse=True)

                        max_mission_len = max(len(mission_id.decode() if isinstance(mission_id, bytes) else mission_id)
                                              for mission_id, _ in sorted_mission_ids)

                        for mission_id, count in sorted_mission_ids:
                            mission_id_str = mission_id.decode().strip() if isinstance(mission_id, bytes) else str(mission_id).strip()
                            click.echo(f"  {mission_id_str:<{max_mission_len}}  : {count:6d} entries")
                            if scans is not None:
                                mid_mask = mission_ids == mission_id
                                unique_scans = sorted(set(int(s) for s in scans[mid_mask]))
                                click.echo(f"  {'':<{max_mission_len}}    scans: {unique_scans}")

                        click.echo("="*70 + "\n")
                except (OSError, ValueError, TypeError):
                    pass

                # Check for TELESCOP column (only if not filtered)
                try:
                    if 'TELESCOP' in hdu.columns.names:
                        telescopes = hdu.data['TELESCOP']

                        telescop_counts = Counter(telescopes)

                        click.echo("="*70)
                        click.echo("UNIQUE BACKENDS/TELESCOPES IN FITS FILE")
                        click.echo("="*70)
                        click.echo(f"Total unique backends: {len(telescop_counts)}\n")

                        sorted_telescopes = sorted(telescop_counts.items(), key=lambda x: x[0])

                        max_tel_len = max(len(t.decode() if isinstance(t, bytes) else t)
                                          for t, _ in sorted_telescopes)

                        for tel, count in sorted_telescopes:
                            tel_str = tel.decode().strip() if isinstance(tel, bytes) else str(tel).strip()
                            click.echo(f"  {tel_str:<{max_tel_len}}  : {count:6d} entries")

                        click.echo("="*70 + "\n")
                except (OSError, ValueError, TypeError):
                    pass

        except (OSError, UnicodeDecodeError, ValueError, TypeError, MemoryError) as e:
            click.echo(f"WARNING: Could not read object information from HDU: {e}")


@click.command()
@click.option(
    "--config",
    type=click.Path(exists=False),
    default=None,
    help="Path to config.toml file (reads fits_file from [input] section)"
)
@click.option(
    "--fits",
    type=click.Path(exists=False),
    default=None,
    help="Path to FITS file to read directly"
)
@click.option(
    "--reduced",
    is_flag=True,
    default=False,
    help="Print info from output.reduced_fits in config.toml"
)
@click.option(
    "--clean",
    is_flag=True,
    default=False,
    help="Print info from output.clean_fits in config.toml"
)
@click.option(
    "--prepared",
    is_flag=True,
    default=False,
    help="Print info from output.prepared_for_pca in config.toml"
)
@click.option(
    "--pcad",
    is_flag=True,
    default=False,
    help="Print info from output.pcad_fits (PCA-corrected) in config.toml"
)
@click.option(
    "--rejected",
    is_flag=True,
    default=False,
    help="Print info from output.rejected_fits in config.toml"
)
def print_fits_info(config: Optional[str], fits: Optional[str], reduced: bool, clean: bool, prepared: bool, pcad: bool, rejected: bool):
    """
    Print basic information about a FITS file.
    
    Can read the FITS file path from a config.toml file or directly specify it.
    
    Use --reduced, --clean, --prepared, --pcad, or --rejected flags to read from output paths in config.toml
    
    Examples:
    
        # Read FITS file from config.toml input section
        print_fits_info --config config.toml
        
        # Read FITS file from output.reduced_fits in config.toml
        print_fits_info --config config.toml --reduced
        
        # Read FITS file from output.clean_fits in config.toml
        print_fits_info --config config.toml --clean
        
        # Read FITS file from output.prepared_for_pca in config.toml
        print_fits_info --config config.toml --prepared
        
        # Read FITS file from output.pcad_fits (PCA-corrected) in config.toml
        print_fits_info --config config.toml --pcad
        
        # Read FITS file from output.rejected_fits in config.toml
        print_fits_info --config config.toml --rejected
        
        # Read FITS file directly
        print_fits_info --fits /path/to/file.fits
        
        # Use default config.toml from current/parent directory
        print_fits_info
    """
    try:
        # Load config to get object filter if available
        object_filter = None
        if config or (reduced or clean or prepared or pcad or rejected):
            try:
                import tomllib
            except ModuleNotFoundError:
                import tomli as tomllib
            
            config_path = config or "config.toml"
            with open(config_path, 'rb') as f:
                cfg = tomllib.load(f)
                parameters_cfg = cfg.get('parameters', {})
                object_filter = parameters_cfg.get('object', None)
                
                # Handle output file flags
                if reduced or clean or prepared or pcad or rejected:
                    output_cfg = cfg.get('output', {})
                    if reduced and 'reduced_fits' in output_cfg:
                        fits = output_cfg['reduced_fits']
                        click.echo(f"Reading from output.reduced_fits: {fits}")
                    elif clean and 'clean_fits' in output_cfg:
                        fits = output_cfg['clean_fits']
                        click.echo(f"Reading from output.clean_fits: {fits}")
                    elif prepared and 'prepared_for_pca' in output_cfg:
                        fits = output_cfg['prepared_for_pca']
                        click.echo(f"Reading from output.prepared_for_pca: {fits}")
                    elif pcad and 'pcad_fits' in output_cfg:
                        fits = output_cfg['pcad_fits']
                        click.echo(f"Reading from output.pcad_fits: {fits}")
                    elif rejected and 'rejected_fits' in output_cfg:
                        fits = output_cfg['rejected_fits']
                        click.echo(f"Reading from output.rejected_fits: {fits}")
                    elif reduced or clean or prepared or pcad or rejected:
                        click.echo(click.style(
                            f"Error: Requested output file not found in config.toml",
                            fg="red"
                        ), err=True)
                        sys.exit(1)
        
        # Read FITS file
        if fits:
            click.echo(f"Reading FITS file: {fits}")
            hdul = read_fits(fits)
        elif config:
            click.echo(f"Reading config from: {config}")
            hdul = read_fits_from_config(config)
        else:
            click.echo("Reading config from default location...")
            hdul = read_fits_from_config()
        
        _print_fits_details(hdul)
        _print_velocity_resolution(hdul)
        _print_object_info(hdul, object_filter=object_filter)
        hdul.close()
        
    except FileNotFoundError as e:
        click.echo(click.style(f"Error: {e}", fg="red"), err=True)
        sys.exit(1)
    except KeyError as e:
        click.echo(click.style(f"Error: {e}", fg="red"), err=True)
        sys.exit(1)


def _create_velocity_axis_from_fits(table_hdu: fits.BinTableHDU, nchans: int) -> Optional[np.ndarray]:
    """
    Create velocity axis on-the-fly from FITS header parameters if available.
    
    Attempts to extract VELOCITY, DELTAV, and CRPIX1 parameters from the FITS
    header and create a velocity axis. Returns None if parameters are missing.
    
    Parameters
    ----------
    table_hdu : fits.BinTableHDU
        FITS binary table HDU containing spectral data
    nchans : int
        Number of spectral channels
    
    Returns
    -------
    np.ndarray or None
        Velocity axis in m/s (shape: nchans), or None if parameters unavailable
    """
    try:
        # Try to get parameters from the table's data first
        if (hasattr(table_hdu, 'data') and table_hdu.data is not None and
            'VELOCITY' in table_hdu.data.dtype.names and 
            'DELTAV' in table_hdu.data.dtype.names):
            velo_ref = float(table_hdu.data['VELOCITY'][0])
            deltav = float(table_hdu.data['DELTAV'][0])
        else:
            return None
        
        # Get reference pixel from header
        crpix1_spec = 1.0
        if 'CRPIX1' in table_hdu.header:
            crpix1_spec = float(table_hdu.header['CRPIX1'])
        
        # Create velocity axis
        channel_indices = np.arange(nchans, dtype=np.float64)
        velocity_axis = velo_ref + (channel_indices - (crpix1_spec - 1.0)) * deltav
        
        return velocity_axis
        
    except (AttributeError, KeyError, TypeError, ValueError):
        return None


@click.command()
@click.option(
    "--config",
    type=click.Path(exists=False),
    default=None,
    help="Path to config.toml file"
)
@click.option(
    "--fits",
    type=click.Path(exists=False),
    default=None,
    help="Path to FITS file to read directly (default: input.fits_file from config)"
)
@click.option(
    "--reduced",
    is_flag=True,
    default=False,
    help="Plot from output.reduced_fits in config.toml"
)
@click.option(
    "--clean",
    is_flag=True,
    default=False,
    help="Plot from output.clean_fits in config.toml"
)
@click.option(
    "--prepared",
    is_flag=True,
    default=False,
    help="Plot from output.prepared_for_pca in config.toml"
)
@click.option(
    "--pcad",
    is_flag=True,
    default=False,
    help="Plot from output.pcad_fits in config.toml"
)
@click.option(
    "--post",
    is_flag=True,
    default=False,
    help="Plot from output.post_filtered_fits in config.toml"
)
@click.option(
    "--object",
    default=None,
    help="Object name to filter (substring match)"
)
@click.option(
    "--num-spectra",
    type=int,
    default=10,
    help="Number of spectra to plot (default: 10)"
)
@click.option(
    "--output",
    type=click.Path(),
    default=None,
    help="Output file for the plot (PNG or PDF). If not specified, show plot."
)
@click.option(
    "--filter-above",
    "filter_above_col",
    default=None,
    help="Column name: only plot spectra with value above --filter-value. "
         "If RMSRATIOB is requested but absent, falls back to RMSRATIO."
)
@click.option(
    "--filter-below",
    "filter_below_col",
    default=None,
    help="Column name: only plot spectra with value below --filter-value. "
         "If RMSRATIOB is requested but absent, falls back to RMSRATIO."
)
@click.option(
    "--filter-value",
    type=float,
    default=None,
    help="Threshold value for --filter-above or --filter-below."
)


def plot_sample_spectra(config: Optional[str], fits: Optional[str], reduced: bool, clean: bool, prepared: bool,
                       pcad: bool, post: bool, object: Optional[str], num_spectra: int, output: Optional[str],
                       filter_above_col: Optional[str], filter_below_col: Optional[str],
                       filter_value: Optional[float]):
    """
    Plot a sample of spectra from a FITS file.
    
    Optionally filter by object name and limit the number of spectra plotted.
    If VELOCITY_AXIS column is not present, attempts to create one on-the-fly
    from FITS header parameters (VELOCITY, DELTAV, CRPIX1).
    
    Examples:
    
        # Plot 20 spectra from the default FITS file in config.toml
        plot_sample_spectra --config config.toml
        
        # Plot only M51 spectra
        plot_sample_spectra --config config.toml --object M51
        
        # Plot 50 spectra and save as PDF
        plot_sample_spectra --fits /path/to/file.fits --num-spectra 50 --output plot.pdf
    """
    try:
        # Load config to get object filter and handle output file flags
        config_data = {}
        try:
            if config or (reduced or clean or prepared or pcad or post):
                try:
                    import tomllib
                except ModuleNotFoundError:
                    import tomli as tomllib
                
                config_path = config or "config.toml"
                with open(config_path, 'rb') as f:
                    cfg = tomllib.load(f)
                    config_data = cfg
                    
                    # Handle output file flags
                    if reduced or clean or prepared or pcad or post:
                        output_cfg = cfg.get('output', {})
                        if reduced and 'reduced_fits' in output_cfg:
                            fits = output_cfg['reduced_fits']
                            click.echo(f"Reading from output.reduced_fits: {fits}")
                        elif clean and 'clean_fits' in output_cfg:
                            fits = output_cfg['clean_fits']
                            click.echo(f"Reading from output.clean_fits: {fits}")
                        elif prepared and 'prepared_for_pca' in output_cfg:
                            fits = output_cfg['prepared_for_pca']
                            click.echo(f"Reading from output.prepared_for_pca: {fits}")
                        elif pcad and 'pcad_fits' in output_cfg:
                            fits = output_cfg['pcad_fits']
                            click.echo(f"Reading from output.pcad_fits: {fits}")
                        elif post and 'post_filtered_fits' in output_cfg:
                            fits = output_cfg['post_filtered_fits']
                            click.echo(f"Reading from output.post_filtered_fits: {fits}")
                        elif reduced or clean or prepared or pcad or post:
                            click.echo(click.style(
                                f"Error: Requested output file not found in config.toml",
                                fg="red"
                            ), err=True)
                            sys.exit(1)
            elif config:
                config_data = get_config(config)
            else:
                config_data = get_config()
        except FileNotFoundError as e:
            click.echo(f"Warning: Could not load config file: {e}")
            config_data = {}
        
        # Read FITS file
        if fits:
            hdul = read_fits(fits)
        elif config:
            hdul = read_fits_from_config(config)
        else:
            hdul = read_fits_from_config()
        
        # Get object name: priority is --object argument, then config file
        if object is None:
            try:
                object = config_data.get("parameters", {}).get("object")
            except (NameError, KeyError):
                pass
        
        if object is None:
            raise ValueError("Object name not specified. Use --object or set 'object' in [parameters] section of config file.")
        
        # Find binary table HDU with spectra
        matrix_hdu = None
        for hdu in hdul:
            if (hasattr(hdu, 'data') and hdu.data is not None
                    and hdu.data.dtype.names is not None
                    and 'SPECTRUM' in hdu.data.dtype.names):
                matrix_hdu = hdu
                break
        
        if matrix_hdu is None:
            raise ValueError("No HDU with SPECTRUM column found")
        
        data = matrix_hdu.data
        
        # Filter by object name (substring match)
        matching_indices = np.where([
            object.lower() in str(obj).lower()
            for obj in data['OBJECT']
        ])[0]
        
        if len(matching_indices) == 0:
            click.echo(click.style(f"No spectra found for object '{object}'", fg="red"), err=True)
            sys.exit(1)

        # Apply --filter-above / --filter-below with --filter-value
        filter_column = filter_above_col or filter_below_col
        filter_direction = 'above' if filter_above_col else ('below' if filter_below_col else None)

        if filter_column is not None and filter_value is not None:
            col = filter_column
            if col not in data.dtype.names:
                if col == 'RMSRATIOB' and 'RMSRATIO' in data.dtype.names:
                    click.echo(f"Warning: RMSRATIOB not found, falling back to RMSRATIO for filtering.")
                    col = 'RMSRATIO'
                else:
                    click.echo(click.style(f"Error: column '{col}' not found in FITS file.", fg="red"), err=True)
                    sys.exit(1)
            col_values = data[col][matching_indices]
            if filter_direction == 'above':
                mask = np.array([np.isfinite(float(v)) and float(v) > filter_value for v in col_values])
                click.echo(f"Filter: {col} > {filter_value} → {np.sum(mask)} spectra remaining")
            else:
                mask = np.array([np.isfinite(float(v)) and float(v) < filter_value for v in col_values])
                click.echo(f"Filter: {col} < {filter_value} → {np.sum(mask)} spectra remaining")
            matching_indices = matching_indices[mask]
            if len(matching_indices) == 0:
                click.echo(click.style(f"No spectra pass the filter.", fg="red"), err=True)
                sys.exit(1)

        # Sample spectra
        sample_size = min(num_spectra, len(matching_indices))
        sampled_indices = np.random.choice(matching_indices, size=sample_size, replace=False)
        sampled_indices = np.sort(sampled_indices)
        original_indices = sampled_indices.tolist()
        
        spectra_to_plot = [data[i] for i in sampled_indices]
        
        click.echo(f"Plotting {sample_size} spectra for object '{object}'")
        click.echo(f"Original FITS row indices (for reference): {original_indices}\n")
        
        # Reconstruct velocity axis from VELOCITY/DELTAV/CRPIX1 columns
        nchans = data['SPECTRUM'][0].shape[0]
        velocity_axis_from_fits = _create_velocity_axis_from_fits(matrix_hdu, nchans)
        if velocity_axis_from_fits is not None:
            click.echo("✓ Velocity axis reconstructed from FITS parameters")
        else:
            click.echo("  No velocity parameters found; using channel indices for x-axis")
        
        # Create plot
        num_to_plot = len(spectra_to_plot)
        cols = 4
        rows = int(np.ceil(num_to_plot / cols))
        
        fig, axes = plt.subplots(rows, cols, figsize=(15, 3*rows))
        
        # Flatten axes array for easier iteration
        if rows == 1 and cols == 1:
            axes = np.array([[axes]])
        elif rows == 1 or cols == 1:
            axes = axes.reshape(rows, cols)
        
        # Plot each spectrum
        for plot_num, (ax, spectrum_data, orig_idx) in enumerate(zip(axes.flat, spectra_to_plot, original_indices)):
            spectrum = spectrum_data['SPECTRUM']
            obj_name = spectrum_data['OBJECT'].decode().strip() if isinstance(spectrum_data['OBJECT'], bytes) else str(spectrum_data['OBJECT']).strip()
            
            x_axis = None
            x_label = "Channel"
            if velocity_axis_from_fits is not None:
                x_axis = velocity_axis_from_fits / 1000.0  # Convert m/s to km/s
                x_label = "Velocity (km/s)"
            
            # Calculate NaN fraction for this spectrum
            nan_mask, nan_frac = detect_nan_channels(spectrum)
            
            # Plot with x-axis (velocity if available, channel index otherwise)
            if x_axis is not None:
                ax.plot(x_axis, spectrum, linewidth=0.8)
            else:
                ax.plot(spectrum, linewidth=0.8)
            
            # Color code the title based on NaN fraction
            if nan_frac == 0:
                title_color = 'green'
                nan_indicator = "✓"
            elif nan_frac < 0.1:
                title_color = 'blue'
                nan_indicator = "~"
            elif nan_frac < 0.3:
                title_color = 'orange'
                nan_indicator = "!"
            else:
                title_color = 'red'
                nan_indicator = "✗"
            
            title = f"{nan_indicator} Row {orig_idx}: {obj_name} ({nan_frac:.1%} NaN)"
            if filter_column is not None:
                col = filter_column if filter_column in data.dtype.names else (
                    'RMSRATIO' if 'RMSRATIO' in data.dtype.names else None)
                if col:
                    try:
                        title += f"\n{col}={float(spectrum_data[col]):.3f}"
                    except (KeyError, TypeError):
                        pass
            ax.set_title(title, fontsize=10, color=title_color, weight='bold')
            ax.set_xlabel(x_label)
            ax.set_ylabel("Intensity")
            ax.grid(True, alpha=0.3)
        
        # Hide unused subplots
        for idx in range(num_to_plot, len(axes.flat)):
            axes.flat[idx].set_visible(False)
        
        # Add a legend explaining the NaN indicators
        legend_text = (
            "Legend: ✓ = No NaNs (green) | ~ = <10% NaNs (blue) | "
            "! = 10-30% NaNs (orange) | ✗ = >30% NaNs (red)"
        )
        fig.text(0.5, 0.02, legend_text, ha='center', fontsize=9, 
                style='italic', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.3))
        
        plt.tight_layout(rect=[0, 0.04, 1, 1])
        
        # Save or show
        if output:
            plt.savefig(output, dpi=150, bbox_inches='tight')
            click.echo(click.style(f"✓ Plot saved to {output}", fg="green"))
        else:
            plt.show()
            click.echo(click.style("✓ Done", fg="green"))
        
        hdul.close()
        
    except FileNotFoundError as e:
        click.echo(click.style(f"Error: {e}", fg="red"), err=True)
        sys.exit(1)
    except ValueError as e:
        click.echo(click.style(f"Error: {e}", fg="red"), err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(click.style(f"Unexpected error: {e}", fg="red"), err=True)
        sys.exit(1)


@click.command()
@click.option(
    "--config",
    type=click.Path(exists=False),
    default=None,
    help="Path to config.toml file"
)
@click.option(
    "--fits",
    type=click.Path(exists=False),
    default=None,
    help="Path to FITS file to read directly (default: input.fits_file from config)"
)
@click.option(
    "--reduced",
    is_flag=True,
    default=False,
    help="Plot from output.reduced_fits in config.toml"
)
@click.option(
    "--clean",
    is_flag=True,
    default=False,
    help="Plot from output.clean_fits in config.toml"
)
@click.option(
    "--prepared",
    is_flag=True,
    default=False,
    help="Plot from output.prepared_for_pca in config.toml"
)
@click.option(
    "--num-spectra",
    type=int,
    default=20,
    help="Number of spectra to plot (default: 20)"
)
@click.option(
    "--output",
    type=click.Path(),
    default=None,
    help="Output file for the plot (PNG or PDF). If not specified, show plot."
)
def plot_skies(config: Optional[str], fits: Optional[str], reduced: bool, clean: bool, prepared: bool,
              num_spectra: int, output: Optional[str]):
    """
    Plot a sample of sky/background spectra (SKYCHOPDIFF or SKY-DIFF observations).
    
    Examples:
    
        plot_skies --config config.toml
        plot_skies --config config.toml --num-spectra 30
        plot_skies --fits /path/to/file.fits --output sky_plot.pdf
    """
    try:
        # Handle output file flags
        if reduced or clean or prepared:
            try:
                import tomllib
            except ModuleNotFoundError:
                import tomli as tomllib
            
            config_path = config or "config.toml"
            with open(config_path, 'rb') as f:
                cfg = tomllib.load(f)
                output_cfg = cfg.get('output', {})
                if reduced and 'reduced_fits' in output_cfg:
                    fits = output_cfg['reduced_fits']
                    click.echo(f"Reading from output.reduced_fits: {fits}")
                elif clean and 'clean_fits' in output_cfg:
                    fits = output_cfg['clean_fits']
                    click.echo(f"Reading from output.clean_fits: {fits}")
                elif prepared and 'prepared_for_pca' in output_cfg:
                    fits = output_cfg['prepared_for_pca']
                    click.echo(f"Reading from output.prepared_for_pca: {fits}")
                elif reduced or clean or prepared:
                    click.echo(click.style(
                        f"Error: Requested output file not found in config.toml",
                        fg="red"
                    ), err=True)
                    sys.exit(1)
        
        # Read FITS file
        if fits:
            hdul = read_fits(fits)
        elif config:
            hdul = read_fits_from_config(config)
        else:
            hdul = read_fits_from_config()
        
        # Find binary table HDU with spectra
        matrix_hdu = None
        for hdu in hdul:
            if (hasattr(hdu, 'data') and hdu.data is not None
                    and hdu.data.dtype.names is not None
                    and 'SPECTRUM' in hdu.data.dtype.names):
                matrix_hdu = hdu
                break
        
        if matrix_hdu is None:
            raise ValueError("No HDU with SPECTRUM column found")
        
        data = matrix_hdu.data
        
        # Filter by OBJECT (look for SKYCHOPDIFF or SKY-DIFF)
        matching_indices = np.where([
            ("SKYCHOPDIFF" in str(obj).upper() or "SKY-DIFF" in str(obj).upper())
            for obj in data['OBJECT']
        ])[0]
        
        if len(matching_indices) == 0:
            click.echo(click.style("No sky spectra (SKYCHOPDIFF or SKY-DIFF) found", fg="red"), err=True)
            sys.exit(1)
        
        # Sample spectra
        sample_size = min(num_spectra, len(matching_indices))
        sampled_indices = np.random.choice(matching_indices, size=sample_size, replace=False)
        sampled_indices = np.sort(sampled_indices)
        original_indices = sampled_indices.tolist()
        
        spectra_to_plot = [data[i] for i in sampled_indices]
        
        click.echo(f"Plotting {sample_size} sky spectra")
        click.echo(f"Original FITS row indices (for reference): {original_indices}\n")
        
        # Create plot
        num_to_plot = len(spectra_to_plot)
        cols = 4
        rows = int(np.ceil(num_to_plot / cols))
        
        fig, axes = plt.subplots(rows, cols, figsize=(15, 3*rows))
        
        # Flatten axes array for easier iteration
        if rows == 1 and cols == 1:
            axes = np.array([[axes]])
        elif rows == 1 or cols == 1:
            axes = axes.reshape(rows, cols)
        
        # Plot each spectrum
        for plot_num, (ax, spectrum_data, orig_idx) in enumerate(zip(axes.flat, spectra_to_plot, original_indices)):
            spectrum = spectrum_data['SPECTRUM']
            obj_name = spectrum_data['OBJECT'].decode().strip() if isinstance(spectrum_data['OBJECT'], bytes) else str(spectrum_data['OBJECT']).strip()
            
            # Calculate NaN fraction for this spectrum
            nan_mask, nan_frac = detect_nan_channels(spectrum)
            
            ax.plot(spectrum, linewidth=0.8)
            
            # Color code the title based on NaN fraction
            if nan_frac == 0:
                title_color = 'green'
                nan_indicator = "✓"
            elif nan_frac < 0.1:
                title_color = 'blue'
                nan_indicator = "~"
            elif nan_frac < 0.3:
                title_color = 'orange'
                nan_indicator = "!"
            else:
                title_color = 'red'
                nan_indicator = "✗"
            
            ax.set_title(
                f"{nan_indicator} Row {orig_idx}: {obj_name} ({nan_frac:.1%} NaN)",
                fontsize=10,
                color=title_color,
                weight='bold'
            )
            ax.set_xlabel("Channel")
            ax.set_ylabel("Intensity")
            ax.grid(True, alpha=0.3)
        
        # Hide unused subplots
        for idx in range(num_to_plot, len(axes.flat)):
            axes.flat[idx].set_visible(False)
        
        # Add a legend explaining the NaN indicators
        legend_text = (
            "Legend: ✓ = No NaNs (green) | ~ = <10% NaNs (blue) | "
            "! = 10-30% NaNs (orange) | ✗ = >30% NaNs (red)"
        )
        fig.text(0.5, 0.02, legend_text, ha='center', fontsize=9, 
                style='italic', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.3))
        
        plt.tight_layout(rect=[0, 0.04, 1, 1])
        
        # Save or show
        if output:
            plt.savefig(output, dpi=150, bbox_inches='tight')
            click.echo(click.style(f"✓ Plot saved to {output}", fg="green"))
        else:
            plt.show()
            click.echo(click.style("✓ Done", fg="green"))
        
        hdul.close()
        
    except FileNotFoundError as e:
        click.echo(click.style(f"Error: {e}", fg="red"), err=True)
        sys.exit(1)
    except ValueError as e:
        click.echo(click.style(f"Error: {e}", fg="red"), err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(click.style(f"Unexpected error: {e}", fg="red"), err=True)
        sys.exit(1)


@click.command()
@click.option(
    "--config",
    type=click.Path(exists=False),
    default=None,
    help="Path to config.toml file"
)
@click.option(
    "--fits",
    type=click.Path(exists=False),
    default=None,
    help="Path to FITS file to read directly"
)
@click.option(
    "--reduced",
    is_flag=True,
    default=False,
    help="Read from output.reduced_fits in config"
)
@click.option(
    "--clean",
    is_flag=True,
    default=False,
    help="Read from output.clean_fits in config"
)
@click.option(
    "--prepared",
    is_flag=True,
    default=False,
    help="Read from output.prepared_for_pca in config"
)
@click.option(
    "--num-plots",
    type=int,
    default=8,
    help="Number of S-H_OBS/S-H_SKY comparisons to show (default: 8)"
)
@click.option(
    "--output",
    type=click.Path(),
    default=None,
    help="Output file for the plot (PNG or PDF). If not specified, show plot."
)
def plot_skyobsfit(config: Optional[str], fits: Optional[str], reduced: bool, clean: bool, 
                   prepared: bool, num_plots: int, output: Optional[str]):
    """
    Compare observed vs fitted sky spectra (S-H_OBS vs S-H_SKY/S-H_FIT).
    
    Displays side-by-side comparisons of observed sky spectra and their fitted 
    counterparts for visual inspection of how well the sky model fits the data.
    
    Examples:
    
        plot_skyobsfit --config config.toml
        plot_skyobsfit --config config.toml --num-plots 12
        plot_skyobsfit --fits /path/to/file.fits --output skyfit_comparison.pdf
    """
    try:
        # Handle output file flags
        if reduced or clean or prepared:
            try:
                import tomllib
            except ModuleNotFoundError:
                import tomli as tomllib
            
            config_path = config or "config.toml"
            with open(config_path, 'rb') as f:
                cfg = tomllib.load(f)
                output_cfg = cfg.get('output', {})
                if reduced and 'reduced_fits' in output_cfg:
                    fits = output_cfg['reduced_fits']
                    click.echo(f"Reading from output.reduced_fits: {fits}")
                elif clean and 'clean_fits' in output_cfg:
                    fits = output_cfg['clean_fits']
                    click.echo(f"Reading from output.clean_fits: {fits}")
                elif prepared and 'prepared_for_pca' in output_cfg:
                    fits = output_cfg['prepared_for_pca']
                    click.echo(f"Reading from output.prepared_for_pca: {fits}")
                elif reduced or clean or prepared:
                    click.echo(click.style(
                        f"Error: Requested output file not found in config.toml",
                        fg="red"
                    ), err=True)
                    sys.exit(1)
        
        # Read FITS file
        if fits:
            hdul = read_fits(fits)
        elif config:
            hdul = read_fits_from_config(config)
        else:
            hdul = read_fits_from_config()
        
        # Find binary table HDU with spectra
        matrix_hdu = None
        for hdu in hdul:
            if (hasattr(hdu, 'data') and hdu.data is not None
                    and hdu.data.dtype.names is not None
                    and 'SPECTRUM' in hdu.data.dtype.names):
                matrix_hdu = hdu
                break
        
        if matrix_hdu is None:
            raise ValueError("No HDU with SPECTRUM column found")
        
        data = matrix_hdu.data
        
        # Find S-H_OBS and S-H_SKY/S-H_FIT spectra
        obs_indices = []
        sky_indices = []
        
        for i, obj in enumerate(data['OBJECT']):
            obj_str = obj.decode().strip() if isinstance(obj, bytes) else str(obj).strip()
            if 'S-H_OBS' in obj_str:
                obs_indices.append(i)
            elif 'S-H_SKY' in obj_str or 'S-H_FIT' in obj_str:
                sky_indices.append(i)
        
        if len(obs_indices) == 0:
            click.echo(click.style("No S-H_OBS (observed sky) spectra found", fg="red"), err=True)
            sys.exit(1)
        
        if len(sky_indices) == 0:
            click.echo(click.style("No S-H_SKY or S-H_FIT (fitted sky) spectra found", fg="red"), err=True)
            sys.exit(1)
        
        # Create pairs: match by proximity or index order
        # Simple strategy: pair consecutive indices
        pairs = []
        for i in range(min(len(obs_indices), len(sky_indices))):
            pairs.append((obs_indices[i], sky_indices[i]))
        
        # Limit to requested number of plots
        num_to_show = min(num_plots, len(pairs))
        pairs = pairs[:num_to_show]
        
        if len(pairs) == 0:
            click.echo(click.style("Could not create any S-H_OBS/S-H_SKY pairs", fg="red"), err=True)
            sys.exit(1)
        
        click.echo(f"Plotting {len(pairs)} S-H_OBS vs S-H_SKY comparisons")
        click.echo(f"(Found {len(obs_indices)} observed and {len(sky_indices)} fitted sky spectra)\n")
        
        # Create plot with side-by-side comparisons
        fig, axes = plt.subplots(len(pairs), 2, figsize=(14, 3.5*len(pairs)))
        
        # Ensure axes is always 2D array
        if len(pairs) == 1:
            axes = axes.reshape(1, 2)
        
        # Plot each pair
        for pair_num, (obs_idx, sky_idx) in enumerate(pairs):
            obs_spectrum = data[obs_idx]['SPECTRUM']
            sky_spectrum = data[sky_idx]['SPECTRUM']
            
            obs_obj = data[obs_idx]['OBJECT'].decode().strip() if isinstance(data[obs_idx]['OBJECT'], bytes) else str(data[obs_idx]['OBJECT']).strip()
            sky_obj = data[sky_idx]['OBJECT'].decode().strip() if isinstance(data[sky_idx]['OBJECT'], bytes) else str(data[sky_idx]['OBJECT']).strip()
            
            obs_nan_mask, obs_nan_frac = detect_nan_channels(obs_spectrum)
            sky_nan_mask, sky_nan_frac = detect_nan_channels(sky_spectrum)
            
            # Left plot: observed sky
            ax_obs = axes[pair_num, 0]
            ax_obs.plot(obs_spectrum, linewidth=0.8, label='S-H_OBS', color='blue')
            
            # Color code by NaN fraction
            if obs_nan_frac == 0:
                title_color_obs = 'green'
                nan_indicator_obs = "✓"
            elif obs_nan_frac < 0.1:
                title_color_obs = 'blue'
                nan_indicator_obs = "~"
            elif obs_nan_frac < 0.3:
                title_color_obs = 'orange'
                nan_indicator_obs = "!"
            else:
                title_color_obs = 'red'
                nan_indicator_obs = "✗"
            
            ax_obs.set_title(
                f"{nan_indicator_obs} Row {obs_idx}: {obs_obj} ({obs_nan_frac:.1%} NaN)",
                fontsize=10,
                color=title_color_obs,
                weight='bold'
            )
            ax_obs.set_xlabel("Channel")
            ax_obs.set_ylabel("Intensity")
            ax_obs.grid(True, alpha=0.3)
            ax_obs.legend(loc='upper right', fontsize=9)
            
            # Right plot: fitted sky
            ax_sky = axes[pair_num, 1]
            ax_sky.plot(sky_spectrum, linewidth=0.8, label='S-H_SKY/FIT', color='orange')
            
            # Color code by NaN fraction
            if sky_nan_frac == 0:
                title_color_sky = 'green'
                nan_indicator_sky = "✓"
            elif sky_nan_frac < 0.1:
                title_color_sky = 'blue'
                nan_indicator_sky = "~"
            elif sky_nan_frac < 0.3:
                title_color_sky = 'orange'
                nan_indicator_sky = "!"
            else:
                title_color_sky = 'red'
                nan_indicator_sky = "✗"
            
            ax_sky.set_title(
                f"{nan_indicator_sky} Row {sky_idx}: {sky_obj} ({sky_nan_frac:.1%} NaN)",
                fontsize=10,
                color=title_color_sky,
                weight='bold'
            )
            ax_sky.set_xlabel("Channel")
            ax_sky.set_ylabel("Intensity")
            ax_sky.grid(True, alpha=0.3)
            ax_sky.legend(loc='upper right', fontsize=9)
        
        # Add legend explaining NaN indicators
        legend_text = (
            "Legend: ✓ = No NaNs (green) | ~ = <10% NaNs (blue) | "
            "! = 10-30% NaNs (orange) | ✗ = >30% NaNs (red)"
        )
        fig.text(0.5, 0.01, legend_text, ha='center', fontsize=9,
                style='italic', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.3))
        
        plt.tight_layout(rect=[0, 0.03, 1, 1])
        
        # Save or show
        if output:
            plt.savefig(output, dpi=150, bbox_inches='tight')
            click.echo(click.style(f"✓ Plot saved to {output}", fg="green"))
        else:
            plt.show()
            click.echo(click.style("✓ Done", fg="green"))
        
        hdul.close()
        
    except FileNotFoundError as e:
        click.echo(click.style(f"Error: {e}", fg="red"), err=True)
        sys.exit(1)
    except ValueError as e:
        click.echo(click.style(f"Error: {e}", fg="red"), err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(click.style(f"Unexpected error: {e}", fg="red"), err=True)
        sys.exit(1)


@click.command()
@click.option(
    "--config",
    type=click.Path(exists=False),
    default=None,
    help="Path to config.toml file"
)
@click.option(
    "--fits",
    type=click.Path(exists=False),
    default=None,
    help="Path to FITS file to read directly"
)
@click.option(
    "--sample-size",
    type=int,
    default=100,
    help="Number of spectra to sample for analysis (default: 100)"
)
def analyze_blanks(config: Optional[str], fits: Optional[str], sample_size: int):
    """
    Analyze spectrum values to identify blank/missing value markers.
    
    This command helps you understand what value is used for blanks in your data.
    
    Examples:
    
        analyze_blanks --config config.toml
        analyze_blanks --fits /path/to/file.fits --sample-size 50
    """
    try:
        # Read FITS file
        if fits:
            hdul = read_fits(fits)
        elif config:
            hdul = read_fits_from_config(config)
        else:
            hdul = read_fits_from_config()
        
        # Analyze
        analysis = analyze_spectrum_values(hdul, sample_size=sample_size)
        
        click.echo("\n" + "="*70)
        click.echo("SPECTRUM VALUE ANALYSIS")
        click.echo("="*70)
        click.echo(f"Sampled {analysis['total_values_analyzed']} values from {sample_size} spectra\n")
        
        click.echo("Statistics:")
        click.echo(f"  Min value:        {analysis['min_value']:.6e}")
        click.echo(f"  Max value:        {analysis['max_value']:.6e}")
        click.echo(f"  Mean value:       {analysis['mean_value']:.6e}")
        click.echo(f"  Median value:     {analysis['median_value']:.6e}")
        click.echo(f"  Std deviation:    {analysis['std_value']:.6e}\n")
        
        click.echo("Special counts:")
        click.echo(f"  Zero (0.0):       {analysis['zero_count']:10d} ({analysis['zero_count']/analysis['total_values_analyzed']*100:.2f}%)")
        click.echo(f"  Negative:         {analysis['negative_count']:10d} ({analysis['negative_count']/analysis['total_values_analyzed']*100:.2f}%)")
        click.echo(f"  Very small (<1e-10): {analysis['very_small_count']:10d} ({analysis['very_small_count']/analysis['total_values_analyzed']*100:.2f}%)\n")
        
        click.echo("Most common values:")
        for value, count in analysis['common_values'][:15]:
            percent = count / analysis['total_values_analyzed'] * 100
            click.echo(f"  {value:20.6e}  :  {count:8d}  ({percent:.2f}%)")
        
        click.echo("="*70 + "\n")
        
        hdul.close()
        
    except FileNotFoundError as e:
        click.echo(click.style(f"Error: {e}", fg="red"), err=True)
        sys.exit(1)
    except ValueError as e:
        click.echo(click.style(f"Error: {e}", fg="red"), err=True)
        sys.exit(1)


@click.command()
@click.option(
    "--config",
    type=click.Path(exists=True),
    default=None,
    help="Path to config.toml file"
)
@click.option(
    "--fits",
    type=click.Path(exists=True),
    default=None,
    help="Path to FITS file to read (default: input.fits_file from config.toml)"
)
@click.option(
    "--object",
    default=None,
    help="Object name to filter (substring match)"
)
@click.option(
    "--nan-threshold",
    type=float,
    default=None,
    help="NaN threshold for filtering (default: from config or 0.20 = 20%)"
)
@click.option(
    "--output-clean",
    type=click.Path(),
    default=None,
    help="Output path for clean FITS file (default: output.clean_fts)"
)
@click.option(
    "--output-rejected",
    type=click.Path(),
    default=None,
    help="Output path for rejected FITS file (default: output.rejected_fits)"
)
@click.option(
    "--remove",
    type=str,
    default=None,
    help="Column name to filter by (e.g., AOR_ID)"
)
@click.option(
    "--remove-values",
    type=str,
    multiple=True,
    help="Values to remove from the column (can be used multiple times). "
         "Example: --remove-values 04_0116_0010 --remove-values 04_0116_0020"
)
@click.option(
    "--apply-only-to-object",
    is_flag=True,
    default=False,
    help="Apply NaN filtering only to target object (not to all spectra)"
)
@click.option(
    "--filter-zero",
    is_flag=True,
    default=False,
    help="Filter out spectra that are all zeros (entire spectrum = 0)"
)
@click.option(
    "--filter-below",
    type=(str, float),
    multiple=True,
    metavar="COLUMN VALUE",
    help="Keep rows where COLUMN <= VALUE (filter out above). "
         "Rows with NaN in the column always pass. Can be repeated."
)
@click.option(
    "--filter-above",
    type=(str, float),
    multiple=True,
    metavar="COLUMN VALUE",
    help="Keep rows where COLUMN >= VALUE (filter out below). "
         "Rows with NaN in the column always pass. Can be repeated."
)
@click.option(
    "--filter-below-adaptive-rmsratio",
    "filter_below_adaptive_rmsratio",
    type=float,
    default=None,
    metavar="PERCENT",
    help="Adaptive per-file RMSRATIOB cut: instead of a fixed --filter-below "
         "RMSRATIOB value, keep the lowest-noise PERCENT% of the spectra. It uses "
         "the smallest RMSRATIOB threshold that still retains at least PERCENT% of "
         "the rows with a finite RMSRATIOB (i.e. the PERCENT-th percentile of "
         "RMSRATIOB, keeping rows <= that), dropping only the noisiest tail. "
         "Computed independently for each input file, so every flight gets its own "
         "cut — useful when a fixed threshold discards too many spectra on the "
         "noisier flights. Rows with NaN RMSRATIOB always pass. Example: "
         "--filter-below-adaptive-rmsratio 80 keeps ~80% of spectra."
)
@click.option(
    "--filter-below-adaptive-rmsratio-limit",
    "filter_below_adaptive_rmsratio_limit",
    type=float,
    default=None,
    metavar="VALUE",
    help="Hard cap on the --filter-below-adaptive-rmsratio threshold: the "
         "effective RMSRATIOB cut is min(adaptive percentile, VALUE), so spectra "
         "with RMSRATIOB > VALUE are always discarded even when a flight is so "
         "noisy that keeping PERCENT% would otherwise retain them (in that case "
         "fewer than PERCENT% are kept). Used together with "
         "--filter-below-adaptive-rmsratio; if given on its own it acts as a plain "
         "hard RMSRATIOB <= VALUE cut. Rows with NaN RMSRATIOB always pass. "
         "Example: --filter-below-adaptive-rmsratio 80 --filter-below-adaptive-rmsratio-limit 1.7."
)
@click.option(
    "--filter-spectrum-peaks",
    "spectrum_peak_threshold",
    type=float,
    default=None,
    metavar="VALUE",
    help="Reject science spectra that contain any channel with |value| > VALUE (e.g. 500)."
)
@click.option(
    "--filter-tau",
    "filter_tau",
    is_flag=True,
    default=False,
    help="Reject science spectra linked (via TAU_SIG_INDEX) to a TAU_SIG spectrum "
         "with any channel outside [0.001, 1.0]. Also removes the bad TAU_SIG rows."
)
@click.option(
    "--filter-flight",
    "filter_flight",
    type=str,
    multiple=True,
    help="Remove all rows whose MISSION_ID contains this string (e.g. F528). "
         "Can be specified multiple times to remove several flights."
)
@click.option(
    "--mission-id",
    "mission_id",
    type=str,
    multiple=True,
    help="Keep only rows whose MISSION_ID contains one of these substrings "
         "(comma-separated, e.g. F296,F528 — or repeat the flag). Applied "
         "before any other filter, so it also reduces memory/time use on "
         "large multi-mission files."
)
@click.option(
    "--aor-id",
    "aor_id",
    type=str,
    multiple=True,
    help="Keep only rows whose AOR_ID contains one of these substrings "
         "(comma-separated, e.g. 04_0116,04_0117 — or repeat the flag). "
         "Applied before any other filter, same as --mission-id."
)
@click.option(
    "--filter-object-exact",
    "filter_object_exact",
    type=str,
    default=None,
    help="Keep only rows whose OBJECT exactly matches one of these comma-separated values "
         "(e.g. --filter-object-exact M82,SKYDIFF)."
)
@click.option(
    "--filter-out-object-exact",
    "filter_out_object_exact",
    type=str,
    default=None,
    help="Remove all rows whose OBJECT exactly matches one of these comma-separated values "
         "(e.g. --filter-out-object-exact SKYCHOPDIFF,TSYS)."
)
@click.option(
    "--exclude-obsmode",
    "exclude_obsmode",
    type=str,
    default=None,
    help="Remove all rows whose OBSMODE exactly matches one of these comma-separated values "
         "(e.g. --exclude-obsmode BSAB,OTFSWB). TSYS and TAU_SIG rows are always kept."
)
@click.option(
    "--couple-tau-tsys",
    "couple_tau_tsys",
    is_flag=True,
    default=False,
    help="Add TSYS_INDEX and TAU_SIG_INDEX columns to the clean output if they are not already "
         "present (same pairing logic used by prepare_for_pca)."
)
@click.option(
    "--velocity-resample",
    "velocity_resample",
    type=float,
    default=None,
    metavar="KM_S",
    help="Resample the clean output spectra onto a single shared velocity grid at this "
         "channel spacing (km/s), via linear interpolation. Applied last, after all other "
         "filters. The grid is anchored to exact multiples of this value relative to 0 km/s "
         "(not to any one row's own reference velocity), so rows with slightly different "
         "native VELOCITY/DELTAV (e.g. different missions) end up on exactly the same grid — "
         "required for combine_fits/gridding to work correctly across them. Channels outside "
         "a row's original coverage become NaN. Does not apply to the rejected file."
)
@click.option(
    "--velocity-resample-range",
    "velocity_resample_range",
    type=(float, float),
    default=None,
    metavar="MIN MAX",
    help="Fix the resampled grid's extent (km/s) instead of deriving it from this file's "
         "own data. Required when running filter_fits separately per flight/mission before "
         "combining them later with combine_fits, since that requires every input file to "
         "have the same SPECTRUM column shape — without a fixed range, each flight's grid "
         "would only span whatever velocities that flight happens to cover. Only meaningful "
         "together with --velocity-resample."
)
@click.option(
    "--extract",
    is_flag=True,
    default=False,
    help="Extract (trim) the clean-output spectra to the velocity range from config "
         "[reduction].extract, exactly as reduce_spectra --extract does. Applied to the "
         "clean output only (the rejected file keeps its original width); every row is "
         "trimmed to the same channel window. Without this flag the spectrum length is "
         "left unchanged. If combined with --velocity-resample, extraction happens first."
)
def filter_fits(config: Optional[str], fits: Optional[str], object: Optional[str],
                nan_threshold: float, output_clean: Optional[str],
                output_rejected: Optional[str], remove: Optional[str],
                remove_values: tuple, apply_only_to_object: bool, filter_zero: bool,
                filter_below: tuple, filter_above: tuple,
                filter_below_adaptive_rmsratio: Optional[float],
                filter_below_adaptive_rmsratio_limit: Optional[float],
                spectrum_peak_threshold: Optional[float], filter_tau: bool,
                filter_flight: tuple, mission_id: tuple, aor_id: tuple, filter_object_exact: Optional[str],
                filter_out_object_exact: Optional[str],
                exclude_obsmode: Optional[str], couple_tau_tsys: bool,
                velocity_resample: Optional[float],
                velocity_resample_range: Optional[Tuple[float, float]],
                extract: bool):
    """
    Filter FITS data by object, NaN content, all-zero spectra, and/or column values.
    
    Creates two FITS files:
    1. Clean file: Spectra passing NaN threshold filter (and not all-zero if --filter-zero)
    2. Rejected file: Spectra failing NaN threshold filter or matching removal criteria
    
    Input file priority: --fits option > config.toml [input][fits_file]
    Output file defaults: From config.toml [output][clean_fits] and [output][rejected_fits]
    
    By default, NaN filtering is applied to ALL spectra. Use --apply-only-to-object
    to filter only the target object and keep all other objects regardless of NaN content.
    
    All-zero spectra filtering:
    - When --filter-zero is used, removes any spectrum where all channels are 0
    - Useful for identifying corrupted or non-existent observations
    - Applied to all spectra regardless of object filtering mode
    
    Examples:
    
        # Use config defaults for input and output files
        filter_fits --config config.toml
        
        # Override input file, use config defaults for output
        filter_fits --config config.toml --fits custom_input.fits
        
        # Filter only M51 spectra, keep all other objects
        filter_fits --config config.toml --apply-only-to-object
        
        # Filter all with custom thresholds
        filter_fits --config config.toml --nan-threshold 0.75
        
        # Filter out all-zero spectra (corrupted/empty data)
        filter_fits --config config.toml --filter-zero
        
        # Filter specific object and specify output files
        filter_fits --config config.toml --object "M51" \\
            --output-clean m51_clean.fits --output-rejected m51_rejected.fits
        
        # Remove specific AOR_ID values AND filter all spectra and zero spectra
        filter_fits --config config.toml --remove AOR_ID \\
            --remove-values 04_0116_0020609 --remove-values 04_0116_0020506 \\
            --filter-zero

        # Keep only spectra with RMSRATIOB <= 2.0
        filter_fits --config config.toml --filter-below RMSRATIOB 2.0

        # Keep only spectra with 1.3 <= RMSRATIOB <= 3.0
        filter_fits --config config.toml --filter-above RMSRATIOB 1.3 --filter-below RMSRATIOB 3.0
    """
    try:
        # Always try to load config
        try:
            if config:
                config_data = get_config(config)
            else:
                config_data = get_config()
        except FileNotFoundError:
            config_data = {}
        
        # Read FITS file - prioritize explicit --fits flag, then config input.fits_file
        if fits:
            click.echo(f"Reading FITS file: {fits}")
            hdul = read_fits(fits)
        else:
            # Try to get from config
            try:
                fits_from_config = config_data.get("input", {}).get("fits_file")
                if fits_from_config:
                    click.echo(f"Reading FITS file from config [input][fits_file]: {fits_from_config}")
                    hdul = read_fits(fits_from_config)
                else:
                    # Fallback to old behavior
                    click.echo(f"Reading config from: {config or 'default'}")
                    hdul = read_fits_from_config(config)
            except (FileNotFoundError, KeyError):
                # Fallback to old behavior
                click.echo(f"Reading config from: {config or 'default'}")
                hdul = read_fits_from_config(config)
        
        # Get object name from argument or config
        if object is None and config:
            try:
                object = config_data.get("parameters", {}).get("object")
            except (NameError, KeyError):
                pass
        
        if object is None and not fits:
            try:
                object = config_data.get("parameters", {}).get("object")
            except (NameError, KeyError):
                pass
        
        if object is None:
            raise ValueError("Object name not specified. Use --object or set in config file.")
        
        # Get NaN threshold from argument, config, or default
        if nan_threshold is None:
            if config:
                try:
                    nan_threshold = config_data.get("filters", {}).get("blank_fraction", 0.20)
                except (NameError, KeyError):
                    nan_threshold = 0.20
            else:
                nan_threshold = 0.20
        
        # Get output paths from config if not specified as arguments
        if output_clean is None and config:
            try:
                output_clean = config_data.get("output", {}).get("clean_fits")
            except (NameError, KeyError):
                pass
        
        # output_rejected is only written when explicitly provided via --output-rejected
        
        click.echo(f"\nFiltering data for object: {object}")
        click.echo(f"NaN threshold: {nan_threshold:.1%}")
        if apply_only_to_object:
            click.echo(f"Applying NaN filter to {object} spectra only")
        else:
            click.echo(f"Applying NaN filter to ALL spectra")
        click.echo(f"Keeping spectra with < {nan_threshold:.1%} NaN channels")
        click.echo(f"Rejecting spectra with >= {nan_threshold:.1%} NaN channels\n")
        
        # Check for zero filtering
        if filter_zero:
            click.echo(f"Also filtering out all-zero spectra\n")
        
        # Check for removal criteria from CLI args
        if remove and remove_values:
            click.echo(f"Also removing rows where {remove} = {', '.join(remove_values)}\n")
        
        # Check for removal criteria from config file
        if not remove and config:
            try:
                remove_aor_id = config_data.get("filters", {}).get("remove_aor_id")
                if remove_aor_id:
                    remove = "AOR_ID"
                    remove_values = (remove_aor_id,)  # Convert to tuple for consistency
                    click.echo(f"Also removing rows where AOR_ID = {remove_aor_id} (from config [filters][remove_aor_id])\n")
            except (NameError, KeyError):
                pass
        
        # Build param_filters from --filter-below / --filter-above
        param_filters = []
        for col, val in filter_below:
            param_filters.append((col, 'below', val))
            click.echo(f"Parameter filter: keep {col} <= {val}")
        for col, val in filter_above:
            param_filters.append((col, 'above', val))
            click.echo(f"Parameter filter: keep {col} >= {val}")

        # Adaptive RMSRATIOB cut: derive the --filter-below threshold from this
        # file's own RMSRATIOB distribution so a fixed number won't over-prune a
        # noisier flight. Keep the lowest-noise PERCENT% via the PERCENT-th
        # percentile of the finite RMSRATIOB values, optionally capped by a hard
        # limit so the effective 'below' threshold is min(percentile, limit).
        _adapt_pct = filter_below_adaptive_rmsratio
        _adapt_limit = filter_below_adaptive_rmsratio_limit
        if _adapt_pct is not None or _adapt_limit is not None:
            if _adapt_pct is not None and not (0 < _adapt_pct <= 100):
                click.echo(click.style(
                    "Error: --filter-below-adaptive-rmsratio must be a percentage in (0, 100]",
                    fg='red'), err=True)
                sys.exit(1)
            if _adapt_limit is not None and _adapt_limit <= 0:
                click.echo(click.style(
                    "Error: --filter-below-adaptive-rmsratio-limit must be positive",
                    fg='red'), err=True)
                sys.exit(1)
            rms_vals = None
            for _hdu in hdul:
                _d = getattr(_hdu, 'data', None)
                if _d is not None and getattr(_d, 'dtype', None) is not None \
                   and _d.dtype.names and 'RMSRATIOB' in _d.dtype.names:
                    rms_vals = np.asarray(_d['RMSRATIOB'], dtype=float)
                    # Restrict to the science rows the 'below' filter will act on:
                    # rows whose OBJECT contains the target object (same rule as
                    # filter_and_save_fits), so "keep PERCENT%" means PERCENT% of
                    # the science spectra, not of the whole (multi-object) table.
                    if object and 'OBJECT' in _d.dtype.names:
                        _objs = np.array([
                            (s.decode() if isinstance(s, bytes) else str(s)).strip().lower()
                            for s in _d['OBJECT']
                        ])
                        _tmask = np.char.find(_objs, object.lower()) >= 0
                        rms_vals = rms_vals[_tmask]
                    break
            if rms_vals is None:
                click.echo(click.style(
                    "Error: --filter-below-adaptive-rmsratio[-limit] requested but no "
                    "RMSRATIOB column was found in the input", fg='red'), err=True)
                sys.exit(1)
            finite = np.isfinite(rms_vals)
            n_finite = int(finite.sum())
            # Percentile part (skipped if RMSRATIOB is all-NaN or pct not given)
            pct_thr = None
            if _adapt_pct is not None:
                if n_finite == 0:
                    click.echo(click.style(
                        "Warning: --filter-below-adaptive-rmsratio requested but RMSRATIOB is "
                        "all-NaN — percentile not computed (all rows pass unless a limit is set)",
                        fg='yellow'), err=True)
                else:
                    pct_thr = float(np.percentile(rms_vals[finite], _adapt_pct))
            # Effective threshold = min(percentile, limit) over whichever are set
            candidates = [t for t in (pct_thr, _adapt_limit) if t is not None]
            if candidates:
                thr = min(candidates)
                param_filters.append(('RMSRATIOB', 'below', thr))
                if pct_thr is not None and _adapt_limit is not None:
                    binds = 'limit' if _adapt_limit <= pct_thr else f'P{_adapt_pct:g}'
                    desc = (f"keep lowest-noise {_adapt_pct:g}% (P{_adapt_pct:g}={pct_thr:.4f}) "
                            f"capped at limit {_adapt_limit:g} [{binds} binds]")
                elif pct_thr is not None:
                    desc = f"keep lowest-noise {_adapt_pct:g}% (P{_adapt_pct:g})"
                else:
                    desc = f"hard limit only"
                kept_frac = (f"{int(np.sum(rms_vals[finite] <= thr))}/{n_finite} = "
                             f"{100.0 * np.sum(rms_vals[finite] <= thr) / n_finite:.1f}% of finite-RMSRATIOB spectra kept"
                             if n_finite > 0 else "no finite-RMSRATIOB science spectra")
                click.echo(f"Adaptive RMSRATIOB filter: {desc} -> RMSRATIOB <= {thr:.4f} ({kept_frac})")
                if any(c == 'RMSRATIOB' for c, _v in filter_below):
                    click.echo(click.style(
                        "  Note: an explicit --filter-below RMSRATIOB is also set; both apply "
                        "(the stricter/lower threshold wins).", fg='yellow'))

        if spectrum_peak_threshold is not None:
            click.echo(f"Spectrum peak filter: reject science spectra with any |channel| > {spectrum_peak_threshold}")
        if filter_tau:
            click.echo("TAU filter: reject science spectra linked to TAU_SIG with channels outside [0.001, 1.0]")
        if filter_flight:
            click.echo(f"Flight filter: removing all rows with MISSION_ID containing: {', '.join(filter_flight)}")

        # Accept both repeated (--mission-id 296 --mission-id 298) and
        # comma-separated (--mission-id 296,298) forms.
        mission_id_list = [s.strip() for m in mission_id for s in m.split(',') if s.strip()]
        if mission_id_list:
            click.echo(f"Mission-ID filter: keeping only rows with MISSION_ID containing: {', '.join(mission_id_list)}")

        # Same dual-form parsing as --mission-id.
        aor_id_list = [s.strip() for a in aor_id for s in a.split(',') if s.strip()]
        if aor_id_list:
            click.echo(f"AOR-ID filter: keeping only rows with AOR_ID containing: {', '.join(aor_id_list)}")

        obj_exact_list = [s.strip() for s in filter_object_exact.split(',')] if filter_object_exact else None
        obj_out_exact_list = [s.strip() for s in filter_out_object_exact.split(',')] if filter_out_object_exact else None

        if obj_exact_list:
            click.echo(f"Object exact keep filter: keeping only OBJECT in: {', '.join(obj_exact_list)}")
        if obj_out_exact_list:
            click.echo(f"Object exact remove filter: removing rows with OBJECT in: {', '.join(obj_out_exact_list)}")

        exclude_obsmode_list = [s.strip() for s in exclude_obsmode.split(',')] if exclude_obsmode else None
        if exclude_obsmode_list:
            click.echo(f"OBSMODE exclusion filter: removing rows with OBSMODE in: {', '.join(exclude_obsmode_list)}")

        if velocity_resample is not None:
            click.echo(f"Velocity resampling: clean output will be resampled to {velocity_resample} km/s/channel")
            if velocity_resample_range is not None:
                click.echo(f"  Fixed grid range: {velocity_resample_range[0]} to {velocity_resample_range[1]} km/s "
                           f"(same for every run — required for combine_fits across separate runs)")
            else:
                click.echo("  WARNING: no --velocity-resample-range given — the grid extent will be "
                            "derived from this file's own data, which can differ run to run "
                            "(e.g. across separate per-flight runs) and break combine_fits later.")

        # Resolve --extract from config [reduction].extract (same as reduce_spectra:
        # flag reads the config range in km/s, converted to m/s here; error if the
        # flag is set but the config key is missing).
        extract_velocity_range_m_s = None
        if extract:
            extract_cfg = (config_data or {}).get('reduction', {}).get('extract', None)
            if extract_cfg is None:
                click.echo(click.style(
                    "Error: --extract specified but [reduction].extract is not defined in config",
                    fg='red'), err=True)
                sys.exit(1)
            try:
                extract_velocity_range_m_s = (float(extract_cfg[0]) * 1000.0,
                                              float(extract_cfg[1]) * 1000.0)
                click.echo(f"Extracting velocity range [{extract_cfg[0]}, {extract_cfg[1]}] km/s")
            except (ValueError, TypeError, IndexError):
                click.echo(click.style(
                    "Error: could not parse [reduction].extract from config (expected [min, max] km/s)",
                    fg='red'), err=True)
                sys.exit(1)

        # Filter and save
        clean_path, rejected_path, stats = filter_and_save_fits(
            hdul,
            object_name=object,
            nan_threshold=nan_threshold,
            output_clean=output_clean,
            output_rejected=output_rejected,
            remove_column=remove,
            remove_values=list(remove_values) if remove_values else None,
            apply_to_all=not apply_only_to_object,
            filter_zero_spectra=filter_zero,
            param_filters=param_filters if param_filters else None,
            spectrum_peak_threshold=spectrum_peak_threshold,
            filter_tau=filter_tau,
            filter_flights=list(filter_flight) if filter_flight else None,
            filter_mission_ids=mission_id_list if mission_id_list else None,
            filter_aor_ids=aor_id_list if aor_id_list else None,
            filter_object_exact=obj_exact_list,
            filter_out_object_exact=obj_out_exact_list,
            exclude_obsmode=exclude_obsmode_list,
            velocity_resample_km_s=velocity_resample,
            velocity_resample_range_km_s=velocity_resample_range,
            extract_velocity_range_m_s=extract_velocity_range_m_s,
        )
        
        # Optionally add TSYS_INDEX / TAU_SIG_INDEX to the clean output
        if couple_tau_tsys:
            from astropy.io import fits as fits_lib
            from .pca_analysis.prepare_for_pca import build_tau_tsys_indices
            with fits_lib.open(clean_path) as hdul_c:
                tbl = hdul_c[1]
                if 'TSYS_INDEX' in tbl.data.dtype.names and 'TAU_SIG_INDEX' in tbl.data.dtype.names:
                    click.echo("--couple-tau-tsys: TSYS_INDEX/TAU_SIG_INDEX already present, skipping.")
                else:
                    click.echo("--couple-tau-tsys: building TSYS_INDEX and TAU_SIG_INDEX columns...")
                    tsys_idx, tau_idx = build_tau_tsys_indices(tbl.data, object)
                    col_list = []
                    for name in tbl.data.dtype.names:
                        col_list.append(fits_lib.Column(
                            name=name,
                            format=tbl.columns[name].format,
                            array=tbl.data[name],
                        ))
                    col_list.append(fits_lib.Column(name='TSYS_INDEX',    format='J', array=tsys_idx))
                    col_list.append(fits_lib.Column(name='TAU_SIG_INDEX', format='J', array=tau_idx))
                    new_tbl = fits_lib.BinTableHDU.from_columns(col_list, header=tbl.header)
                    fits_lib.HDUList([fits_lib.PrimaryHDU(), new_tbl]).writeto(
                        clean_path, overwrite=True, checksum=False
                    )
                    n_sci = int(np.sum(tsys_idx >= 0))
                    click.echo(f"  Coupled {n_sci} science spectra to TSYS/TAU_SIG.")

        # Get statistics before closing
        from astropy.io import fits as fits_lib
        hdul_clean = fits_lib.open(clean_path)
        n_clean = len(hdul_clean[1].data) if len(hdul_clean) > 1 else 0
        hdul_clean.close()

        hdul.close()

        # Print results
        click.echo("="*70)
        click.echo(click.style("✓ Successfully created FITS files", fg="green"))
        click.echo("="*70)
        click.echo(f"\nClean FITS file: {clean_path}")
        click.echo(f"  Records: {n_clean}")
        click.echo(f"  (Science spectra passing all filters + calibration rows not rejected by --filter-tau)")

        n_rejected = 0
        if rejected_path is not None:
            hdul_rejected = fits_lib.open(rejected_path)
            n_rejected = len(hdul_rejected[1].data) if len(hdul_rejected) > 1 and hdul_rejected[1].data is not None else 0
            hdul_rejected.close()
            click.echo(f"\nRejected FITS file: {rejected_path}")
            click.echo(f"  Records: {n_rejected}")

        # Print detailed rejection statistics
        n_sci_rejected = (stats['rejected_nan'] + stats['rejected_zero'] + stats['rejected_removed']
                          + stats['rejected_param'] + stats.get('rejected_peaks', 0)
                          + stats.get('rejected_tau_science', 0))
        n_cal_rejected = stats.get('rejected_tau_tau', 0)
        click.echo(f"\n  Rejection breakdown (science = {object}, filters applied to science rows only):")
        click.echo(f"    - NaN threshold violations: {stats['rejected_nan']}")
        if filter_zero:
            click.echo(f"    - All-zero spectra: {stats['rejected_zero']}")
        if stats['rejected_removed'] > 0:
            click.echo(f"    - Removed by --remove criteria: {stats['rejected_removed']}")
        if stats.get('rejected_flights', 0) > 0:
            click.echo(f"    - Removed by --filter-flight: {stats['rejected_flights']} rows")
            for mid in stats.get('flight_removed_missions', []):
                click.echo(f"        {mid}")
        if stats['rejected_param'] > 0:
            click.echo(f"    - Removed by --filter-below/--filter-above: {stats['rejected_param']} total")
            for col, direction, value, n in stats.get('param_filter_details', []):
                op = '<=' if direction == 'below' else '>='
                click.echo(f"        kept {col} {op} {value}: {n} removed individually")
        if spectrum_peak_threshold is not None:
            click.echo(f"    - Removed by --filter-spectrum-peaks {spectrum_peak_threshold}: {stats.get('rejected_peaks', 0)}")
        if filter_tau:
            n_sci = stats.get('rejected_tau_science', 0)
            n_tau = stats.get('rejected_tau_tau', 0)
            click.echo(f"    - Removed by --filter-tau: {n_sci} science spectra, {n_tau} calibration (TAU_SIG) spectra")
            for orig_idx_val, n_sci_removed in stats.get('tau_filter_details', []):
                if n_sci_removed > 0:
                    click.echo(f"        TAU_SIG original index {orig_idx_val}: removed {n_sci_removed} science spectra")
        click.echo(f"    ----------------------------------------")
        click.echo(f"    Total science rejected : {n_sci_rejected}")
        click.echo(f"    Total calibration rejected: {n_cal_rejected}")
        click.echo(f"    Total rejected         : {n_sci_rejected + n_cal_rejected}  (file records: {n_rejected})")
        click.echo(f"    Clean records          : {n_clean}")
        click.echo(f"    Input total (est.)     : {n_clean + n_rejected}")

        click.echo()
        
    except FileNotFoundError as e:
        click.echo(click.style(f"Error: {e}", fg="red"), err=True)
        sys.exit(1)
    except ValueError as e:
        click.echo(click.style(f"Error: {e}", fg="red"), err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(click.style(f"Unexpected error: {e}", fg="red"), err=True)
        sys.exit(1)


@click.command()
@click.option(
    "--fits",
    type=click.Path(exists=True),
    required=True,
    help="Input FITS file to filter."
)
@click.option(
    "--output",
    type=click.Path(),
    required=True,
    help="Output FITS file (filtered)."
)
@click.option(
    "--yaml",
    "yaml_file",
    type=click.Path(exists=True),
    default=None,
    help="Path to mission parameters YAML file (default: mission_id_parameters.yml bundled with the package)."
)
def filter_missions(fits, output, yaml_file):
    """
    Drop spectra from a FITS file according to the drop rules in mission_id_parameters.yml.

    For each mission, three drop rules are supported:

    \b
      drop:
        telescope:            # drop ALL rows for these telescopes in this mission
          - LFAV_3
        scans:
          complete:           # drop ALL rows for these scans (any telescope)
            - 14381
          telescope:          # drop these scans only for the named telescope
            LFAV_0:
              - 18611

    Rows that match any rule are removed from the output; all other rows are kept unchanged.

    Examples:

        filter_missions --fits clean.fits --output clean_filtered.fits

        filter_missions --fits clean.fits --output clean_filtered.fits \\
            --yaml /path/to/custom_parameters.yml
    """
    import yaml as _yaml
    from pathlib import Path as _Path
    from astropy.io import fits as _fits

    # Resolve YAML file
    if yaml_file is None:
        yaml_file = _Path(__file__).parent / "pca_analysis" / "mission_id_parameters.yml"
    else:
        yaml_file = _Path(yaml_file)

    if not yaml_file.exists():
        click.echo(click.style(f"Error: YAML file not found: {yaml_file}", fg="red"), err=True)
        sys.exit(1)

    with open(yaml_file, "r") as f:
        mission_params = _yaml.safe_load(f) or {}

    click.echo(f"Loaded {len(mission_params)} mission entries from {yaml_file}")

    # Load FITS
    click.echo(f"Reading {fits} ...")
    with _fits.open(fits) as hdul:
        primary_hdu = hdul[0].copy()
        # Find table HDU
        table_hdu = None
        for hdu in hdul[1:]:
            if (hasattr(hdu, "data") and hdu.data is not None
                    and hasattr(hdu.data, "dtype")
                    and hdu.data.dtype.names is not None
                    and "MISSION_ID" in hdu.data.dtype.names):
                    table_hdu = hdu
                    break
        if table_hdu is None:
            click.echo(click.style("Error: no HDU with MISSION_ID column found.", fg="red"), err=True)
            sys.exit(1)

        data = table_hdu.data
        header = table_hdu.header.copy()
        n_total = len(data)
        click.echo(f"  {n_total} rows total")

        def _to_str(x):
            if isinstance(x, bytes):
                return x.decode().strip()
            return str(x).strip()

        mission_id_col = np.array([_to_str(x) for x in data["MISSION_ID"]])
        telescop_col   = np.array([_to_str(x) for x in data["TELESCOP"]])
        scan_col       = np.array([int(x) for x in data["SCAN"]]) if "SCAN" in data.dtype.names else None

        keep = np.ones(n_total, dtype=bool)
        n_dropped_total = 0

        for mission_id, params in mission_params.items():
            if not params:
                continue
            drop_cfg = params.get("drop") or {}
            if not drop_cfg:
                continue

            mission_mask = mission_id_col == mission_id

            # drop.telescope: remove all rows for these telescopes in this mission
            for tele in (drop_cfg.get("telescope") or []):
                affected = mission_mask & (telescop_col == tele)
                n = int(np.sum(affected))
                if n:
                    keep &= ~affected
                    n_dropped_total += n
                    click.echo(f"  {mission_id} / {tele}: dropped {n} rows (drop.telescope)")

            # drop.scans
            scans_cfg = drop_cfg.get("scans") or {}
            if scans_cfg and scan_col is not None:

                # drop.scans.complete: any telescope
                for scan_num in (scans_cfg.get("complete") or []):
                    affected = mission_mask & (scan_col == int(scan_num))
                    n = int(np.sum(affected))
                    if n:
                        keep &= ~affected
                        n_dropped_total += n
                        click.echo(f"  {mission_id} / scan {scan_num} (all telescopes): dropped {n} rows")

                # drop.scans.telescope.<TELE>: specific telescope
                for tele, scan_list in ((scans_cfg.get("telescope") or {}).items()):
                    for scan_num in (scan_list or []):
                        affected = mission_mask & (telescop_col == tele) & (scan_col == int(scan_num))
                        n = int(np.sum(affected))
                        if n:
                            keep &= ~affected
                            n_dropped_total += n
                            click.echo(f"  {mission_id} / {tele} / scan {scan_num}: dropped {n} rows")

        n_kept = int(np.sum(keep))
        click.echo(f"\nDropped {n_dropped_total} rows, keeping {n_kept} / {n_total}")

        filtered_data = data[keep]
        new_hdu = _fits.BinTableHDU(data=filtered_data, header=header)
        new_hdu.name = table_hdu.name
        _fits.HDUList([primary_hdu, new_hdu]).writeto(output, overwrite=True)

    click.echo(click.style(f"\nWritten to {output}", fg="green"))


@click.command()
@click.option('--config', type=click.Path(exists=True),
              help='Path to config.toml file.')
@click.option('--fits', type=click.Path(exists=True),
              help='Path to FITS file (overrides config).')
@click.option('--order', type=int, default=None,
              help='Baseline polynomial order (default: read from config or 1).')
@click.option('--window', type=(int, int), default=None,
              help='Channel range [start end] to exclude from baseline fitting. '
                   'Example: --window 100 110')
@click.option('--output', type=click.Path(), default=None,
              help='Output FITS file path (default: read from config or reduced_data.fits).')
def apply_baseline(config, fits, order, window, output):
    """
    Apply polynomial baseline subtraction to FITS spectra.
    
    Reads baseline order and window from config.toml [reduction] section,
    or accepts them via command-line options (CLI overrides config).
    
    The baseline subtraction:
    - Fits a polynomial to each spectrum
    - Ignores NaN channels
    - Optionally ignores a specified channel range (for strong emission lines)
    - Uses sigma-clipping to avoid bright features affecting the baseline
    
    Example:
        apply_baseline --config config.toml
        apply_baseline --config config.toml --order 2 --window 100 120
        apply_baseline --fits clean_data.fits --order 1 --output reduced.fits
    """
    try:
        # Load config
        config_path = config if config else None
        cfg = get_config(config_path) if config_path else {}
        
        # Determine FITS file to process
        if fits:
            # Use provided FITS file
            hdul = read_fits(fits)
            if output is None:
                output = cfg.get('output', {}).get('reduced_fits', 'reduced_data.fits')
        else:
            # Read from config
            hdul = None  # Let apply_baseline_from_config handle reading
            if output is None:
                output = cfg.get('output', {}).get('reduced_fits', 'reduced_data.fits')
        
        # Determine baseline order
        if order is None:
            reduction_cfg = cfg.get('reduction', {})
            order = reduction_cfg.get('baseline_order', reduction_cfg.get('baseline', 1))
            try:
                order = int(order)
            except Exception:
                order = 1
        
        # Use apply_baseline_from_config
        output_path = apply_baseline_from_config(
            config_path=config_path,
            hdul=hdul,
            output_path=output,
            overwrite=True,
            window=window
        )
        
        click.echo(f"\n✓ Baseline subtraction complete.")
        click.echo(f"  Output: {output_path}")
        click.echo(f"  Order: {order}")
        if window:
            click.echo(f"  Ignored window: channels {window[0]}-{window[1]}")
        click.echo()
        
    except FileNotFoundError as e:
        click.echo(click.style(f"Error: {e}", fg="red"), err=True)
        sys.exit(1)
    except ValueError as e:
        click.echo(click.style(f"Error: {e}", fg="red"), err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(click.style(f"Unexpected error: {e}", fg="red"), err=True)
        import traceback
        traceback.print_exc()
        sys.exit(1)


@click.command()
@click.option('--config', type=click.Path(exists=True), 
              help='Path to config.toml file.')
@click.option('--fits', type=click.Path(exists=True), 
              help='Path to FITS file (overrides config).')
@click.option('--clean', is_flag=True, default=False,
              help='Use clean_data.fits from config [output] instead of input file.')
@click.option('--unblank', is_flag=True, default=False,
              help='Fill NaN values using linear interpolation.')
@click.option('--fill-telluric-noise', 'fill_telluric_noise', is_flag=True, default=False,
              help="Fill each mission's telluric line window with Gaussian noise "
                   '(matched to that spectrum\'s own continuum std) BEFORE baselining. '
                   'Always runs before --baseline regardless of flag order on the command '
                   'line, since a real un-excluded telluric feature would otherwise bias '
                   'the baseline fit for missions whose telluric position falls outside a '
                   'single fixed --baseline-window. Per-mission telluric_line_center/width '
                   'are read from the mission parameters YAML (--mission-parameters, or '
                   'config [pca]/[input].mission_parameters, falling back to the bundled file).')
@click.option('--mission-parameters', 'mission_parameters', type=click.Path(exists=True), default=None,
              help='Path to mission parameters YAML file for --fill-telluric-noise '
                   '(overrides config [pca]/[input].mission_parameters).')
@click.option('--baseline', is_flag=True, default=False,
              help='Apply baseline subtraction.')
@click.option('--baseline-order', type=int, default=None,
              help='Polynomial order for baseline (default: read from config or 1).')
@click.option('--baseline-window', type=(int, int), default=None,
              help='Channel range [start end] to exclude from baseline fitting (default: read from config [reduction].line_window or None to use all channels).')
@click.option('--smooth', type=int, default=None, is_flag=False, flag_value=5,
              help='Apply boxcar smoothing. Optional window size (default: 5). E.g. --smooth or --smooth 7.')
@click.option('--decimate', is_flag=True, default=False,
              help='Decimate spectra by taking every Nth channel after smoothing, using the --smooth window size as the decimation factor.')
@click.option('--extract', is_flag=True, default=False,
              help='Extract velocity range from config [reduction].extract before baselining. Without this flag the spectrum length is left unchanged.')
@click.option('--output', type=click.Path(), default=None,
              help='Output FITS file path.')
def reduce_spectra_cmd(config, fits, clean, unblank, fill_telluric_noise, mission_parameters,
                       baseline, baseline_order, baseline_window,
                       smooth, decimate, extract, output):
    """
    Perform spectral reduction with selected methods.

    Applies reduction methods in sequence:
    1. Unblank (--unblank): Fill NaN values using interpolation
    2. Extract (--extract): Trim spectrum to velocity range from config [reduction].extract
    3. Fill telluric noise (--fill-telluric-noise): Overwrite each mission's telluric
       window with Gaussian noise, before baselining.
    4. Baseline subtraction (--baseline): Remove polynomial baseline
    5. Smoothing (--smooth [N]): Apply boxcar smoothing with window N (default: 5).
    6. Decimation (--decimate): Keep every Nth channel using the --smooth window size as N; VELOCITY_AXIS is updated accordingly.

    --extract, --smooth, and --decimate must be given explicitly. They are NOT applied
    automatically from config, so re-running reduce_spectra on an already-processed file
    (e.g. to re-baseline) will not change the spectrum length or channel count.

    Input file:
    - By default: Uses [input].fits_file from config
    - With --fits: Uses specified FITS file
    - With --clean: Uses [output].clean_fits from config

    If output path is not specified:
    - Uses config [output].reduced_fits if available
    - Falls back to 'reduced_data.fits' with a warning

    Examples:
        reduce_spectra --config config.toml --baseline
        reduce_spectra --config config.toml --fits clean_1.fits --baseline --output clean_2.fits
        reduce_spectra --config config.toml --extract --unblank --baseline --smooth --decimate
        reduce_spectra --config config.toml --unblank --baseline --baseline-order 2 --baseline-window 100 120
        reduce_spectra --fits clean_data.fits --unblank --baseline --smooth --output my_reduced.fits
    """
    try:
        # Load config
        config_path = config if config else None
        cfg = get_config(config_path) if config_path else {}
        
        # Determine FITS file to process
        if fits:
            hdul = read_fits(fits)
        elif clean:
            # Use clean_data.fits from config
            output_cfg = cfg.get('output', {})
            clean_fits_path = output_cfg.get('clean_fits', None)
            if not clean_fits_path:
                click.echo(click.style("Error: --clean flag specified but [output].clean_fits not defined in config", fg="red"), err=True)
                sys.exit(1)
            click.echo(f"Using clean FITS file: {clean_fits_path}")
            hdul = read_fits(clean_fits_path)
        else:
            hdul = None  # Let reduce_spectra_from_config handle reading
        
        # Build methods dictionary based on flags
        methods = {}
        
        if unblank:
            methods['unblank'] = {}

        reduction_cfg = cfg.get('reduction', {})

        # Fill telluric noise must run before baselining (see reduce_spectra() in
        # reduction/core.py). Resolve mission_parameters the same way prepare_for_pca does.
        fill_telluric_from_config = reduction_cfg.get('fill_telluric_noise', False)
        if fill_telluric_noise or fill_telluric_from_config:
            eff_mission_parameters = (mission_parameters
                                      or cfg.get('pca', {}).get('mission_parameters')
                                      or cfg.get('input', {}).get('mission_parameters'))
            methods['fill_telluric_noise'] = {'mission_parameters': eff_mission_parameters}

        # Extract velocity range only when --extract is explicitly requested.
        # Without the flag the spectrum length is left unchanged, so re-running
        # reduce_spectra on an already-processed file does not re-trim the channels.
        if extract:
            extract_cfg = reduction_cfg.get('extract', None)
            if extract_cfg is None:
                click.echo(click.style(
                    "Error: --extract specified but [reduction].extract is not defined in config",
                    fg='red'), err=True)
                sys.exit(1)
            try:
                if isinstance(extract_cfg, (list, tuple)) and len(extract_cfg) == 2:
                    extract_km_s = [float(extract_cfg[0]), float(extract_cfg[1])]
                    extract_m_s = [extract_km_s[0] * 1000.0, extract_km_s[1] * 1000.0]
                    methods['extract'] = extract_m_s
                    methods['extract_mode'] = 'velocity'
                    click.echo(f"Extracting velocity range [{extract_cfg[0]}, {extract_cfg[1]}] km/s")
            except (ValueError, TypeError, IndexError):
                click.echo(click.style("Warning: could not parse [reduction].extract from config", fg='yellow'), err=True)
        
        # Check if baseline should be applied (either from --baseline flag or from config)
        baseline_from_config = reduction_cfg.get('baseline', None)
        baseline_window_from_config = reduction_cfg.get('line_window', None)
        
        if baseline or baseline_from_config:
            # Get baseline order
            if baseline_order is None:
                if baseline:
                    # --baseline given explicitly without --baseline-order → default order 1
                    baseline_order = 1
                else:
                    # Baseline triggered only by config: use the config's baseline value
                    # as the order (falling back to 1 if it isn't a valid integer).
                    try:
                        baseline_order = int(baseline_from_config)
                    except (TypeError, ValueError):
                        baseline_order = 1
            
            # Get baseline window from line_window in config if not specified via CLI
            if baseline_window is None and baseline_window_from_config is not None:
                try:
                    if isinstance(baseline_window_from_config, (list, tuple)) and len(baseline_window_from_config) == 2:
                        # Window from config is in km/s
                        baseline_window = baseline_window_from_config
                except (ValueError, TypeError, IndexError):
                    baseline_window = None
            
            methods['baseline'] = {
                'order': baseline_order,
                'window': baseline_window
            }
            
            # If window is in km/s (from config), convert to m/s for internal use
            if baseline_window is not None and len(baseline_window) == 2:
                try:
                    # Check if values look like km/s (< 1000 typically means km/s)
                    if baseline_window[0] < 1000 and baseline_window[1] < 1000:
                        # Convert km/s to m/s
                        window_ms = [baseline_window[0] * 1000.0, baseline_window[1] * 1000.0]
                        methods['baseline']['window_m_s'] = window_ms
                        methods['baseline']['window_km_s'] = baseline_window
                except (ValueError, TypeError):
                    pass
        
        smooth_from_config = reduction_cfg.get('smooth', None)
        decimate_from_config = reduction_cfg.get('decimate', False)

        _smooth_config_val = int(smooth_from_config) if smooth_from_config is not None else None
        effective_smooth = smooth if smooth is not None else (_smooth_config_val if _smooth_config_val and _smooth_config_val > 1 else None)
        effective_decimate = decimate or bool(decimate_from_config)

        if effective_smooth is not None:
            methods['smooth'] = {
                'window_size': effective_smooth
            }

        if effective_decimate:
            if effective_smooth is None:
                click.echo(click.style("Error: --decimate requires --smooth to be set.", fg='red'), err=True)
                sys.exit(1)
            methods['decimate'] = {'factor': effective_smooth}

        if not methods:
            click.echo(click.style("No reduction methods selected. Use --unblank, --baseline, and/or --smooth [N].",
                                  fg="yellow"), err=True)
            sys.exit(1)
        
        # Apply reduction
        output_path = reduce_spectra_from_config(
            config_path=config_path,
            hdul=hdul,
            output_path=output,
            overwrite=True,
            methods=methods
        )
        
        click.echo(f"\n✓ Spectral reduction complete.")
        click.echo(f"  Output: {output_path}")
        click.echo(f"  Methods applied:")
        
        methods_applied = False
        if unblank:
            click.echo(f"    - Unblank (fill NaN values with linear interpolation)")
            methods_applied = True

        if 'fill_telluric_noise' in methods:
            source = " (from config)" if not fill_telluric_noise and fill_telluric_from_config else ""
            click.echo(f"    - Fill telluric noise (per-mission window, before baselining){source}")
            methods_applied = True

        # Check if baseline is in methods dict (either from --baseline flag or auto-applied from config)
        if 'baseline' in methods:
            baseline_info = methods['baseline']
            baseline_order = baseline_info.get('order', 1)
            baseline_window = baseline_info.get('window', None)
            source = " (from config)" if not baseline and baseline_from_config else ""
            if baseline_window:
                click.echo(f"    - Baseline subtraction (order={baseline_order}, window={baseline_window[0]}-{baseline_window[1]}){source}")
            else:
                click.echo(f"    - Baseline subtraction (order={baseline_order}){source}")
            methods_applied = True
        
        if effective_smooth is not None:
            source = " (from config)" if smooth is None and smooth_from_config is not None else ""
            click.echo(f"    - Smoothing (window={effective_smooth}){source}")
            methods_applied = True

        if effective_decimate:
            source = " (from config)" if not decimate and decimate_from_config else ""
            click.echo(f"    - Decimation (factor={effective_smooth}){source}")
            methods_applied = True

        if not methods_applied:
            click.echo(f"    (Extraction only, if configured; no other processing applied)")
        
        click.echo()
        
    except FileNotFoundError as e:
        click.echo(click.style(f"Error: {e}", fg="red"), err=True)
        sys.exit(1)
    except ValueError as e:
        click.echo(click.style(f"Error: {e}", fg="red"), err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(click.style(f"Unexpected error: {e}", fg="red"), err=True)
        import traceback
        traceback.print_exc()
        sys.exit(1)


@click.command()
@click.option('--config', type=click.Path(exists=True),
              help='Path to config.toml file.')
@click.option('--fits', 'fits_input', type=click.Path(exists=True), default=None,
              help='Input FITS file to average.')
@click.option('--reduced', is_flag=True, default=False,
              help='Use reduced_data.fits from config [output].reduced_fits.')
@click.option('--output', type=click.Path(), default=None,
              help='Output FITS file path for the averaged spectrum. If not given, no file is written.')
@click.option('--object', type=str, multiple=True, default=None,
              help='Object substring to include (e.g. "M51CENTER"). Can be repeated. '
                   'Defaults to [parameters].object from config if not specified.')
@click.option('--no-group', is_flag=True, default=False,
              help='Ignore OBJECT grouping and average all selected spectra into one.')
@click.option('--plot', type=click.Path(), default=None,
              help='Save plot to this path. If not specified, plot is shown interactively.')
def average_cmd(config, fits_input, reduced, output, object, no_group, plot):
    """
    Average spectra from a FITS file and optionally plot the result.

    Input priority: --fits > --reduced > config [input].fits_file.
    By default spectra are grouped by OBJECT; use --no-group for a single global average.
    Object filter defaults to [parameters].object from config if --object is not given.

    No output file is written unless --output is given.

    Examples:
        average --fits postpd.fits --plot avg.png
        average --fits postpd.fits --object M51CENTER --output avg.fits
        average --config config.toml --reduced --no-group --plot avg.png
    """
    import matplotlib.pyplot as plt
    try:
        from .reduction.core import average_spectra, average_spectra_from_config
        from .reduction.core import _extract_spectral_params, _create_velocity_axis
    except ImportError:
        from oi_zeigt.reduction.core import average_spectra, average_spectra_from_config
        from oi_zeigt.reduction.core import _extract_spectral_params, _create_velocity_axis

    try:
        config_path = config if config else None
        cfg = get_config(config_path) if config_path else {}

        # --- Load FITS ---
        if fits_input:
            hdul = read_fits(fits_input)
        elif reduced:
            reduced_fits_path = cfg.get('output', {}).get('reduced_fits')
            if not reduced_fits_path:
                click.echo(click.style("Error: --reduced specified but [output].reduced_fits not in config", fg='red'), err=True)
                sys.exit(1)
            hdul = read_fits(reduced_fits_path)
        else:
            hdul = read_fits_from_config(config_path)

        # --- Object filter ---
        object_filter = list(object) if object else []
        if not object_filter:
            obj_from_cfg = cfg.get('parameters', {}).get('object')
            if obj_from_cfg:
                object_filter = [obj_from_cfg]

        if object_filter:
            matrix_hdu = next(
                (hdu for hdu in hdul
                 if hasattr(hdu, 'data') and hdu.data is not None
                 and hdu.data.dtype.names is not None
                 and 'SPECTRUM' in hdu.data.dtype.names),
                None
            )
            if matrix_hdu is None:
                click.echo(click.style("Error: no SPECTRUM column found in FITS file", fg='red'), err=True)
                sys.exit(1)
            data_raw = matrix_hdu.data
            obj_col = data_raw['OBJECT'] if 'OBJECT' in data_raw.dtype.names else None
            if obj_col is None:
                click.echo(click.style("Warning: OBJECT column not found, ignoring --object filter", fg='yellow'), err=True)
            else:
                mask = np.zeros(len(obj_col), dtype=bool)
                for req in object_filter:
                    mask |= np.array([
                        req.upper() in (o.decode().strip() if isinstance(o, bytes) else str(o).strip()).upper()
                        for o in obj_col
                    ])
                if not np.any(mask):
                    click.echo(click.style(f"Error: no spectra matching {object_filter} found", fg='red'), err=True)
                    sys.exit(1)
                new_tbl = fits.BinTableHDU(data_raw[mask])
                new_tbl.name = matrix_hdu.name
                for k in matrix_hdu.header:
                    if k not in ('NAXIS1', 'NAXIS2', 'TFIELDS', ''):
                        try:
                            new_tbl.header[k] = matrix_hdu.header[k]
                        except (ValueError, KeyError):
                            pass
                hdul = fits.HDUList([hdul[0].copy(), new_tbl])
                click.echo(f"Object filter: {object_filter} → {int(mask.sum())} spectra")

        group_by_col = None if no_group else 'OBJECT'

        # --- Average ---
        avg_results = average_spectra(hdul, group_by=group_by_col)
        click.echo(f"\n✓ Averaging complete — grouped by: {'OBJECT' if group_by_col else 'none (global)'}")

        # --- Build velocity axis if possible ---
        try:
            sp = _extract_spectral_params(hdul)
            vel_ms = _create_velocity_axis(sp['velo_ref'], sp['deltav'], sp['crpix1_spec'], sp['nchans'])
            xaxis = vel_ms / 1000.0
            xlabel = 'Velocity (km/s)'
        except Exception:
            xaxis = None
            xlabel = 'Channel'

        # --- Plot ---
        fig, ax = plt.subplots(figsize=(12, 5))

        if isinstance(avg_results, dict) and 'avg_spectrum' in avg_results:
            # Single global average
            sp_arr = avg_results['avg_spectrum']
            std_arr = avg_results['std_spectrum']
            x = xaxis if xaxis is not None else np.arange(len(sp_arr))
            ax.plot(x, sp_arr, linewidth=1.5, color='steelblue', label=f"Average (N={avg_results['count']})")
            ax.legend(fontsize=9)
            click.echo(f"  N spectra averaged: {avg_results['count']}")
            click.echo(f"  Mean RMS: {avg_results['rms']:.4f} K")
        else:
            # Per-group averages
            for group_name, gdata in sorted(avg_results.items()):
                sp_arr = gdata['avg_spectrum']
                x = xaxis if xaxis is not None else np.arange(len(sp_arr))
                ax.plot(x, sp_arr, linewidth=1.2, alpha=0.85,
                        label=f"{group_name} (N={gdata['count']})")
                click.echo(f"  {group_name}: N={gdata['count']}, RMS={gdata['rms']:.4f} K")
            ax.legend(fontsize=8, loc='best')

        ax.set_xlabel(xlabel)
        ax.set_ylabel('T$_A^*$ (K)')
        ax.axhline(0, color='k', linewidth=0.5, linestyle='--')
        ax.grid(True, alpha=0.3)
        ax.set_title('Averaged spectrum' + (f' — {", ".join(object_filter)}' if object_filter else ''))
        fig.tight_layout()

        if plot:
            fig.savefig(plot, dpi=150)
            click.echo(f"  Plot saved: {plot}")
        else:
            plt.show()
        plt.close(fig)

        # --- Write output only if requested ---
        if output:
            average_spectra_from_config(
                config_path=config_path,
                hdul=hdul,
                output_fits=output,
                group_by=group_by_col,
                overwrite=True,
            )
            click.echo(f"  Output FITS: {output}")

    except FileNotFoundError as e:
        click.echo(click.style(f"Error: {e}", fg='red'), err=True)
        sys.exit(1)
    except ValueError as e:
        click.echo(click.style(f"Error: {e}", fg='red'), err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(click.style(f"Unexpected error: {e}", fg='red'), err=True)
        import traceback
        traceback.print_exc()
        sys.exit(1)


@click.command()
@click.option('--config', type=click.Path(exists=True), 
              help='Path to config.toml file.')
@click.option('--fits', type=click.Path(exists=True), 
              help='Path to FITS file (overrides config [input].fits_file).')
@click.option('--reduced', is_flag=True, default=False,
              help='Use reduced_data.fits instead of the input data.')
@click.option('--clean', is_flag=True, default=False,
              help='Use clean_data.fits from config [output].clean_fits.')
@click.option('--rejected', is_flag=True, default=False,
              help='Use rejected_data.fits from config [output].rejected_fits.')
@click.option('--prepared', is_flag=True, default=False,
              help='Use prepared_for_pca.fits from config [output].prepared_for_pca.')
@click.option('--pcad', is_flag=True, default=False,
              help='Use pcad.fits from config [output].pcad_fits.')
@click.option('--postfiltered', is_flag=True, default=False,
              help='Use post_filtered.fits from config [output].post_filtered_fits.')
@click.option('--postprocessed', is_flag=True, default=False,
              help='Use post_processed_fits from config [output].post_processed_fits.')
@click.option('--object', type=str, default=None,
              help='Filter by object name. If not specified, uses "object" from config.toml if available.')
@click.option('--all', 'all_metrics', is_flag=True, default=False,
              help='Include all available quality metrics.')
@click.option('--rmsratio', is_flag=True, default=False,
              help='Include RMSRATIO metric.')
@click.option('--rmsratiob', is_flag=True, default=False,
              help='Include RMSRATIOB (radiometer-based RMS ratio) metric.')
@click.option('--squality', is_flag=True, default=False,
              help='Include signal quality metric.')
@click.option('--roll-rms-n', is_flag=True, default=False,
              help='Include roll RMS metric.')
@click.option('--mh2o', is_flag=True, default=False,
              help='Include water vapor (MH2O) metric.')
@click.option('--tsys', is_flag=True, default=False,
              help='Include system temperature (Tsys) metric.')
@click.option('--tau-atm', is_flag=True, default=False,
              help='Include atmospheric optical depth (tau-atm) metric.')
@click.option('--chi-sqr', is_flag=True, default=False,
              help='Include chi-square metric.')
@click.option('--err-pwv', is_flag=True, default=False,
              help='Include PWV error metric.')
@click.option('--rms-baseline', is_flag=True, default=False,
              help='Include RMS baseline metric.')
@click.option('--bins', type=int, default=30,
              help='Number of histogram bins (default: 30).')
@click.option('--plot', type=click.Path(), default=None,
              help='Output path for plot. If not specified, plot is shown but not saved.')
def spechistogram_cmd(config, fits, reduced, clean, rejected, prepared, pcad, postfiltered, postprocessed, object, all_metrics, rmsratio, rmsratiob, squality, roll_rms_n, mh2o,
                      tsys, tau_atm, chi_sqr, err_pwv, rms_baseline, bins, plot):
    """
    Generate histograms of multiple spectral quality metrics.
    
    Creates a combined multi-panel figure showing histograms of selected quality metrics.
    Each metric shows mean and median lines.
    
    Input file priority: --fits option > output flags (--clean, --rejected, --prepared, --pcad, --postfiltered) > 
                        config [input].fits_file > --reduced flag > config [output].reduced_fits
    
    Available metrics:
    - --rmsratio: RMS ratio quality metric
    - --rmsratiob: RMS ratio vs radiometer equation (from post_process_data)
    - --squality: Signal quality flag
    - --roll-rms-n: RMS rolloff (mean across channels)
    - --mh2o: Water vapor column density
    - --tsys: System temperature
    - --tau-atm: Atmospheric optical depth
    - --chi-sqr: Chi-square fit value
    - --err-pwv: PWV error
    - --rms-baseline: RMS baseline metric
    
    Use --all to include all available metrics.
    
    Examples:
        spechistogram --config config.toml --rmsratio --tsys
        spechistogram --config config.toml --fits custom_data.fits --rmsratio
        spechistogram --config config.toml --clean --all
        spechistogram --config config.toml --prepared --rmsratio
        spechistogram --config config.toml --reduced --all
        spechistogram --config config.toml --all --object M51 --plot metrics.png
    """
    try:
        # Import fits module with alias to avoid collision with 'fits' parameter
        from astropy.io import fits as fits_module
        
        # Load config
        config_path = config if config else None
        cfg = get_config(config_path) if config_path else {}
        
        # Determine FITS file to process
        # Priority: --fits option > output flags (--clean, --rejected, --prepared, --pcad, --postfiltered) > 
        #           config [input].fits_file > --reduced flag > config [output].reduced_fits
        fits_file = None
        
        if fits:
            fits_file = fits
            file_description = f"custom FITS file"
        elif clean:
            output_cfg = cfg.get('output', {})
            fits_file = output_cfg.get('clean_fits')
            file_description = f"clean data"
            if not fits_file:
                click.echo(click.style("Error: --clean flag specified but [output].clean_fits not defined in config", fg="red"), err=True)
                sys.exit(1)
        elif rejected:
            output_cfg = cfg.get('output', {})
            fits_file = output_cfg.get('rejected_fits')
            file_description = f"rejected data"
            if not fits_file:
                click.echo(click.style("Error: --rejected flag specified but [output].rejected_fits not defined in config", fg="red"), err=True)
                sys.exit(1)
        elif prepared:
            output_cfg = cfg.get('output', {})
            fits_file = output_cfg.get('prepared_for_pca')
            file_description = f"prepared for PCA"
            if not fits_file:
                click.echo(click.style("Error: --prepared flag specified but [output].prepared_for_pca not defined in config", fg="red"), err=True)
                sys.exit(1)
        elif pcad:
            output_cfg = cfg.get('output', {})
            fits_file = output_cfg.get('pcad_fits')
            file_description = f"PCAD"
            if not fits_file:
                click.echo(click.style("Error: --pcad flag specified but [output].pcad_fits not defined in config", fg="red"), err=True)
                sys.exit(1)
        elif postfiltered:
            output_cfg = cfg.get('output', {})
            fits_file = output_cfg.get('post_filtered_fits')
            file_description = f"post-filtered"
            if not fits_file:
                click.echo(click.style("Error: --postfiltered flag specified but [output].post_filtered_fits not defined in config", fg="red"), err=True)
                sys.exit(1)
        elif postprocessed:
            output_cfg = cfg.get('output', {})
            fits_file = output_cfg.get('post_processed_fits')
            file_description = f"post-processed"
            if not fits_file:
                click.echo(click.style("Error: --postprocessed flag specified but [output].post_processed_fits not defined in config", fg="red"), err=True)
                sys.exit(1)
        elif reduced:
            output_cfg = cfg.get('output', {})
            fits_file = output_cfg.get('reduced_fits')
            file_description = f"reduced data"
            if not fits_file:
                click.echo(click.style("Error: --reduced flag specified but [output].reduced_fits not defined in config", fg="red"), err=True)
                sys.exit(1)
        else:
            # Try to get from config input.fits_file
            try:
                fits_from_config = cfg.get("input", {}).get("fits_file")
                if fits_from_config:
                    fits_file = fits_from_config
                    file_description = f"input data"
                else:
                    # Fallback to read_fits_from_config
                    click.echo("Reading FITS file from config (default input)...")
                    hdul = read_fits_from_config(config_path)
                    fits_file = None
            except (FileNotFoundError, KeyError):
                # Fallback to old behavior
                click.echo("Reading FITS file from config (default input)...")
                hdul = read_fits_from_config(config_path)
                fits_file = None
        
        # Open FITS file if we have a path
        if fits_file:
            click.echo(f"Analyzing {file_description}: {fits_file}")
            hdul = fits_module.open(fits_file)
        
        # Determine object filter
        if object is None:
            # Try to get from config
            parameters_cfg = cfg.get('parameters', {})
            obj_from_config = parameters_cfg.get('object', None)
            object_filter = obj_from_config
        else:
            object_filter = object
        
        # Build list of metrics to plot
        metrics = []
        
        # If --all is specified, enable all metrics
        if all_metrics:
            metrics = ['rmsratio', 'rmsratiob', 'squality', 'roll_rms_n', 'mh2o', 'tsys', 'tau_atm', 'chi_sqr', 'err_pwv', 'rms_baseline']
        else:
            # Build list from individual flags
            if rmsratio:
                metrics.append('rmsratio')
            if rmsratiob:
                metrics.append('rmsratiob')
            if squality:
                metrics.append('squality')
            if roll_rms_n:
                metrics.append('roll_rms_n')
            if mh2o:
                metrics.append('mh2o')
            if tsys:
                metrics.append('tsys')
            if tau_atm:
                metrics.append('tau_atm')
            if chi_sqr:
                metrics.append('chi_sqr')
            if err_pwv:
                metrics.append('err_pwv')
            if rms_baseline:
                metrics.append('rms_baseline')
        
        if not metrics:
            click.echo(click.style("Error: No metrics selected. Use --all or at least one of: --rmsratio, --rmsratiob, --squality, --roll-rms-n, --mh2o, --tsys, --tau-atm, --chi-sqr, --err-pwv, --rms-baseline", fg="red"), err=True)
            sys.exit(1)
        
        # Generate histograms
        fig, stats_dict = get_spechistogram(
            hdul,
            metrics=metrics,
            object_filter=object_filter,
            bins=bins
        )
        
        # Display statistics
        click.echo("\n" + "="*90)
        click.echo("SPECTRAL QUALITY METRICS STATISTICS")
        if object_filter:
            click.echo(f"Object filter: {object_filter}")
        click.echo("="*90)
        
        for metric, stats in stats_dict.items():
            click.echo(f"\n{metric.upper()}:")
            click.echo(f"  Count:        {stats['count']}")
            click.echo(f"  Mean:         {stats['mean']:.6f}")
            click.echo(f"  Median:       {stats['median']:.6f}")
            click.echo(f"  Std Dev:      {stats['std']:.6f}")
            click.echo(f"  Range:        [{stats['min']:.6f}, {stats['max']:.6f}]")
            click.echo(f"  5-95 %ile:    [{stats['percentile_5']:.6f}, {stats['percentile_95']:.6f}]")
        
        click.echo("\n" + "="*90 + "\n")
        
        # Save plot if requested
        if plot:
            plt.savefig(plot, dpi=150)
            click.echo(f"✓ Plot saved: {plot}\n")
        
        # Always show the plot
        plt.show()
        plt.close()
        
        hdul.close()
        
    except FileNotFoundError as e:
        click.echo(click.style(f"Error: {e}", fg="red"), err=True)
        sys.exit(1)
    except ValueError as e:
        click.echo(click.style(f"Error: {e}", fg="red"), err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(click.style(f"Unexpected error: {e}", fg="red"), err=True)
        import traceback
        traceback.print_exc()
        sys.exit(1)


@click.command()
@click.option('--config', type=click.Path(exists=True),
              help='Path to config.toml file.')
@click.option('--fits', 'fits_file', type=click.Path(exists=True), default=None,
              help='Path to FITS file to analyse directly. Overrides all dataset flags.')
@click.option('--reduced', is_flag=True, default=False,
              help='Use reduced_data.fits instead of the input data.')
@click.option('--clean', is_flag=True, default=False,
              help='Use clean_fits from config.toml instead of the input data.')
@click.option('--rejected', is_flag=True, default=False,
              help='Use rejected_fits from config.toml instead of the input data.')
@click.option('--prepared', is_flag=True, default=False,
              help='Use prepared_for_pca from config.toml instead of the input data.')
@click.option('--pcad', is_flag=True, default=False,
              help='Use pcad_fits (PCA-corrected) from config.toml instead of the input data.')
@click.option('--object', type=str, default=None,
              help='Filter by object name (e.g., "M51"). Defaults to "object" from config.toml if not specified. If not in config either, all spectra are analyzed.')
@click.option('--bins', type=int, default=30,
              help='Number of histogram bins (default: 30).')
@click.option('--plot', type=click.Path(), default=None,
              help='Output path for plot (e.g., rmsratio_hist.png). If not specified, plot is shown but not saved.')
def rmsratio_cmd(config, fits_file, reduced, clean, rejected, prepared, pcad, object, bins, plot):
    """
    Analyze RMSRATIO quality metric and generate histogram.
    
    Reads input file from (in order of priority):
    1. --fits if specified
    2. [output].reduced_fits if --reduced flag is specified
    3. [output].clean_fits if --clean flag is specified
    4. [output].rejected_fits if --rejected flag is specified
    5. [output].prepared_for_pca if --prepared flag is specified
    6. [output].pcad_fits if --pcad flag is specified
    7. [input].fits_file from config.toml otherwise
    
    Computes statistics and creates a histogram of RMSRATIO values
    for all spectra or a specific object.
    
    Object name defaults to [parameters].object from config.toml if not specified.
    
    RMSRATIO is a quality metric where lower values indicate better quality.
    
    Outputs:
    - Console: Statistics (mean, median, std, min, max, percentiles)
    - Plot: Histogram with mean and median lines (always displayed, optionally saved)
    
    Examples:
        rmsratio --config config.toml
        rmsratio --config config.toml --reduced
        rmsratio --config config.toml --clean
        rmsratio --config config.toml --rejected
        rmsratio --config config.toml --prepared
        rmsratio --config config.toml --pcad
        rmsratio --config config.toml --object M51
        rmsratio --config config.toml --bins 50 --plot hist.png
    """
    try:
        # Load config
        config_path = config if config else None
        cfg = get_config(config_path) if config_path else {}
        
        # Determine FITS file to process
        fits_path = None
        if fits_file:
            hdul = read_fits(fits_file)
        elif reduced or clean or rejected or prepared or pcad:
            output_cfg = cfg.get('output', {})
            if reduced:
                fits_path = output_cfg.get('reduced_fits', None)
                flag_name = "reduced_fits"
            elif clean:
                fits_path = output_cfg.get('clean_fits', None)
                flag_name = "clean_fits"
            elif rejected:
                fits_path = output_cfg.get('rejected_fits', None)
                flag_name = "rejected_fits"
            elif prepared:
                fits_path = output_cfg.get('prepared_for_pca', None)
                flag_name = "prepared_for_pca"
            elif pcad:
                fits_path = output_cfg.get('pcad_fits', None)
                flag_name = "pcad_fits"
            
            if not fits_path:
                click.echo(click.style(f"Error: --{flag_name.replace('_', '-')} flag specified but [output].{flag_name} not defined in config", fg="red"), err=True)
                sys.exit(1)
            hdul = fits.open(fits_path)
        else:
            # Use input file from config
            hdul = read_fits_from_config(config_path)
        
        # Determine object filter
        if object is None:
            # Try to get from config
            parameters_cfg = cfg.get('parameters', {})
            object = parameters_cfg.get('object', None)
        
        # Compute statistics
        stats = rmsratio_statistics(hdul, object_filter=object)
        
        # Display statistics
        click.echo("\n" + "="*70)
        click.echo("RMSRATIO STATISTICS")
        if object:
            click.echo(f"Object filter: {object}")
        click.echo("="*70)
        click.echo(f"Count:        {stats['count']}")
        click.echo(f"Mean:         {stats['mean']:.6f}")
        click.echo(f"Median:       {stats['median']:.6f}")
        click.echo(f"Std Dev:      {stats['std']:.6f}")
        click.echo(f"Min:          {stats['min']:.6f}")
        click.echo(f"Max:          {stats['max']:.6f}")
        click.echo(f"5th percentile: {stats['percentile_5']:.6f}")
        click.echo(f"95th percentile: {stats['percentile_95']:.6f}")
        click.echo("="*70 + "\n")
        
        # Generate histogram
        fig, hist, bin_edges = get_rmsratio_histogram(
            hdul,
            object_filter=object,
            bins=bins,
            figsize=(12, 6)
        )
        
        # Save plot if requested
        if plot:
            plt.savefig(plot, dpi=150)
            click.echo(f"✓ Plot saved: {plot}\n")
        
        # Always show the plot
        plt.show()
        plt.close()
        
        hdul.close()
        
    except FileNotFoundError as e:
        click.echo(click.style(f"Error: {e}", fg="red"), err=True)
        sys.exit(1)
    except ValueError as e:
        click.echo(click.style(f"Error: {e}", fg="red"), err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(click.style(f"Unexpected error: {e}", fg="red"), err=True)
        import traceback
        traceback.print_exc()
        sys.exit(1)


@click.command()
@click.option('--config', type=click.Path(exists=True), 
              help='Path to config.toml file.')
@click.option('--reduced', is_flag=True, default=False,
              help='Use reduced_data.fits instead of the input data.')
@click.option('--clean', is_flag=True, default=False,
              help='Use clean_fits instead of the input data.')
@click.option('--object', type=str, default=None,
              help='Filter by object name (partial match). If not specified, uses "object" from config.toml if available.')
@click.option('--column', type=str, default=None,
              help='Column to map (e.g., "TSYS", "TAU-ATM", "RMSRATIO"). Required.')
@click.option('--beamsize', type=float, default=None,
              help='Beam size in degrees for gridding kernel. If not specified, reads from config [gridding].beamsize_arcsec.')
@click.option('--pixsize', type=float, default=None,
              help='Map pixel size in degrees. If not specified, uses beamsize/3 or config [gridding].pixel_size_arcsec.')
@click.option('--scatter', is_flag=True, default=False,
              help='Show observation points as scatter plot on map.')
@click.option('--plot', type=click.Path(), default=None,
              help='Output path for plot (e.g., map.png). If not specified, plot is shown but not saved.')
@click.option('--fits-output', type=click.Path(), default=None,
              help='Output path for FITS file (e.g., map.fits). If not specified, FITS is not saved.')
def map_column_cmd(config, reduced, clean, object, column, beamsize, pixsize, scatter, plot, fits_output):
    """
    Create a spatial map from a quality metric or observational parameter.
    
    Maps a single column (e.g., TSYS, TAU-ATM, RMSRATIO) spatially using proper WCS 
    coordinates with cygrid-based gridding for optimal Gaussian kernel weighting.
    
    Reads input file from config.toml [input].fits_file by default, or uses
    [output].reduced_fits if --reduced flag is specified.

    Gridding parameters read from config.toml [gridding] section if not specified on command line.

    Examples:
        map_column --config config.toml --column TSYS
        map_column --config config.toml --reduced --column TAU-ATM --object M51
        map_column --config config.toml --column RMSRATIO --beamsize 0.3 --plot map.png
        map_column --config config.toml --column MH2O --pixsize 0.05 --object M51
    """
    try:
        from oi_zeigt.mapping.gridding import get_gridding_params_from_config
        
        if not column:
            click.echo(click.style("Error: --column is required", fg="red"), err=True)
            sys.exit(1)

        from .mapping.gridding import create_map_from_column

        # Load config
        config_path = config if config else None
        cfg = get_config(config_path) if config_path else {}

        # Get beamsize and pixsize from config if not provided on CLI
        beamsize_deg, pixsize_deg = get_gridding_params_from_config(
            config_path=config_path,
            beamsize_deg=beamsize,
            pixsize_deg=pixsize
        )

        # Determine FITS file to process
        if reduced:
            # Use reduced_data.fits from config
            output_cfg = cfg.get('output', {})
            reduced_fits_path = output_cfg.get('reduced_fits', None)
            if not reduced_fits_path:
                click.echo(click.style("Error: --reduced flag specified but [output].reduced_fits not defined in config", fg="red"), err=True)
                sys.exit(1)
            click.echo(f"Reading FITS file: {reduced_fits_path}")
            hdul = fits.open(reduced_fits_path)
        elif clean:
            # Use clean_fits from config
            output_cfg = cfg.get('output', {})
            clean_fits_path = output_cfg.get('clean_fits', None)
            if not clean_fits_path:
                click.echo(click.style("Error: --clean flag specified but [output].clean_fits not defined in config", fg="red"), err=True)
                sys.exit(1)
            click.echo(f"Reading FITS file: {clean_fits_path}")
            hdul = fits.open(clean_fits_path)
        else:
            # Use input file from config
            hdul = read_fits_from_config(config_path)

        # Determine object filter
        if object is None:
            # Try to get from config
            parameters_cfg = cfg.get('parameters', {})
            obj_from_config = parameters_cfg.get('object', None)
            object_filter = obj_from_config
        else:
            object_filter = object

        # Create map
        grid_map, wcs_header, fig = create_map_from_column(
            hdul,
            column_name=column.upper(),
            object_filter=object_filter,
            beamsize_deg=beamsize_deg,
            pixsize=pixsize_deg,
            show_scatter=scatter
        )

        click.echo(f"\n✓ Map created from column: {column.upper()}")
        click.echo(f"  Beam size: {beamsize_deg:.4f}°")
        if pixsize_deg:
            click.echo(f"  Pixel size: {pixsize_deg:.6f}°")
        if object_filter:
            click.echo(f"  Object filter: {object_filter}")
        click.echo(f"  Map dimensions: {grid_map.shape[1]} × {grid_map.shape[0]}")
        
        # Save plot if requested
        if plot:
            plt.savefig(plot, dpi=150)
            click.echo(f"  Plot saved: {plot}")
        
        # Save FITS file if requested
        if fits_output:
            from oi_zeigt.mapping.gridding import save_map_to_fits
            save_map_to_fits(grid_map, wcs_header, fits_output, 
                           beam_maj_deg=beamsize, overwrite=True)
            click.echo(f"  FITS saved: {fits_output}")
        
        click.echo()
        
        # Always show the plot
        plt.show()
        plt.close()

        hdul.close()

    except FileNotFoundError as e:
        click.echo(click.style(f"Error: {e}", fg="red"), err=True)
        sys.exit(1)
    except ValueError as e:
        click.echo(click.style(f"Error: {e}", fg="red"), err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(click.style(f"Unexpected error: {e}", fg="red"), err=True)
        import traceback
        traceback.print_exc()
        sys.exit(1)


@click.command()
@click.option('--config', type=click.Path(exists=True), 
              help='Path to config.toml file.')
@click.option('--fits', type=click.Path(exists=True), default=None,
              help='Path to FITS file to process. Overrides [input].fits_file from config.')
@click.option('--reduced', is_flag=True, default=False,
              help='Use reduced_data.fits instead of the input data.')
@click.option('--clean', is_flag=True, default=False,
              help='Use clean_fits instead of the input data.')
@click.option('--pcad', is_flag=True, default=False,
              help='Use pca_corrected.fits (output.pcad_fits from config).')
@click.option('--rejected', is_flag=True, default=False,
              help='Use rejected_data.fits (output.rejected_fits from config).')
@click.option('--postfiltered', is_flag=True, default=False,
              help='Use post_filtered_fits (output.post_filtered_fits from config).')
@click.option('--prepared', is_flag=True, default=False,
              help='Use prepared_for_pca.fits (output.prepared_for_pca from config).')
@click.option('--object', type=str, default=None,
              help='Filter by object name (partial match). If not specified, uses "object" from config.toml if available.')
@click.option('--beamsize', type=float, default=None,
              help='Beam size in degrees for gridding kernel. If not specified, reads from config [gridding].beamsize_arcsec.')
@click.option('--pixsize', type=float, default=None,
              help='Map pixel size in degrees. If not specified, uses beamsize/3 or config [gridding].pixel_size_arcsec.')
@click.option('--velocity-range', type=(float, float), default=None, nargs=2,
              help='Velocity range in km/s (e.g., --velocity-range 450 500). If not specified, integrates entire spectrum.')
@click.option('--weight-column', type=str, default=None,
              help='Column name for per-spectrum weighting (e.g., RMSRATIO). Overrides --weight-spectra.')
@click.option('--weight-spectra', is_flag=True, default=False,
              help='Weight each spectrum by its RMSRATIOB value (falls back to RMSRATIO if absent). '
                   'Shorthand for --weight-column RMSRATIOB.')
@click.option('--weight-channels', is_flag=True, default=False,
              help='Weight each channel by exp(-tau)/T_sys using the TSYS_INDEX/TAU_SIG_INDEX columns '
                   '(added by prepare_for_pca). Falls back to uniform if columns are absent.')
@click.option('--scatter', is_flag=True, default=False,
              help='Show observation points as scatter plot on map.')
@click.option('--plot', type=click.Path(), default=None,
              help='Output path for plot (e.g., integrated_map.png). If not specified, plot is shown but not saved.')
@click.option('--fits-output', type=click.Path(), default=None,
              help='Output path for FITS file (e.g., integrated_map.fits). If not specified, FITS is not saved.')
def map_integrated_cmd(config, fits, reduced, clean, pcad, rejected, postfiltered, prepared, object, beamsize, pixsize, velocity_range, weight_column, weight_spectra, weight_channels, scatter, plot, fits_output):
    """
    Create a spatial map of integrated spectral intensity with optional per-spectrum weighting.
    
    Integrates the spectrum across all frequency channels (or specified velocity range)
    for each observation, then creates a WCS-based spatial map with cygrid gridding.
    Optionally weights each integrated spectrum by a quality metric (e.g., RMSRATIO).
    
    Reads input file from (in order of priority):
    1. --fits parameter if specified (overrides config)
    2. --pcad flag (uses output.pcad_fits from config)
    3. --rejected flag (uses output.rejected_fits from config)
    4. --postfiltered flag (uses output.post_filtered_fits from config)
    5. --prepared flag (uses output.prepared_for_pca from config)
    6. --reduced flag (uses output.reduced_fits from config)
    7. --clean flag (uses output.clean_fits from config)
    8. [input].fits_file from config.toml otherwise
    
    Gridding parameters read from config.toml [gridding] section if not specified on command line.

    Examples:
        map_integrated --config config.toml
        map_integrated --config config.toml --fits /path/to/data.fits
        map_integrated --config config.toml --pcad --object M51
        map_integrated --config config.toml --reduced --object M51
        map_integrated --config config.toml --prepared --beamsize 0.3 --plot integrated.png
        map_integrated --config config.toml --rejected --velocity-range 450 500
        map_integrated --config config.toml --postfiltered --velocity-range -50 50 --plot map_narrow_range.png
        map_integrated --config config.toml --pcad --object M51 --weight-column RMSRATIO --plot weighted_map.png
    """
    try:
        from oi_zeigt.mapping.gridding import get_gridding_params_from_config, create_integrated_map

        # Load config
        config_path = config if config else None
        cfg = get_config(config_path) if config_path else {}

        # Get beamsize and pixsize from config if not provided on CLI
        beamsize_deg, pixsize_deg = get_gridding_params_from_config(
            config_path=config_path,
            beamsize_deg=beamsize,
            pixsize_deg=pixsize
        )

        # Determine FITS file to process
        if fits:
            # Use directly specified FITS file (overrides everything)
            click.echo(f"Reading FITS file: {fits}")
            hdul = read_fits(fits)
        elif pcad:
            # Use pca_corrected.fits from config
            output_cfg = cfg.get('output', {})
            pcad_fits_path = output_cfg.get('pcad_fits', None)
            if not pcad_fits_path:
                click.echo(click.style("Error: --pcad flag specified but [output].pcad_fits not defined in config", fg="red"), err=True)
                sys.exit(1)
            click.echo(f"Reading FITS file: {pcad_fits_path}")
            hdul = read_fits(pcad_fits_path)
        elif rejected:
            # Use rejected_data.fits from config
            output_cfg = cfg.get('output', {})
            rejected_fits_path = output_cfg.get('rejected_fits', None)
            if not rejected_fits_path:
                click.echo(click.style("Error: --rejected flag specified but [output].rejected_fits not defined in config", fg="red"), err=True)
                sys.exit(1)
            click.echo(f"Reading FITS file: {rejected_fits_path}")
            hdul = read_fits(rejected_fits_path)
        elif postfiltered:
            # Use post_filtered_fits from config
            output_cfg = cfg.get('output', {})
            postfiltered_fits_path = output_cfg.get('post_filtered_fits', None)
            if not postfiltered_fits_path:
                click.echo(click.style("Error: --postfiltered flag specified but [output].post_filtered_fits not defined in config", fg="red"), err=True)
                sys.exit(1)
            click.echo(f"Reading FITS file: {postfiltered_fits_path}")
            hdul = read_fits(postfiltered_fits_path)
        elif prepared:
            # Use prepared_for_pca.fits from config
            output_cfg = cfg.get('output', {})
            prepared_fits_path = output_cfg.get('prepared_for_pca', None)
            if not prepared_fits_path:
                click.echo(click.style("Error: --prepared flag specified but [output].prepared_for_pca not defined in config", fg="red"), err=True)
                sys.exit(1)
            click.echo(f"Reading FITS file: {prepared_fits_path}")
            hdul = read_fits(prepared_fits_path)
        elif reduced:
            # Use reduced_data.fits from config
            output_cfg = cfg.get('output', {})
            reduced_fits_path = output_cfg.get('reduced_fits', None)
            if not reduced_fits_path:
                click.echo(click.style("Error: --reduced flag specified but [output].reduced_fits not defined in config", fg="red"), err=True)
                sys.exit(1)
            click.echo(f"Reading FITS file: {reduced_fits_path}")
            hdul = read_fits(reduced_fits_path)
        elif clean:
            # Use clean_fits from config
            output_cfg = cfg.get('output', {})
            clean_fits_path = output_cfg.get('clean_fits', None)
            if not clean_fits_path:
                click.echo(click.style("Error: --clean flag specified but [output].clean_fits not defined in config", fg="red"), err=True)
                sys.exit(1)
            click.echo(f"Reading FITS file: {clean_fits_path}")
            hdul = read_fits(clean_fits_path)
        else:
            # Use input file from config
            hdul = read_fits_from_config(config_path)

        # Determine object filter
        if object is None:
            # Try to get from config
            parameters_cfg = cfg.get('parameters', {})
            obj_from_config = parameters_cfg.get('object', None)
            object_filter = obj_from_config
        else:
            object_filter = object

        # Resolve per-spectrum weight column (map_integrated_cmd block)
        eff_weight_column = weight_column
        if eff_weight_column is None and weight_spectra:
            eff_weight_column = 'RMSRATIOB'  # fallback to RMSRATIO handled inside create_integrated_map

        # Create integrated map
        grid_map, wcs_header, fig = create_integrated_map(
            hdul,
            object_filter=object_filter,
            beamsize_deg=beamsize_deg,
            pixsize=pixsize_deg,
            show_scatter=scatter,
            velocity_range=velocity_range,
            weight_column=eff_weight_column,
            channel_weights=weight_channels,
        )

        click.echo(f"\n✓ Integrated intensity map created")
        click.echo(f"  Beam size: {beamsize_deg*3600:.2f}\"")
        if pixsize_deg:
            click.echo(f"  Pixel size: {pixsize_deg*3600:.3f}\"")
        if object_filter:
            click.echo(f"  Object filter: {object_filter}")
        if eff_weight_column:
            click.echo(f"  Per-spectrum weighting: {eff_weight_column}")
        if weight_channels:
            click.echo(f"  Per-channel weighting: exp(-tau)/T_sys from TSYS/TAU_SIG calibration spectra")
        if velocity_range:
            click.echo(f"  Velocity range: {velocity_range[0]:.1f} - {velocity_range[1]:.1f} km/s")
        click.echo(f"  Map dimensions: {grid_map.shape[1]} × {grid_map.shape[0]}")

        # Save plot if requested
        if plot:
            plt.savefig(plot, dpi=150)
            click.echo(f"  Plot saved: {plot}")

        # Save FITS file if requested
        if fits_output:
            from oi_zeigt.mapping.gridding import save_map_to_fits
            save_map_to_fits(grid_map, wcs_header, fits_output, 
                           beam_maj_deg=beamsize_deg, overwrite=True)
            click.echo(f"  FITS saved: {fits_output}")

        click.echo()

        # Always show the plot
        plt.show()
        plt.close()

        hdul.close()

    except FileNotFoundError as e:
        click.echo(click.style(f"Error: {e}", fg="red"), err=True)
        sys.exit(1)
    except ValueError as e:
        click.echo(click.style(f"Error: {e}", fg="red"), err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(click.style(f"Unexpected error: {e}", fg="red"), err=True)
        import traceback
        traceback.print_exc()
        sys.exit(1)


@click.command()
@click.option('--config', type=click.Path(exists=True),
              help='Path to config.toml file.')
@click.option('--fits', 'fits_files', type=click.Path(exists=True), multiple=True,
              help='Explicit FITS files to include (repeat for each). Labels taken from basenames.')
@click.option('--clean',    is_flag=True, default=False, help='Include clean dataset from config (output.clean_fits).')
@click.option('--reduced',  is_flag=True, default=False, help='Include reduced dataset from config (output.reduced_fits).')
@click.option('--prepared', is_flag=True, default=False, help='Include prepared-for-PCA dataset from config (output.prepared_for_pca).')
@click.option('--pcad',     is_flag=True, default=False, help='Include PCA-corrected dataset from config (output.pcad_fits).')
@click.option('--object', type=str, default=None,
              help='Filter by object name. If not specified, uses "object" from config or auto-excludes calibration rows (TSYS, TAU_SIG, etc.).')
@click.option('--beamsize', type=float, default=None,
              help='Beam size in degrees for gridding kernel. If not specified, reads from config [gridding].beamsize_arcsec (default 14.1″).')
@click.option('--pixsize', type=float, default=None,
              help='Map pixel size in degrees. If not specified, uses beamsize/3 or config [gridding].pixel_size_arcsec.')
@click.option('--velocity-range', type=(float, float), default=None, nargs=2,
              help='Velocity range in km/s (e.g., --velocity-range 450 500). If not specified, integrates entire spectrum.')
@click.option('--weight-column', type=str, default=None,
              help='Column name for per-spectrum weighting (e.g., RMSRATIO). Overrides --weight-spectra.')
@click.option('--weight-spectra', is_flag=True, default=False,
              help='Weight each spectrum by its RMSRATIOB value (falls back to RMSRATIO if absent).')
@click.option('--weight-channels', is_flag=True, default=False,
              help='Weight each channel by exp(-tau)/T_sys using TSYS_INDEX/TAU_SIG_INDEX columns. '
                   'Falls back to uniform if columns are absent.')
@click.option('--scatter', is_flag=True, default=False,
              help='Show observation points as scatter plot on map.')
@click.option('--plot', type=click.Path(), default=None,
              help='Output path for plot (e.g., compare_map.png). If not specified, plot is shown but not saved.')
def compare_map_integrated_cmd(config, fits_files, clean, reduced, prepared, pcad, object, beamsize, pixsize, velocity_range, weight_column, weight_spectra, weight_channels, scatter, plot):
    """
    Compare integrated intensity maps side by side (max 3 per row).

    Datasets can come from explicit --fits files and/or config-defined datasets
    selected with --clean / --reduced / --prepared / --pcad flags. Both can be
    combined freely. When no flags and no --fits are given, defaults to showing
    reduced + prepared + PCA-corrected from config.

    Examples:
        compare_map_integrated --config config.toml
        compare_map_integrated --config config.toml --reduced --pcad
        compare_map_integrated --config config.toml --fits postpcarmsr_1.3.fits --fits postpcarmsr_1.5.fits --fits postpcarmsr_2.0.fits --object M51CENTER
        compare_map_integrated --config config.toml --pcad --fits postpcarmsr_1.3.fits
    """
    try:
        from oi_zeigt.mapping.gridding import (
            get_gridding_params_from_config, create_integrated_map, format_ra_hms, format_dec_dms, HAS_CYGRID
        )
        from matplotlib.ticker import FuncFormatter, MaxNLocator
        from pathlib import Path
        import math

        config_path = config if config else None
        cfg = get_config(config_path) if config_path else {}

        beamsize_deg, pixsize_deg = get_gridding_params_from_config(
            config_path=config_path,
            beamsize_deg=beamsize,
            pixsize_deg=pixsize,
        )

        if object is None:
            object_filter = cfg.get('parameters', {}).get('object', None)
        else:
            object_filter = object

        # --- Build dataset list ---
        # Start with config-flag datasets (in a logical order), then append explicit --fits
        datasets = []
        output_cfg = cfg.get('output', {})
        # Default (no flags, no --fits): show reduced + prepared + pcad from config
        use_defaults = not fits_files and not any([clean, reduced, prepared, pcad])

        if clean:
            p = output_cfg.get('clean_fits')
            if not p:
                click.echo(click.style("Error: output.clean_fits not in config", fg='red'), err=True); sys.exit(1)
            datasets.append(('Clean', p))
        if reduced or use_defaults:
            p = output_cfg.get('reduced_fits')
            if reduced and not p:
                click.echo(click.style("Error: output.reduced_fits not in config", fg='red'), err=True); sys.exit(1)
            if p: datasets.append(('Reduced', p))
        if prepared or use_defaults:
            p = output_cfg.get('prepared_for_pca')
            if prepared and not p:
                click.echo(click.style("Error: output.prepared_for_pca not in config", fg='red'), err=True); sys.exit(1)
            if p: datasets.append(('Prepared', p))
        if pcad or use_defaults:
            p = output_cfg.get('pcad_fits')
            if pcad and not p:
                click.echo(click.style("Error: output.pcad_fits not in config", fg='red'), err=True); sys.exit(1)
            if p: datasets.append(('PCA corrected', p))
        for fp in fits_files:
            datasets.append((Path(fp).stem, fp))

        if not datasets:
            click.echo(click.style("Error: no datasets to compare. Use --fits, --reduced, --prepared, --pcad, or --clean.", fg='red'), err=True)
            sys.exit(1)

        # Resolve per-spectrum weight column
        eff_weight_column = weight_column
        if eff_weight_column is None and weight_spectra:
            eff_weight_column = 'RMSRATIOB'

        # --- Grid each dataset ---
        maps = []
        for label, path in datasets:
            click.echo(f"Reading {label}: {path}")
            hdul = read_fits(path)
            grid_map, wcs_header, _fig = create_integrated_map(
                hdul,
                object_filter=object_filter,
                beamsize_deg=beamsize_deg,
                pixsize=pixsize_deg,
                show_scatter=False,
                velocity_range=velocity_range,
                weight_column=eff_weight_column,
                channel_weights=weight_channels,
            )
            plt.close(_fig)
            hdul.close()
            maps.append((label, grid_map, wcs_header))

        # --- Build figure: max 3 panels per row ---
        n = len(maps)
        ncols = min(n, 3)
        nrows = math.ceil(n / 3)
        fig, axes = plt.subplots(nrows, ncols, figsize=(6 * ncols, 6 * nrows),
                                 squeeze=False)
        # Flatten axes and hide any unused panels
        axes_flat = [axes[r][c] for r in range(nrows) for c in range(ncols)]
        for ax in axes_flat[n:]:
            ax.set_visible(False)

        for ax, (label, grid_map, wcs_header) in zip(axes_flat, maps):
            crval1 = wcs_header['CRVAL1']
            crval2 = wcs_header['CRVAL2']
            cdelt1 = wcs_header['CDELT1']
            cdelt2 = wcs_header['CDELT2']
            crpix1 = wcs_header['CRPIX1']
            crpix2 = wcs_header['CRPIX2']
            naxis1 = wcs_header['NAXIS1']
            naxis2 = wcs_header['NAXIS2']

            ra_max = crval1 + (1 - crpix1) * cdelt1
            ra_min = crval1 + (naxis1 - crpix1) * cdelt1
            dec_min = crval2 + (1 - crpix2) * cdelt2
            dec_max = crval2 + (naxis2 - crpix2) * cdelt2

            grid_map_display = np.fliplr(grid_map)
            im = ax.imshow(grid_map_display, origin='lower',
                           extent=[ra_min, ra_max, dec_min, dec_max],
                           cmap='viridis', aspect='auto')

            ra_formatter = FuncFormatter(lambda x, p: format_ra_hms(x))
            dec_formatter = FuncFormatter(lambda x, p: format_dec_dms(x))
            ax.xaxis.set_major_formatter(ra_formatter)
            ax.yaxis.set_major_formatter(dec_formatter)
            ax.xaxis.set_major_locator(MaxNLocator(nbins=4, integer=False))
            ax.invert_xaxis()

            ax.set_xlabel('RA (h:m:s)', fontsize=10)
            ax.set_ylabel('Dec (°:′:″)', fontsize=10)
            ax.set_title(label, fontsize=12, fontweight='bold')

            cbar = plt.colorbar(im, ax=ax)
            cbar.set_label('Int. Intensity (K·m/s)', fontsize=9)

        gridding_method = 'cygrid' if HAS_CYGRID else 'scipy'
        weight_info = f', weighted by {eff_weight_column}' if eff_weight_column else ''
        if weight_channels:
            weight_info += ', chan-weighted by exp(-τ)/T_sys'
        suptitle = f'Integrated Intensity Comparison ({gridding_method}, beam={beamsize_deg*3600:.1f}″{weight_info})'
        if object_filter:
            suptitle += f' — {object_filter}'
        if velocity_range:
            suptitle += f' [{velocity_range[0]:.0f}–{velocity_range[1]:.0f} km/s]'
        fig.suptitle(suptitle, fontsize=13, fontweight='bold')
        fig.tight_layout()

        click.echo(f"\n✓ {n} comparison maps created ({nrows}×{ncols} grid)")
        click.echo(f"  Beam size: {beamsize_deg*3600:.2f}\"")
        if object_filter:
            click.echo(f"  Object filter: {object_filter}")
        if eff_weight_column:
            click.echo(f"  Per-spectrum weighting: {eff_weight_column}")
        if weight_channels:
            click.echo(f"  Per-channel weighting: exp(-tau)/T_sys from TSYS/TAU_SIG calibration spectra")
        if velocity_range:
            click.echo(f"  Velocity range: {velocity_range[0]:.1f} - {velocity_range[1]:.1f} km/s")

        if plot:
            plt.savefig(plot, dpi=150)
            click.echo(f"  Plot saved: {plot}")
        else:
            plt.show()

        click.echo()
        plt.close()

    except FileNotFoundError as e:
        click.echo(click.style(f"Error: {e}", fg="red"), err=True)
        sys.exit(1)
    except ValueError as e:
        click.echo(click.style(f"Error: {e}", fg="red"), err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(click.style(f"Unexpected error: {e}", fg="red"), err=True)
        import traceback
        traceback.print_exc()
        sys.exit(1)


@click.command()
@click.option('--config', type=click.Path(exists=True),
              help='Path to config.toml file.')
@click.option('--fits', 'fits_file', type=click.Path(exists=True), default=None,
              help='Path to an arbitrary input FITS file (takes priority over dataset flags).')
@click.option('--reduced', is_flag=True, default=False,
              help='Use reduced_data.fits (output.reduced_fits from config).')
@click.option('--pcad', is_flag=True, default=False,
              help='Use pca_corrected.fits (output.pcad_fits from config).')
@click.option('--prepared', is_flag=True, default=False,
              help='Use prepared_for_pca.fits (output.prepared_for_pca from config).')
@click.option('--object', type=str, default=None,
              help='Filter by object name. If not specified, uses "object" from config.toml if available.')
@click.option('--beamsize', type=float, default=None,
              help='Beam size in arcseconds for gridding kernel (e.g. --beamsize 14.1). '
                   'If not specified, reads from config [gridding].beamsize_arcsec.')
@click.option('--pixsize', type=float, default=None,
              help='Map pixel size in arcseconds. If not specified, uses beamsize/3 or config [gridding].pixel_size_arcsec.')
@click.option('--pixel-size-arcsec', type=float, default=None,
              help='Map pixel size in arcseconds (e.g. --pixel-size-arcsec 4.7). '
                   'Overrides --pixsize and config. Default is beamsize/3 ≈ 4.7″ for a 14.1″ beam.')
@click.option('--output', type=click.Path(), default=None,
              help='Output FITS file path. If not specified, uses "datacube" from config.toml or ./datacube.fits.')
@click.option('--plot', type=click.Path(), default=None,
              help='Output path for diagnostic plot (e.g., datacube_slices.png). If not specified, plot is shown but not saved.')
@click.option('--n-jobs', type=int, default=-1,
              help='Number of parallel workers for channel gridding. -1 = all CPUs (default), 1 = sequential.')
@click.option('--weight-spectra', type=str, default=None, metavar='COLUMN',
              help='Column name to use for per-spectrum weighting during gridding '
                   '(e.g. --weight-spectra RMS, --weight-spectra RMS_THEORETICAL, or '
                   '--weight-spectra RMSRATIOB). For RMSRATIO* columns a Gaussian transform '
                   'exp(-(v-1)²/0.5²) is applied; for all other columns inverse-variance '
                   'weighting w=1/value² is used (e.g. RMS gives GILDAS-style SIGMA weighting). '
                   'A WEIGHT_MAP extension is written to the output FITS file.')
@click.option('--weight-channels', is_flag=True, default=False,
              help='Weight each channel by exp(-tau)/T_sys from TSYS/TAU_SIG calibration spectra.')
@click.option('--create-weights-datacube', is_flag=True, default=False,
              help='Only effective together with --weight-channels: writes a full 3D '
                   'WEIGHT_CUBE extension (one weight-sum plane per channel) instead of the '
                   'usual 2D WEIGHT_MAP, since per-channel weights vary channel to channel '
                   '(GILDAS\'s xy_map never needs this — its per-spectrum weights are constant '
                   'across channels, hence its .wei file is always 2D). Without '
                   '--weight-channels this flag has no effect: the existing 2D WEIGHT_MAP/'
                   'COVERAGE extensions are written exactly as before whenever --weight-spectra '
                   'is set.')
@click.option('--kernel-fwhm', type=float, default=None,
              help='Gridding kernel FWHM in arcseconds. If not specified, uses the beam size '
                   '(more smoothing, better sensitivity for noise-limited data). Use a value '
                   'smaller than the beam (e.g. beam/3, matching GILDAS xy_map\'s default) to '
                   'trade sensitivity for resolution on well-sampled, high-S/N data instead '
                   '(e.g., --kernel-fwhm 4.7 for a 14.1″ beam).')
@click.option('--telescop', type=str, default=None,
              help='Telescope name written to the FITS header. '
                   'If not specified, reads from config [gridding].telescop (default: IRAM-30M).')
def create_datacube_cmd(config, fits_file, reduced, pcad, prepared, object, beamsize, pixsize, pixel_size_arcsec, output, plot, n_jobs, weight_spectra, weight_channels, create_weights_datacube, kernel_fwhm, telescop):
    """
    Create a full 3D spectral datacube by gridding spectra across spatial and spectral axes.
    
    Grids each velocity channel separately onto a spatial map, creating a proper 3D datacube
    (nvel, dec, ra) with WCS headers. Uses cygrid for optimal Gaussian kernel gridding
    with scipy fallback.
    
    Reads input file from (in order of priority):
    1. --fits path (arbitrary file, takes priority over all flags)
    2. --pcad flag (uses output.pcad_fits from config)
    3. --prepared flag (uses output.prepared_for_pca from config)
    4. --reduced flag (uses output.reduced_fits from config)
    5. [input].fits_file from config.toml otherwise

    Output file path is determined by (in order of priority):
    1. --output command line option
    2. [output].datacube from config.toml
    3. Default: ./datacube.fits

    Gridding parameters read from config.toml [gridding] section if not specified on command line.

    Examples:
        create_datacube --config config.toml
        create_datacube --config config.toml --fits postpcarmsr_1.3.fits --object M51
        create_datacube --config config.toml --pcad --object M51 --weight-spectra RMSRATIOB --weight-channels
        create_datacube --config config.toml --fits postprocessed.fits --weight-spectra RMS_BASELINE_2
        create_datacube --config config.toml --prepared --object M51
        create_datacube --config config.toml --reduced --object M51
        create_datacube --config config.toml --output my_datacube.fits --plot slices.png
        create_datacube --config config.toml --beamsize 14.1 --pixsize 4.7
        create_datacube --fits postpd.fits --beamsize 14.1 --object M51CENTER --telescop IRAM-30M --output cube.fits
    """
    try:
        from oi_zeigt.mapping.gridding import get_gridding_params_from_config
        
        # Load config
        config_path = config if config else None
        cfg = get_config(config_path) if config_path else {}

        # Get beamsize and pixsize from config if not provided on CLI
        # Both --beamsize and --pixsize are in arcseconds; convert to degrees here.
        # --pixel-size-arcsec is an explicit arcsec override for pixsize.
        effective_beamsize = beamsize / 3600.0 if beamsize is not None else None
        if pixel_size_arcsec is not None:
            effective_pixsize = pixel_size_arcsec / 3600.0
        elif pixsize is not None:
            effective_pixsize = pixsize / 3600.0
        else:
            effective_pixsize = None
        beamsize_deg, pixsize_deg = get_gridding_params_from_config(
            config_path=config_path,
            beamsize_deg=effective_beamsize,
            pixsize_deg=effective_pixsize
        )

        # Determine FITS file to process
        output_cfg = cfg.get('output', {})
        if fits_file:
            click.echo(f"Reading FITS file: {fits_file}")
            hdul = fits.open(fits_file)
        elif pcad:
            fits_path = output_cfg.get('pcad_fits', None)
            if not fits_path:
                click.echo(click.style("Error: --pcad flag specified but [output].pcad_fits not defined in config", fg="red"), err=True)
                sys.exit(1)
            click.echo(f"Reading FITS file: {fits_path}")
            hdul = fits.open(fits_path)
        elif prepared:
            fits_path = output_cfg.get('prepared_for_pca', None)
            if not fits_path:
                click.echo(click.style("Error: --prepared flag specified but [output].prepared_for_pca not defined in config", fg="red"), err=True)
                sys.exit(1)
            click.echo(f"Reading FITS file: {fits_path}")
            hdul = fits.open(fits_path)
        elif reduced:
            fits_path = output_cfg.get('reduced_fits', None)
            if not fits_path:
                click.echo(click.style("Error: --reduced flag specified but [output].reduced_fits not defined in config", fg="red"), err=True)
                sys.exit(1)
            click.echo(f"Reading FITS file: {fits_path}")
            hdul = fits.open(fits_path)
        else:
            if not config_path:
                click.echo(click.style("Error: no input file specified. Use --fits or --config.", fg='red'), err=True)
                sys.exit(1)
            hdul = read_fits_from_config(config_path)

        # Determine object filter
        if object is None:
            # Try to get from config
            parameters_cfg = cfg.get('parameters', {})
            obj_from_config = parameters_cfg.get('object', None)
            object_filter = obj_from_config
        else:
            object_filter = object

        # Determine output file path
        if output is None:
            # Try to get from config
            output_cfg = cfg.get('output', {})
            output_from_config = output_cfg.get('datacube', None)
            output_file = output_from_config if output_from_config else './datacube.fits'
        else:
            output_file = output

        # Extract unique PCAPARAM values from input before gridding
        _pcapar_values = []
        for _hdu in hdul:
            if (hasattr(_hdu, 'data') and _hdu.data is not None
                    and hasattr(_hdu.data, 'dtype')
                    and _hdu.data.dtype.names is not None
                    and 'PCAPARAM' in _hdu.data.dtype.names):
                def _s(v): return v.decode('utf-8').strip() if isinstance(v, bytes) else str(v).strip()
                _raw = [_s(v) for v in _hdu.data['PCAPARAM']]
                _pcapar_values = sorted(set(_raw))
                break

        # Create datacube
        from oi_zeigt.mapping.gridding import create_spectral_datacube

        gridding_cfg = cfg.get('gridding', {})
        telescop = telescop if telescop is not None else gridding_cfg.get('telescop', 'IRAM-30M')

        eff_weight_column = weight_spectra  # None or column name string

        datacube, wcs_header, fig = create_spectral_datacube(
            hdul,
            beamsize_deg=beamsize_deg,
            pixsize=pixsize_deg,
            object_filter=object_filter,
            output_file=output_file,
            telescop=telescop,
            n_jobs=n_jobs,
            weight_column=eff_weight_column,
            channel_weights=weight_channels,
            kernel_fwhm_arcsec=kernel_fwhm,
            create_weights_datacube=create_weights_datacube,
        )

        # Write PCA parameters into the cube primary header
        if _pcapar_values:
            with fits.open(output_file, mode='update') as _out:
                for _i, _val in enumerate(_pcapar_values):
                    _key = f'PCAPAR{_i}'
                    _out[0].header[_key] = (_val, 'PCA correction parameters')
                _out.flush()
            click.echo(f"  PCA parameters written to cube header ({len(_pcapar_values)} unique value(s))")

        click.echo(f"\n✓ Spectral datacube created")
        click.echo(f"  Beam size: {beamsize_deg:.4f}°")
        if pixsize_deg:
            click.echo(f"  Pixel size: {pixsize_deg:.6f}°")
        if object_filter:
            click.echo(f"  Object filter: {object_filter}")
        if eff_weight_column:
            click.echo(f"  Per-spectrum weighting: {eff_weight_column}")
        if weight_channels:
            click.echo(f"  Per-channel weighting: exp(-tau)/T_sys from TSYS/TAU_SIG calibration spectra")
            if create_weights_datacube:
                click.echo(f"  Weight output: 3D WEIGHT_CUBE (per-channel weight sum)")
            elif eff_weight_column:
                click.echo(f"  Weight output: 2D WEIGHT_MAP (per-spectrum only — pass --create-weights-datacube for the full 3D cube)")
        elif eff_weight_column:
            click.echo(f"  Weight output: 2D WEIGHT_MAP")
        click.echo(f"  Datacube shape: {datacube.shape[0]} channels × {datacube.shape[1]} × {datacube.shape[2]} pixels")
        click.echo(f"  Output file: {output_file}")

        # Save plot if requested
        if plot:
            plt.savefig(plot, dpi=150)
            click.echo(f"  Plot saved: {plot}")
        plt.close()

        click.echo()

        hdul.close()

    except FileNotFoundError as e:
        click.echo(click.style(f"Error: {e}", fg="red"), err=True)
        sys.exit(1)
    except ValueError as e:
        click.echo(click.style(f"Error: {e}", fg="red"), err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(click.style(f"Unexpected error: {e}", fg="red"), err=True)
        import traceback
        traceback.print_exc()
        sys.exit(1)


@click.command()
@click.option("--fits", "fits_path", type=click.Path(exists=True), required=True,
              help="Input FITS file to split.")
@click.option("--output-dir", type=click.Path(), required=True,
              help="Directory to write chunk files into.")
@click.option("--manifest", "manifest_path", type=click.Path(), default=None,
              help="Write a manifest listing all chunk paths (for combine_fits --input).")
@click.option("--manifest-science", "manifest_science_path", type=click.Path(), default=None,
              help="Write a manifest listing only science chunks (excludes TREC-only scans).")
@click.option("--no-scan-split", is_flag=True, default=False,
              help="Split by (MISSION_ID, TELESCOP) only, skipping the SCAN dimension.")
def split_fits_cmd(fits_path, output_dir, manifest_path, manifest_science_path, no_scan_split):
    """
    Split a large FITS file into one chunk per (MISSION_ID, TELESCOP, SCAN).

    TELESCOP is the receiver/pixel identifier (e.g. LFAV_PX00_S).  Each scan
    contains its own complete calibration rows (TSYS, TAU_SIG, SKYCHOPDIFF, …)
    alongside science rows, so scan-level splitting is safe.

    Two manifests can be written:

    \b
      --manifest         all chunks (including TREC-only scans)
      --manifest-science science chunks only (skip TREC-only scans in pipeline)

    Use --no-scan-split for one chunk per flight × receiver instead.

    Examples:

        split_fits --fits big.fits --output-dir chunks/ \\
            --manifest chunks/manifest.txt \\
            --manifest-science chunks/manifest_science.txt
    """
    try:
        split_by_scan = not no_scan_split
        label = "MISSION_ID / TELESCOP / SCAN" if split_by_scan else "MISSION_ID / TELESCOP"
        click.echo(f"Splitting {fits_path} by {label} ...")
        all_paths, science_paths = split_fits_by_mission(
            fits_path, output_dir,
            manifest_path=manifest_path,
            manifest_science_path=manifest_science_path,
            split_by_scan=split_by_scan,
        )
        n_trec = len(all_paths) - len(science_paths)
        click.echo(click.style(f"✓ Written {len(all_paths)} chunk files to {output_dir}", fg="green"))
        click.echo(f"  Science chunks : {len(science_paths)}")
        click.echo(f"  TREC-only chunks: {n_trec}")
        if manifest_path:
            click.echo(f"  All-chunks manifest    : {manifest_path}")
        if manifest_science_path:
            click.echo(f"  Science manifest       : {manifest_science_path}")
    except Exception as e:
        click.echo(click.style(f"Error: {e}", fg="red"), err=True)
        import traceback; traceback.print_exc()
        sys.exit(1)


@click.command()
@click.option(
    "--input",
    "input_paths",
    type=click.Path(exists=True),
    required=True,
    multiple=True,
    help="Either a single text file listing FITS file paths (one per line), "
         "or one or more FITS files given directly (repeat --input for each "
         "file)."
)
@click.option(
    "--output",
    type=click.Path(),
    required=True,
    help="Output path for the combined FITS file"
)
@click.option(
    "--single-hdu",
    is_flag=True,
    default=False,
    help="Combine all spectra into a single binary table HDU instead of keeping separate extension HDUs"
)
def combine_fits(input_paths, output, single_hdu):
    """
    Combine multiple FITS files into a single FITS file.

    --input accepts either a text file with a list of FITS file paths (one
    per line), or the FITS files themselves, given directly and repeated
    once per file.

    Examples:

        # From a list file: keep all HDUs as separate extensions
        combine_fits --input list.txt --output combined.fits

        # FITS files given directly
        combine_fits --input a.fits --input b.fits --input c.fits --output combined.fits

        # Merge all data into single HDU
        combine_fits --input list.txt --output combined.fits --single-hdu

    A list file can contain comments (lines starting with '#')
    and blank lines, which will be ignored.
    """
    try:
        click.echo("\n" + "="*70)
        click.echo("COMBINING FITS FILES")
        click.echo("="*70 + "\n")

        output_path = Path(output)

        # A single --input whose own extension isn't .fits/.fit is treated as
        # a list file (legacy behavior); anything else (multiple --input, or
        # a single one that's itself a .fits file) is taken as literal FITS
        # file paths to combine directly.
        list_file_path = None
        if len(input_paths) == 1 and Path(input_paths[0]).suffix.lower() not in ('.fits', '.fit'):
            list_file_path = Path(input_paths[0])
            click.echo(f"Reading FITS file list from: {list_file_path}\n")
            with open(list_file_path, 'r') as f:
                files = [line.strip() for line in f if line.strip() and not line.strip().startswith('#')]
        else:
            files = list(input_paths)

        click.echo(f"Found {len(files)} FITS files to combine:\n")
        for i, f in enumerate(files, 1):
            click.echo(f"  {i:2d}. {f}")
        click.echo()

        # Display combination mode
        if single_hdu:
            click.echo(click.style("Mode: Single HDU (all data merged into one binary table)", fg="cyan"))
        else:
            click.echo(click.style("Mode: Multiple HDUs (separate extension HDUs)", fg="cyan"))
        click.echo()

        # Combine the files
        if list_file_path is not None:
            combined_hdul = combine_fits_from_list(list_file_path, output_path, single_hdu=single_hdu)
        else:
            combined_hdul = combine_fits_files(files, single_hdu=single_hdu)
            output_path.parent.mkdir(parents=True, exist_ok=True)
            combined_hdul.writeto(output_path, overwrite=True)

        # Display results
        click.echo(click.style("✓ Successfully combined FITS files!", fg="green"))
        click.echo(f"\nOutput file: {output_path}")
        click.echo(f"Total HDUs: {len(combined_hdul)}\n")

        combined_hdul.info()
        click.echo("\n" + "="*70 + "\n")

        combined_hdul.close()
        
    except FileNotFoundError as e:
        click.echo(click.style(f"Error: {e}", fg="red"), err=True)
        sys.exit(1)
    except ValueError as e:
        click.echo(click.style(f"Error: {e}", fg="red"), err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(click.style(f"Unexpected error: {e}", fg="red"), err=True)
        import traceback
        traceback.print_exc()
        sys.exit(1)


@click.command()
@click.option(
    "--input",
    "input_fits",
    type=click.Path(exists=True),
    required=True,
    help="Path to FITS file containing spectra."
)
@click.option(
    "--output-folder",
    type=click.Path(),
    default=None,
    help="Directory to save cascade plots into (created if missing). One "
         "subdirectory per MISSION_ID is created inside it. If omitted, "
         "nothing is saved — all groups are shown on screen instead, as "
         "panels in a single figure."
)
@click.option(
    "--hdu-index",
    type=int,
    default=None,
    help="HDU index to read (default: first HDU with a SPECTRUM column)."
)
@click.option(
    "--max-figsize-height",
    type=float,
    default=30.0,
    show_default=True,
    help="Cap on figure height in inches, regardless of how many spectra are in a group."
)
@click.option(
    "--object",
    "object_filter",
    default=None,
    help="Only include spectra whose OBJECT contains this substring (case-insensitive). "
         "Default: no filtering, all rows included regardless of OBJECT."
)
def cascade_plots(input_fits, output_folder, hdu_index, max_figsize_height, object_filter):
    """
    Draw cascade (waterfall) plots of spectra, one plot per
    MISSION_ID / SCAN / TELESCOP group.

    Each group is rendered as an image — one row per spectrum, color encoding
    amplitude, with a colorbar (same color-coding convention as pca_correct's
    diagnostic plots) — useful for visually scanning many individual spectra
    in a scan for line shape, strength, or artefacts, without rows overlapping
    regardless of how many spectra are in the group.

    Examples:

        # Save one PNG per group
        cascade_plots --input input.fits --output-folder plots/dir

        # No --output-folder: show all groups as panels on screen instead
        cascade_plots --input input.fits

        cascade_plots --input input.fits --output-folder plots/dir --object M51CENTER
    """
    try:
        from .cascade_plots import draw_cascade_plots

        click.echo("\n" + "="*70)
        click.echo("CASCADE PLOTS")
        click.echo("="*70 + "\n")

        saved_paths = draw_cascade_plots(
            input_fits=input_fits,
            output_folder=output_folder,
            hdu_index=hdu_index,
            max_figsize_height=max_figsize_height,
            object_filter=object_filter,
        )

        if output_folder is not None:
            click.echo(click.style(f"\n✓ Saved {len(saved_paths)} cascade plot(s) to {output_folder}", fg="green"))
        else:
            click.echo(click.style("\n✓ Displayed cascade panels on screen", fg="green"))

    except FileNotFoundError as e:
        click.echo(click.style(f"Error: {e}", fg="red"), err=True)
        sys.exit(1)
    except ValueError as e:
        click.echo(click.style(f"Error: {e}", fg="red"), err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(click.style(f"Unexpected error: {e}", fg="red"), err=True)
        import traceback
        traceback.print_exc()
        sys.exit(1)


@click.command()
@click.option(
    "--config",
    type=click.Path(exists=False),
    default=None,
    help="Path to config.toml file (reads input/output paths and settings)"
)
@click.option(
    "--fits",
    type=click.Path(exists=False),
    default=None,
    help="Input FITS file to prepare. If not given, uses output.reduced_fits from config.toml."
)
@click.option(
    "--output",
    type=click.Path(),
    default=None,
    help="Output FITS file path (default: output.prepared_for_pca from config.toml)"
)
@click.option(
    "--pca-source",
    type=str,
    default=None,
    help="PCA source to filter for (default: pca.pca_source from config.toml, e.g., SKYCHOPDIFF)"
)
@click.option(
    "--object",
    type=str,
    default=None,
    help="Object substring to filter for (default: parameters.object from config.toml, e.g., M51CENTER)"
)
@click.option(
    "--mission-id",
    type=str,
    default=None,
    help="Filter to specific mission ID (e.g., 2016-05-12_GR_F296) for faster testing"
)
@click.option(
    "--scan",
    type=int,
    default=None,
    help="Filter to specific SCAN number (e.g., 13686) for faster testing"
)
@click.option(
    "--aor-id",
    "aor_id",
    type=str,
    default=None,
    help="Keep only spectra whose AOR_ID contains one of these substrings (comma-separated)."
)
@click.option(
    "--fill-telluric-with-noise",
    "fill_noise",
    is_flag=True,
    default=False,
    help="Fill telluric line channels with Gaussian noise (default: off)"
)
@click.option(
    "--baseline",
    "baseline_order",
    type=int,
    default=None,
    metavar="ORDER",
    help="Subtract a polynomial baseline of this order from the science spectra "
         "(OBJECT == the --object value) before writing. The order must be given "
         "here on the command line — there is no default and it is NEVER read from "
         "config.toml; omit the flag to skip baselining. The PCA reference "
         "(--pca-source, e.g. SKYCHOPDIFF) and TSYS/TAU_SIG rows are left untouched. "
         "The emission line is protected by iterative sigma-clipping (no line window "
         "is read from config)."
)
@click.option(
    "--filter-missions",
    "filter_missions",
    is_flag=True,
    default=False,
    help="Apply drop rules from mission_id_parameters.yml (telescope/scan exclusions per mission)."
)
@click.option(
    "--filter-flight",
    "filter_flight",
    type=str,
    multiple=True,
    help="Remove all rows whose MISSION_ID contains this string (e.g. F528). Can be repeated."
)
@click.option(
    "--filter-below",
    "filter_below",
    type=(str, float),
    multiple=True,
    metavar="COLUMN VALUE",
    help="Keep spectra where COLUMN <= VALUE (drop the noisier ones above it), "
         "same semantics as filter_fits --filter-below. Applies to science rows "
         "(OBJECT contains --object) AND the PCA reference (--pca-source, e.g. "
         "SKYCHOPDIFF) so noisy references can be dropped from the basis; TSYS and "
         "TAU_SIG calibration rows are ALWAYS exempt. Rows with NaN in the column "
         "always pass. RMSRATIO and RMSRATIOB fall back to each other if only one is "
         "present. NOTE: SKYCHOPDIFF has NaN RMSRATIO (RMS_THEORETICAL=0), so a "
         "RMSRATIO cut only affects science — use RMS to also cut noisy references. "
         "Can be repeated. Example: --filter-below RMSRATIO 1.3"
)
@click.option(
    "--filter-above",
    "filter_above",
    type=(str, float),
    multiple=True,
    metavar="COLUMN VALUE",
    help="Keep spectra where COLUMN >= VALUE (mirror of --filter-below). "
         "Same row scope (science + PCA reference; calibration exempt) and "
         "NaN/fallback rules. Can be repeated."
)
@click.option(
    "--mission-parameters",
    "mission_parameters",
    type=click.Path(exists=True),
    default=None,
    help="Path to mission parameters YAML file (overrides pca.mission_parameters in config)."
)
def prepare_for_pca(config: Optional[str], fits: Optional[str], output: Optional[str],
                   pca_source: Optional[str], object: Optional[str], mission_id: Optional[str],
                   scan: Optional[int], aor_id: Optional[str], fill_noise: bool,
                   baseline_order: Optional[int],
                   filter_missions: bool, filter_flight: tuple,
                   filter_below: tuple, filter_above: tuple,
                   mission_parameters: Optional[str]):
    """
    Prepare FITS data for PCA analysis.
    
    Filters data to include only specified PCA source and object,
    then fills telluric line regions with Gaussian noise.
    
    Uses configuration from config.toml by default:
    - Input:       output.reduced_fits  (override with --fits)
    - Output:      output.prepared_for_pca  (override with --output)
    - PCA source:  pca.pca_source  (override with --pca-source, e.g. SKYCHOPDIFF, SKYDIFF)
    - Object filter: parameters.object  (override with --object)
    
    Examples:
    
        # Use config.toml defaults
        prepare_for_pca --config config.toml
        
        # Specify custom files
        prepare_for_pca --fits input.fits --output output.fits
        
        # Override config settings
        prepare_for_pca --config config.toml --pca-source SKYCHOPDIFF --object M51CENTER
        
        # Fast testing with single mission (7x speedup)
        prepare_for_pca --config config.toml --mission-id 2016-05-12_GR_F296
        
        # Very fast testing with single scan (38x speedup)
        prepare_for_pca --config config.toml --mission-id 2016-05-12_GR_F296 --scan 13686
    """
    try:
        from .pca_analysis.prepare_for_pca import prepare_for_pca as prepare_func
        
        click.echo("="*70)
        click.echo("Preparing data for PCA analysis")
        click.echo("="*70 + "\n")
        
        prepare_func(
            fits_file=fits,
            output_fits=output,
            config=config,
            pca_source=pca_source,
            object_filter=object,
            mission_id=mission_id,
            scan=scan,
            aor_id=aor_id,
            fill_noise=fill_noise,
            filter_missions=filter_missions,
            filter_flight=list(filter_flight) if filter_flight else None,
            filter_below=list(filter_below) if filter_below else None,
            filter_above=list(filter_above) if filter_above else None,
            mission_params_file=mission_parameters,
            baseline_order=baseline_order,
        )
        
        click.echo("\n" + "="*70)
        click.echo(click.style("✓ Data preparation complete!", fg="green"))
        click.echo("="*70 + "\n")
        
    except FileNotFoundError as e:
        click.echo(click.style(f"Error: {e}", fg="red"), err=True)
        sys.exit(1)
    except ValueError as e:
        click.echo(click.style(f"Error: {e}", fg="red"), err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(click.style(f"Unexpected error: {e}", fg="red"), err=True)
        import traceback
        traceback.print_exc()
        sys.exit(1)


@click.command()
@click.option('--config', type=click.Path(exists=True),
              help='Path to config.toml file.')
@click.option('--pcad', is_flag=True, default=False,
              help='Use pca_corrected.fits from config (default if no dataset flag given).')
@click.option('--clean', is_flag=True, default=False,
              help='Use clean_fits from config instead of pcad.')
@click.option('--prepared', is_flag=True, default=False,
              help='Use prepared_for_pca from config instead of pcad.')
@click.option('--reduced', is_flag=True, default=False,
              help='Use reduced_fits from config instead of pcad.')
@click.option('--fits', 'fits_input', type=click.Path(exists=True), default=None,
              help='Path to input FITS file. Overrides --pcad/--clean/--prepared/--reduced and config.')
@click.option('--output', type=click.Path(), default=None,
              help='Output FITS file path. Defaults to output.post_processed_fits from config.')
@click.option('--refill-telluric-noise', 'refill_telluric', is_flag=True, default=False,
              help='Refill telluric line channels with Gaussian noise at the post-PCA noise level.')
@click.option('--baseline', 'baseline_order', type=str, default=None, metavar='ORDER|auto',
              help='Apply polynomial baseline subtraction to science spectra before computing RMS. '
                   'Pass an integer polynomial order, or "auto" to choose the order per spectrum '
                   '(capped at 3) by the radiometer stopping rule: raise the order until the '
                   'residual line-free RMS drops to ~the theoretical (radiometer) RMS, then stop. '
                   'The line window from config [reduction].line_window is excluded from the fit. '
                   '"auto" also writes a BLORDER column and prints a per-flight order histogram. '
                   'If not given, no baseline is applied.')
@click.option('--smooth', 'smooth_window', type=int, default=None, is_flag=False, flag_value=5,
              help='Apply boxcar smoothing to science spectra (same kernel as reduce_spectra). '
                   'Optional window size (default: 5). E.g. --smooth or --smooth 7. Applied after '
                   '--baseline and before the RMS is measured, so the RMS reflects the smoothed spectra.')
@click.option('--filter-flight', 'filter_flight', type=str, multiple=True,
              help='Remove all rows whose MISSION_ID contains this string (e.g. F528). Can be repeated.')
@click.option('--aor-id', 'aor_id', type=str, default=None,
              help='Keep only rows whose AOR_ID contains one of these substrings (comma-separated, e.g. 506,507).')
@click.option('--filter-missions', 'filter_missions', is_flag=True, default=False,
              help='Apply telescope/scan drop rules from mission_id_parameters.yml '
                   '(same rules as prepare_for_pca --filter-missions).')
@click.option('--mission-parameters', 'mission_parameters', type=click.Path(), default=None,
              help='Path to mission_id_parameters YAML. Overrides config [pca]/[input][mission_parameters].')
@click.option('--whiteness-check', 'whiteness_check', is_flag=True, default=False,
              help='Flag problematic spectra by an Allan-variance whiteness test: decimate the '
                   'line-free baseline channels by a factor k and check that the RMS drops as √k '
                   '(pure thermal noise) rather than slower (correlated baseline ripple/fringes/'
                   'standing waves, which survive binning). Metric R = √k·σ_binned/σ_baseline is ~1 '
                   'for good spectra, >1 for bad. Writes a WHITENESS column (R) and a WHITEFLAG '
                   'column (1=bad, 0=good, -1=non-science), and prints how many spectra were flagged.')
@click.option('--whiteness-bin', 'whiteness_bin', type=int, default=4, show_default=True, metavar='K',
              help='Decimation factor k for the whiteness test.')
@click.option('--whiteness-sigma', 'whiteness_sigma', type=float, default=3.0, show_default=True, metavar='Z',
              help='Empirical cut: flag spectra with R > median(R) + Z·MAD(R) over all science '
                   'spectra (data-adaptive, robust to any mild non-whiteness common to all spectra).')
def post_process_data_cmd(config, pcad, clean, prepared, reduced, fits_input, output, refill_telluric,
                          baseline_order, smooth_window, filter_flight, aor_id, filter_missions, mission_parameters,
                          whiteness_check, whiteness_bin, whiteness_sigma):
    """
    Post-process spectra: optionally baseline-subtract, then compute per-spectrum RMS.

    Reads spectra from the PCA-corrected file (default) or an alternative dataset.
    Optionally applies a polynomial baseline subtraction (--baseline N) and/or boxcar
    smoothing (--smooth [N]) to science spectra (baseline first, then smoothing, as in
    reduce_spectra), excluding the line window from the baseline fit.  Then computes the
    noise RMS in channels outside [reduction].line_window (km/s) from config.toml, and
    writes the result as RMS / RMS_THEORETICAL / RMSRATIOB columns in the output FITS file.

    Input priority (first matching flag wins):
      --fits > --clean > --prepared > --reduced > --pcad (default)

    Examples:
        post_process_data --config config.toml
        post_process_data --config config.toml --pcad --baseline 3
        post_process_data --config config.toml --pcad --baseline auto
        post_process_data --config config.toml --pcad --baseline 1 --smooth 5
        post_process_data --config config.toml --reduced
        post_process_data --config config.toml --output my_post.fits

    Use filter_fits to filter spectra by RMSRATIOB after post-processing.
    """
    try:
        from astropy.table import Table
        from .reduction.core import _extract_spectral_params, _create_velocity_axis

        config_path = config if config else None
        cfg = get_config(config_path) if config_path else {}
        output_cfg = cfg.get('output', {})
        reduction_cfg = cfg.get('reduction', {})
        object_filter = cfg.get('parameters', {}).get('object', None)

        # --baseline accepts an integer polynomial order, or the literal 'auto'
        # for per-spectrum order selection (capped at 3) via the radiometer
        # stopping rule (raise order until residual line-free RMS ≈ theoretical RMS).
        baseline_auto = False
        if baseline_order is not None:
            if str(baseline_order).strip().lower() == 'auto':
                baseline_auto = True
                baseline_order = None
            else:
                try:
                    baseline_order = int(baseline_order)
                except (TypeError, ValueError):
                    click.echo(click.style(
                        f"Error: --baseline must be an integer order or 'auto' (got {baseline_order!r})",
                        fg='red'), err=True)
                    sys.exit(1)

        # --- Determine input file ---
        if fits_input:
            fits_path = fits_input
            label = 'custom'
        elif clean:
            fits_path = output_cfg.get('clean_fits')
            if not fits_path:
                click.echo(click.style("Error: output.clean_fits not in config", fg='red'), err=True)
                sys.exit(1)
            label = 'clean'
        elif prepared:
            fits_path = output_cfg.get('prepared_for_pca')
            if not fits_path:
                click.echo(click.style("Error: output.prepared_for_pca not in config", fg='red'), err=True)
                sys.exit(1)
            label = 'prepared'
        elif reduced:
            fits_path = output_cfg.get('reduced_fits')
            if not fits_path:
                click.echo(click.style("Error: output.reduced_fits not in config", fg='red'), err=True)
                sys.exit(1)
            label = 'reduced'
        else:
            # Default: pcad
            fits_path = output_cfg.get('pcad_fits')
            if not fits_path:
                click.echo(click.style("Error: output.pcad_fits not in config", fg='red'), err=True)
                sys.exit(1)
            label = 'pcad'

        click.echo(f"Reading {label}: {fits_path}")
        hdul = read_fits(fits_path)

        # --- Determine output file ---
        if output is None:
            output = output_cfg.get('post_processed_fits')
            if not output:
                click.echo(click.style("Error: no --output given and output.post_processed_fits not in config", fg='red'), err=True)
                sys.exit(1)

        # --- Get line_window from config (km/s) ---
        line_window = reduction_cfg.get('line_window', None)
        if line_window is None or len(line_window) != 2:
            click.echo(click.style("Error: [reduction].line_window not defined in config (expected [v_min, v_max] in km/s)", fg='red'), err=True)
            sys.exit(1)
        v_min_kms, v_max_kms = float(line_window[0]), float(line_window[1])
        click.echo(f"Line window: [{v_min_kms}, {v_max_kms}] km/s — RMS computed outside this range")

        # --- Build velocity axis ---
        spectral_params = _extract_spectral_params(hdul)

        # Prefer per-spectrum VELOCITY_AXIS column if present (set by reduce_spectra)
        matrix_hdu = None
        for hdu in hdul:
            if (hasattr(hdu, 'data') and hdu.data is not None
                    and hdu.data.dtype.names is not None
                    and 'SPECTRUM' in hdu.data.dtype.names):
                matrix_hdu = hdu
                break
        if matrix_hdu is None:
            raise ValueError("No SPECTRUM column found in FITS file")

        data = matrix_hdu.data

        if filter_flight and 'MISSION_ID' in data.dtype.names:
            mid_col = np.array([
                s.decode().strip() if isinstance(s, bytes) else str(s).strip()
                for s in data['MISSION_ID']
            ])
            keep_mask = np.array([not any(f in mid for f in filter_flight) for mid in mid_col])
            removed_missions = sorted(set(mid_col[~keep_mask].tolist()))
            for mid in removed_missions:
                click.echo(f"  Removed flight: {mid}")
            click.echo(f"Flight filter: removed {int(np.sum(~keep_mask))} rows ({len(data)} → {int(np.sum(keep_mask))})")
            data = data[keep_mask]

        if aor_id is not None and 'AOR_ID' in data.dtype.names:
            terms = [t.strip() for t in aor_id.split(',') if t.strip()]
            aor_id_col = np.array([
                s.decode().strip() if isinstance(s, bytes) else str(s).strip()
                for s in data['AOR_ID']
            ])
            aor_mask = np.array([any(t in a for t in terms) for a in aor_id_col])
            matched = sorted(set(aor_id_col[aor_mask].tolist()))
            click.echo(f"AOR_ID filter {terms}: {int(np.sum(aor_mask))}/{len(data)} rows kept, matched: {matched}")
            data = data[aor_mask]

        if filter_missions:
            import yaml as _yaml
            from pathlib import Path as _Path
            from .pca_analysis.prepare_for_pca import _default_mission_params_file

            if mission_parameters:
                _yml_path = _Path(mission_parameters)
            else:
                _yml_path = (cfg.get('pca', {}).get('mission_parameters') or
                             cfg.get('input', {}).get('mission_parameters'))
                _yml_path = _Path(_yml_path) if _yml_path else _default_mission_params_file()

            if not _yml_path.exists():
                click.echo(click.style(f"Warning: mission YAML not found: {_yml_path} — skipping --filter-missions", fg='yellow'))
            else:
                click.echo(f"Applying mission drop rules from {_yml_path}")
                with open(_yml_path, 'r') as _f:
                    _mparams = _yaml.safe_load(_f) or {}

                mission_ids_col = np.array([
                    s.decode().strip() if isinstance(s, bytes) else str(s).strip()
                    for s in data['MISSION_ID']
                ]) if 'MISSION_ID' in data.dtype.names else None
                scan_col = np.array(data['SCAN']) if 'SCAN' in data.dtype.names else None
                telescop_col = np.array([
                    s.decode().strip() if isinstance(s, bytes) else str(s).strip()
                    for s in data['TELESCOP']
                ]) if 'TELESCOP' in data.dtype.names else None

                if mission_ids_col is not None:
                    keep_mask = np.ones(len(data), dtype=bool)
                    n_before = len(data)

                    def _tele_match(tele):
                        return np.array([t == tele or t.startswith(tele + '_')
                                         for t in telescop_col])

                    for mid, params in _mparams.items():
                        if not params:
                            continue
                        drop_val = params.get('drop')
                        mid_mask = mission_ids_col == mid

                        # drop: flight — remove every row for this mission/flight
                        if drop_val == 'flight' or (
                                isinstance(drop_val, dict) and drop_val.get('flight')):
                            n = int(np.sum(mid_mask & keep_mask))
                            if n:
                                click.echo(f"  drop.flight: {mid} — dropping entire flight "
                                           f"({n} rows)")
                                keep_mask &= ~mid_mask
                            continue

                        drop_cfg = drop_val if isinstance(drop_val, dict) else {}
                        if not drop_cfg:
                            continue
                        if not np.any(mid_mask & keep_mask):
                            continue

                        for tele in (drop_cfg.get('telescope') or []):
                            if telescop_col is not None:
                                drop = keep_mask & mid_mask & _tele_match(tele)
                                n = int(np.sum(drop))
                                if n:
                                    click.echo(f"  drop.telescope: {mid} / {tele} — removing {n} rows")
                                    keep_mask &= ~drop

                        scans_cfg = drop_cfg.get('scans') or {}
                        if scans_cfg and scan_col is not None:
                            for scan_num in (scans_cfg.get('complete') or []):
                                drop = keep_mask & mid_mask & (scan_col == int(scan_num))
                                n = int(np.sum(drop))
                                if n:
                                    click.echo(f"  drop.scans.complete: {mid} / scan {scan_num} — removing {n} rows")
                                    keep_mask &= ~drop
                            for tele, scan_list in ((scans_cfg.get('telescope') or {}).items()):
                                if telescop_col is not None:
                                    for scan_num in (scan_list or []):
                                        drop = keep_mask & mid_mask & _tele_match(tele) & (scan_col == int(scan_num))
                                        n = int(np.sum(drop))
                                        if n:
                                            click.echo(f"  drop.scans.telescope: {mid} / {tele} / scan {scan_num} — removing {n} rows")
                                            keep_mask &= ~drop

                    data = data[keep_mask]
                    click.echo(f"Mission drop rules: {n_before} → {len(data)} rows "
                               f"(removed {n_before - len(data)})")

        n_spectra = len(data)

        velocity_axis_ms = _create_velocity_axis(
            spectral_params['velo_ref'],
            spectral_params['deltav'],
            spectral_params['crpix1_spec'],
            spectral_params['nchans'],
        )
        velocity_axis_kms = velocity_axis_ms / 1000.0

        outside_mask = (velocity_axis_kms < v_min_kms) | (velocity_axis_kms > v_max_kms)
        n_outside = int(np.sum(outside_mask))
        click.echo(f"Channels outside window: {n_outside} / {len(velocity_axis_kms)}")

        # --- Compute measured and theoretical RMS per spectrum ---
        C_MS = 299792458.0  # speed of light in m/s

        spectra       = np.array(data['SPECTRUM'],  dtype=np.float64)

        # --- Per-spectrum theoretical (radiometer) RMS ---
        # Computed up-front (it depends only on the metadata columns, not on the
        # spectrum values) so --baseline auto can use it as the per-spectrum
        # stopping target while selecting the polynomial order.
        tsys_vals     = np.array(data['TSYS'],      dtype=np.float64)   # K
        deltav_vals   = np.array(data['DELTAV'],    dtype=np.float64)   # m/s (actual channel width after decimation)
        restfreq_vals = np.array(data['RESTFREQ'],  dtype=np.float64)   # Hz
        spectime_vals = np.array(data['SPECTIME'],  dtype=np.float64)   # s, on-source
        reftime_vals  = np.array(data['REFTIME'],   dtype=np.float64)   # s, off-source
        tau_vals      = np.array(data['TAU-ATM'],   dtype=np.float64)   # opacity
        elev_vals     = np.array(data['ELEVATIO'],  dtype=np.float64)   # degrees

        # Theoretical radiometer RMS:
        #   σ = T_sys / sqrt(Δν_Hz) * sqrt(1/t_sig + 1/t_ref) * exp(tau/sin(elev))
        # Δν is computed from DELTAV (actual channel width after any decimation) and RESTFREQ.
        # FREQRES reflects the original pre-decimation channel width and must not be used here.
        # Opacity correction skipped when TAU-ATM <= 0 (unphysical values).
        delta_nu = np.abs(deltav_vals) / C_MS * restfreq_vals          # Hz, shape (n_spectra,)
        elev_rad = elev_vals * (np.pi / 180.0)
        sin_elev = np.sin(elev_rad)
        opacity  = np.where(
            (tau_vals > 0) & (sin_elev > 0),
            np.exp(tau_vals / np.where(sin_elev > 0, sin_elev, 1.0)),
            1.0
        )
        valid_base = (delta_nu > 0) & (spectime_vals > 0)
        safe_tref  = np.where(reftime_vals > 0, reftime_vals, np.inf)  # inf → 1/t_ref = 0
        with np.errstate(invalid='ignore', divide='ignore'):
            rms_theoretical = np.where(
                valid_base,
                tsys_vals / np.sqrt(delta_nu) * np.sqrt(1.0 / spectime_vals + 1.0 / safe_tref) * opacity,
                np.nan
            ).astype(np.float32)

        # --- Optional baseline subtraction (before RMS computation) ---
        rms_baseline_values = None
        chosen_orders = None
        if baseline_auto or baseline_order is not None:
            from .reduction.core import _reduce_baseline
            # Only subtract from science rows
            NON_SCIENCE_BL = {'TSYS', 'TAU_SIG', 'SKYCHOPDIFF'}
            if 'OBJECT' in data.dtype.names:
                obj_strs_bl = np.array([
                    s.decode().strip() if isinstance(s, bytes) else str(s).strip()
                    for s in data['OBJECT']
                ])
                sci_bl_mask = np.array([o not in NON_SCIENCE_BL for o in obj_strs_bl])
            else:
                sci_bl_mask = np.ones(len(spectra), dtype=bool)
            # Line window → channel exclusion range for the baseline fit
            line_ch = np.where(~outside_mask)[0]
            bl_window = (int(line_ch[0]), int(line_ch[-1])) if len(line_ch) > 0 else None
            sci_idx = np.where(sci_bl_mask)[0]

            # Channels used to measure the post-baseline residual RMS (line-free).
            bl_outside = np.ones(spectra.shape[1], dtype=bool)
            if bl_window is not None:
                bl_outside[bl_window[0]:bl_window[1] + 1] = False

            rms_baseline_values = np.full(len(spectra), np.nan, dtype=np.float32)

            if baseline_auto:
                # --- Per-spectrum order selection (radiometer stopping rule) ---
                # Fit every science row at orders 0..CAP once, keep each order's
                # residual line-free RMS, then pick the LOWEST order whose residual
                # RMS has dropped to ~the theoretical radiometer RMS (tol×). Picking
                # the lowest such order — rather than the minimum-RMS order — avoids
                # the overfitting trap (residual RMS falls monotonically with order).
                BL_CAP = 3
                BL_TOL = 1.15          # accept an order once residual RMS <= tol × theoretical
                BL_MIN_GAIN = 0.02     # fallback (no theoretical RMS): stop when relative RMS gain < 2%

                sci_spectra = spectra[sci_idx]
                n_sci = len(sci_idx)
                order_stack = []                                   # order_stack[k] = baselined sci spectra at order k
                resid_rms = np.full((BL_CAP + 1, n_sci), np.nan)   # residual line-free RMS per (order, sci row)
                for k in range(BL_CAP + 1):
                    bl_k = _reduce_baseline(sci_spectra, order=k, window=bl_window)
                    order_stack.append(bl_k)
                    ch = bl_k[:, bl_outside]
                    nv = np.sum(~np.isnan(ch), axis=1)
                    with np.errstate(invalid='ignore'):
                        resid_rms[k] = np.where(nv >= 2, np.nanstd(ch, axis=1), np.nan)

                rt_sci = rms_theoretical[sci_idx].astype(np.float64)
                valid_rt = np.isfinite(rt_sci) & (rt_sci > 0)
                # Radiometer rule (vectorised): first order at/under tol × theoretical.
                under = np.isfinite(resid_rms) & (resid_rms <= (BL_TOL * rt_sci)[None, :])
                any_under = under.any(axis=0)
                first_under = under.argmax(axis=0)                 # first True order, else 0
                pick = np.where(valid_rt & any_under, first_under, BL_CAP).astype(np.int16)
                # Fallback for rows lacking a valid theoretical RMS: diminishing-returns rule.
                for j in np.where(~valid_rt)[0]:
                    rr = resid_rms[:, j]
                    p = BL_CAP
                    for k in range(1, BL_CAP + 1):
                        if not (np.isfinite(rr[k - 1]) and rr[k - 1] > 0 and np.isfinite(rr[k])):
                            p = k - 1
                            break
                        if (rr[k - 1] - rr[k]) / rr[k - 1] < BL_MIN_GAIN:
                            p = k - 1
                            break
                    pick[j] = p

                # Apply the chosen order per row and record it.
                for k in range(BL_CAP + 1):
                    sel = np.where(pick == k)[0]
                    if sel.size:
                        spectra[sci_idx[sel]] = order_stack[k][sel]
                chosen_orders = np.full(len(spectra), -1, dtype=np.int16)   # -1 = non-science row
                chosen_orders[sci_idx] = pick
                rms_baseline_values[sci_idx] = resid_rms[pick, np.arange(n_sci)].astype(np.float32)
                data['SPECTRUM'][:] = spectra.astype(data['SPECTRUM'].dtype)

                win_str = f" (excluding channels {bl_window[0]}–{bl_window[1]})" if bl_window else ""
                click.echo(f"Baseline AUTO (cap {BL_CAP}, stop at {BL_TOL}×theoretical RMS) "
                           f"applied to {n_sci} science spectra{win_str}")

                # Per-flight histogram of the selected orders (homogeneity diagnostic).
                gcounts = [int(np.sum(pick == k)) for k in range(BL_CAP + 1)]
                click.echo("  Selected order — overall: " + "  ".join(
                    f"ord{k}: {gcounts[k]} ({100.0 * gcounts[k] / max(n_sci, 1):.1f}%)"
                    for k in range(BL_CAP + 1)))
                if 'MISSION_ID' in data.dtype.names:
                    mid_sci = np.array([
                        s.decode().strip() if isinstance(s, bytes) else str(s).strip()
                        for s in data['MISSION_ID'][sci_idx]
                    ])
                    click.echo("  Selected order per flight:")
                    click.echo("    {:<14}".format("MISSION_ID")
                               + "".join(f"  ord{k}".rjust(7) for k in range(BL_CAP + 1))
                               + "   median")
                    for fl in sorted(set(mid_sci)):
                        oc = pick[mid_sci == fl]
                        counts = [int(np.sum(oc == k)) for k in range(BL_CAP + 1)]
                        med = int(np.median(oc)) if oc.size else 0
                        click.echo("    {:<14}".format(fl)
                                   + "".join(f"{c:7d}" for c in counts)
                                   + f"   {med:4d}")
            else:
                # --- Fixed order (existing behaviour) ---
                _CHUNK = 1000
                for _start in range(0, len(sci_idx), _CHUNK):
                    _idx = sci_idx[_start:_start + _CHUNK]
                    spectra[_idx] = _reduce_baseline(spectra[_idx], order=baseline_order, window=bl_window)
                n_bl = int(np.sum(sci_bl_mask))
                win_str = f" (excluding channels {bl_window[0]}–{bl_window[1]})" if bl_window else ""
                click.echo(f"Baseline order {baseline_order} applied to {n_bl} science spectra{win_str}")
                # Write back so the output table contains the baselined spectra
                data['SPECTRUM'][:] = spectra.astype(data['SPECTRUM'].dtype)
                for _i in sci_idx:
                    _ch = spectra[_i][bl_outside]
                    _nv = int(np.sum(~np.isnan(_ch)))
                    if _nv >= 2:
                        rms_baseline_values[_i] = np.nanstd(_ch.astype(np.float64))

        # --- Optional boxcar smoothing (after baseline, before RMS) ---
        # Same kernel as reduce_spectra (_reduce_smooth → smooth_spectrum, a
        # normalised boxcar via np.convolve). Applied to science rows only —
        # TSYS / TAU_SIG / SKYCHOPDIFF calibration rows are left untouched, matching
        # the --baseline block above. Runs before the RMS is measured so RMS /
        # RMSRATIOB reflect the smoothed spectra.
        if smooth_window is not None:
            from .reduction.core import _reduce_smooth
            NON_SCIENCE_SM = {'TSYS', 'TAU_SIG', 'SKYCHOPDIFF'}
            if 'OBJECT' in data.dtype.names:
                obj_strs_sm = np.array([
                    s.decode().strip() if isinstance(s, bytes) else str(s).strip()
                    for s in data['OBJECT']
                ])
                sci_sm_mask = np.array([o not in NON_SCIENCE_SM for o in obj_strs_sm])
            else:
                sci_sm_mask = np.ones(len(spectra), dtype=bool)
            sci_sm_idx = np.where(sci_sm_mask)[0]
            _CHUNK = 1000
            for _start in range(0, len(sci_sm_idx), _CHUNK):
                _idx = sci_sm_idx[_start:_start + _CHUNK]
                spectra[_idx] = _reduce_smooth(spectra[_idx], window_size=smooth_window)
            click.echo(f"Boxcar smoothing (window {smooth_window}) applied to "
                       f"{len(sci_sm_idx)} science spectra")
            # Write back so the output table and RMS below use the smoothed spectra.
            data['SPECTRUM'][:] = spectra.astype(data['SPECTRUM'].dtype)

        # --- Measured RMS (fully vectorised) on the (possibly baselined) spectra ---
        # Theoretical (radiometer) RMS was already computed up-front, before the
        # baseline block, so --baseline auto could use it as its stopping target.
        outside_spectra = spectra[:, outside_mask]          # (n_spectra, n_outside)
        n_valid_per = np.sum(~np.isnan(outside_spectra), axis=1)
        with np.errstate(invalid='ignore'):
            rms_measured = np.where(
                n_valid_per >= 2,
                np.nanstd(outside_spectra, axis=1),
                np.nan
            ).astype(np.float32)

        with np.errstate(invalid='ignore', divide='ignore'):
            rms_ratio = np.where(
                np.isfinite(rms_measured) & np.isfinite(rms_theoretical) & (rms_theoretical > 0),
                rms_measured / rms_theoretical,
                np.nan
            ).astype(np.float32)

        n_valid_rms = int(np.sum(~np.isnan(rms_measured)))

        # --- Optionally refill telluric channels with post-PCA noise ---
        if refill_telluric:
            from pathlib import Path as _TPath
            from .pca_analysis.prepare_for_pca import load_mission_parameters, _default_mission_params_file
            click.echo("Refilling telluric channels with post-PCA Gaussian noise...")
            np.random.seed(42)

            # Resolve mission parameters YAML — same priority as filter_missions block
            if mission_parameters:
                _tel_yml = _TPath(mission_parameters)
            else:
                _tel_yml_raw = (cfg.get('pca', {}).get('mission_parameters') or
                                cfg.get('input', {}).get('mission_parameters'))
                _tel_yml = _TPath(_tel_yml_raw) if _tel_yml_raw else _default_mission_params_file()
            click.echo(f"  Mission parameters: {_tel_yml}")

            if 'MISSION_ID' not in data.dtype.names:
                click.echo(click.style("Warning: MISSION_ID column not found — skipping telluric refill", fg='yellow'), err=True)
            else:
                mission_ids = np.array([
                    s.decode().strip() if isinstance(s, bytes) else str(s).strip()
                    for s in data['MISSION_ID']
                ])
                objects_col = np.array([
                    s.decode().strip() if isinstance(s, bytes) else str(s).strip()
                    for s in data['OBJECT']
                ])
                NON_SCIENCE = {'TSYS', 'TAU_SIG', 'SKYCHOPDIFF'}
                science_mask = np.array([o not in NON_SCIENCE for o in objects_col])

                unique_missions = sorted(set(mission_ids))
                mission_params_cache = {mid: load_mission_parameters(mid, yaml_file=_tel_yml)
                                        for mid in unique_missions}

                n_refilled = 0

                # Batch by mission_id: telluric mask is the same for all spectra
                # sharing the same mission_id and velocity axis.
                # Only science spectra are refilled — TSYS, TAU_SIG, SKYCHOPDIFF are left unchanged.
                for mid in unique_missions:
                    params = mission_params_cache.get(mid, {})
                    if not params or 'telluric_line_center' not in params:
                        continue

                    idx = np.where((mission_ids == mid) & science_mask)[0]

                    vel_kms = velocity_axis_kms

                    # Compute telluric mask directly from cached params (correct YAML already used)
                    _center = params['telluric_line_center']
                    _width  = params.get('telluric_line_width', 30)
                    telluric_mask = ((vel_kms >= _center - _width / 2.0) &
                                     (vel_kms <= _center + _width / 2.0))
                    if not np.any(telluric_mask):
                        continue

                    n_tel = int(np.sum(telluric_mask))
                    click.echo(f"  {mid}: center={_center} km/s  width={_width} km/s  "
                               f"({n_tel} channels,  {len(idx)} science spectra)")
                    clean_mask = outside_mask & ~telluric_mask

                    # Vectorised: noise level per spectrum from clean channels
                    batch = spectra[idx]                            # (n_batch, n_chan)
                    if np.any(clean_mask):
                        with np.errstate(invalid='ignore'):
                            noise_levels = np.nanstd(batch[:, clean_mask], axis=1)
                    else:
                        with np.errstate(invalid='ignore'):
                            noise_levels = np.nanstd(batch[:, outside_mask], axis=1)

                    for j, (row_idx, nl) in enumerate(zip(idx, noise_levels)):
                        if np.isfinite(nl) and nl > 0:
                            spectra[row_idx][telluric_mask] = np.random.normal(0, nl, n_tel)
                            n_refilled += 1

                n_cal = int(np.sum(~science_mask))
                click.echo(f"  Refilled telluric channels in {n_refilled} science spectra "
                           f"({n_cal} calibration rows skipped)")
                # Write modified spectra back into the data array
                data['SPECTRUM'][:] = spectra.astype(data['SPECTRUM'].dtype)

        # --- Restrict summary to science rows only ---
        # Science rows are those whose OBJECT *contains* object_filter (substring
        # match, consistent with the rest of the codebase, e.g. _print_object_info).
        # This lets a single value like "M51" cover both M51CENTER and M51EDGE for
        # whole-map runs, while a more specific value like "M51CENTER" still
        # matches only the center rows. If object_filter is not set, fall back to
        # excluding known calibration row types (TSYS, TAU_SIG, SKYCHOPDIFF) so
        # calibration spectra do not skew the statistics.
        CAL_TYPES = {'TSYS', 'TAU_SIG', 'SKYCHOPDIFF'}
        objects_col_all = np.array([
            s.decode().strip() if isinstance(s, bytes) else str(s).strip()
            for s in data['OBJECT']
        ])
        if object_filter:
            is_science_all = np.array([object_filter in o for o in objects_col_all], dtype=bool)
        else:
            is_science_all = np.array([o not in CAL_TYPES for o in objects_col_all])
        sci_rms    = rms_measured[is_science_all]
        sci_rmsth  = rms_theoretical[is_science_all]
        sci_ratio  = rms_ratio[is_science_all]
        n_science  = int(np.sum(is_science_all))

        # --- Whiteness (Allan-variance) quality check -------------------------
        # White thermal noise averages down as sqrt(N) under binning; correlated
        # baseline structure (ripple, standing waves, fringes, drifts) does not.
        # Decimate the line-free baseline channels by a factor k and compare the
        # RMS before/after: R = sqrt(k) * sigma_binned / sigma_baseline is ~1 for
        # a white spectrum and > 1 for a spectrum with residual structure. The cut
        # is EMPIRICAL: median(R) + z * MAD(R) over the science population, robust
        # and self-calibrating against any mild non-whiteness common to all rows.
        whiteness_R    = np.full(n_spectra, np.nan, dtype=np.float64)
        whiteflag      = np.full(n_spectra, -1, dtype=np.int16)  # -1 non-science
        n_white_bad    = 0
        white_R_crit   = np.nan
        white_lf       = None   # line-free baseline, kept for the sweep table
        white_nbl      = 0
        if whiteness_check:
            k   = max(2, int(whiteness_bin))
            nbl = int(np.sum(outside_mask))
            lf  = spectra[:, outside_mask].astype(np.float64)     # post-refill baseline
            white_lf, white_nbl = lf, nbl

            def _whiteness_R_for_k(kk):
                """R = sqrt(kk)*sigma_binned/sigma_baseline per spectrum, or None."""
                mm = (lf.shape[1] // kk) * kk
                if mm < kk or nbl < 2 * kk:
                    return None
                b = np.nanmean(lf[:, :mm].reshape(lf.shape[0], -1, kk), axis=2)
                with np.errstate(invalid='ignore', divide='ignore'):
                    _s0 = np.nanstd(lf, axis=1)
                    _sk = np.nanstd(b,  axis=1)
                    return np.where(_s0 > 0, np.sqrt(kk) * _sk / _s0, np.nan)

            m   = (lf.shape[1] // k) * k
            if m >= k and nbl >= 2 * k:
                R = _whiteness_R_for_k(k)
                whiteness_R = R
                # Empirical robust threshold over science spectra only.
                Rsci = R[is_science_all]
                Rsci = Rsci[np.isfinite(Rsci)]
                if len(Rsci) >= 10:
                    med = float(np.median(Rsci))
                    mad = 1.4826 * float(np.median(np.abs(Rsci - med)))
                    white_R_crit = med + whiteness_sigma * mad
                    # Analytic white-null reference (for context only).
                    r_analytic = 1.0 + whiteness_sigma * np.sqrt(k / (2.0 * nbl))
                    bad = is_science_all & np.isfinite(R) & (R > white_R_crit)
                    whiteflag[is_science_all & np.isfinite(R)] = 0
                    whiteflag[bad] = 1
                    n_white_bad = int(np.sum(bad))
                else:
                    click.echo("  [whiteness] too few valid science spectra (<10); skipping cut")
            else:
                click.echo(f"  [whiteness] not enough line-free channels "
                           f"(nbl={nbl}, k={k}); skipping")

        sep = "─" * 60
        click.echo(f"\n{sep}")
        click.echo(f"  Post-process summary  ({n_valid_rms}/{n_spectra} spectra computed, {n_science} science rows)")
        click.echo(sep)

        def _row(label, arr):
            valid = arr[np.isfinite(arr)]
            if len(valid) == 0:
                click.echo(f"  {label:<18}  (no valid values)")
                return
            click.echo(
                f"  {label:<18}"
                f"  min={np.min(valid):8.4f}"
                f"  p25={np.percentile(valid, 25):8.4f}"
                f"  med={np.median(valid):8.4f}"
                f"  p75={np.percentile(valid, 75):8.4f}"
                f"  max={np.max(valid):8.4f}"
                f"  mean={np.mean(valid):8.4f}"
            )

        _row("RMS measured [K]",  sci_rms)
        _row("RMS theoretical [K]", sci_rmsth)
        _row("RMSRATIOB",          sci_ratio)

        # RMSRATIOB distribution bins
        thresholds = [1.0, 1.3, 1.5, 2.0, 3.0]
        valid_ratio = sci_ratio[np.isfinite(sci_ratio)]
        n_valid_ratio = len(valid_ratio)
        click.echo(f"\n  RMSRATIOB distribution ({n_valid_ratio} science spectra with valid ratio):")
        prev = 0.0
        for thr in thresholds:
            n_bin = int(np.sum(valid_ratio <= thr)) - int(np.sum(valid_ratio <= prev))
            cum   = int(np.sum(valid_ratio <= thr))
            pct   = 100.0 * cum / n_valid_ratio if n_valid_ratio > 0 else 0.0
            click.echo(f"    <= {thr:.1f} : {cum:6d} cumulative ({pct:5.1f}%)   +{n_bin} in this bin")
            prev = thr
        n_above = int(np.sum(valid_ratio > thresholds[-1]))
        click.echo(f"    >  {thresholds[-1]:.1f} : {n_above:6d} spectra")

        if whiteness_check and np.isfinite(white_R_crit):
            k   = max(2, int(whiteness_bin))
            nbl = int(np.sum(outside_mask))
            _row("WHITENESS R", whiteness_R[is_science_all])
            r_analytic = 1.0 + whiteness_sigma * np.sqrt(k / (2.0 * nbl))
            click.echo(
                f"\n  Whiteness test (bin k={k}, {nbl} line-free channels):"
                f"\n    R_crit (empirical, med + {whiteness_sigma:g}·MAD) = {white_R_crit:.4f}"
                f"\n    R_crit (analytic white-null, ref only)     = {r_analytic:.4f}"
            )

            # --- Sensitivity sweep: how the flagged fraction moves with the two
            #     knobs, so the tradeoff is visible without a manual re-run. ---
            def _flagged(Rarr, z):
                rs = Rarr[is_science_all]; rs = rs[np.isfinite(rs)]
                if len(rs) < 10:
                    return None, np.nan, 0
                md = float(np.median(rs))
                m9 = 1.4826 * float(np.median(np.abs(rs - md)))
                crit = md + z * m9
                nb = int(np.sum(rs > crit))
                return len(rs), crit, nb

            # Z-sweep at the current k (reuses the already-computed R).
            click.echo(f"\n  Strictness vs --whiteness-sigma  (current k={k}, ★ = current z={whiteness_sigma:g}):")
            click.echo(f"    {'z':>4}  {'R_crit':>8}  {'flagged':>8}  {'pct':>6}")
            for z in [1.5, 2.0, 2.5, 3.0, 3.5, 4.0]:
                nsci, crit, nb = _flagged(whiteness_R, z)
                if nsci is None:
                    continue
                star = " ★" if abs(z - whiteness_sigma) < 1e-9 else ""
                click.echo(f"    {z:>4.1f}  {crit:>8.4f}  {nb:>8d}  {100.0*nb/nsci:>5.1f}%{star}")

            # k-sweep at the current z (recomputes R per k; keep N_bl/k >= 8).
            click.echo(f"\n  Timescale vs --whiteness-bin  (current z={whiteness_sigma:g}, ★ = current k={k}):")
            click.echo(f"    {'k':>4}  {'N_bl/k':>7}  {'R_crit':>8}  {'flagged':>8}  {'pct':>6}")
            for kk in [2, 4, 8, 16, 32]:
                if white_nbl < 2 * kk or (white_lf.shape[1] // kk) < 1:
                    continue
                Rk = _whiteness_R_for_k(kk)
                if Rk is None:
                    continue
                nsci, crit, nb = _flagged(Rk, whiteness_sigma)
                if nsci is None:
                    continue
                ratio = white_nbl / kk
                warn  = " (low)" if ratio < 8 else ""
                star  = " ★" if kk == k else ""
                click.echo(f"    {kk:>4d}  {ratio:>7.1f}  {crit:>8.4f}  {nb:>8d}  "
                           f"{100.0*nb/nsci:>5.1f}%{warn}{star}")
            click.echo("    (N_bl/k < 8 → R gets noisy per-spectrum; treat those rows as indicative)")

        click.echo(sep)

        # --- Build output table: same columns as input, add/replace columns ---
        table = Table(data)
        table['RMS']             = rms_measured
        table['RMS_THEORETICAL'] = rms_theoretical
        table['RMSRATIOB']       = rms_ratio    # ratio of measured RMS to theoretical radiometer RMS

        if whiteness_check:
            w_name = 'WHITENESS'
            _wi = 2
            while w_name in table.colnames:
                w_name = f'WHITENESS_{_wi}'; _wi += 1
            table[w_name] = whiteness_R.astype(np.float32)
            wf_name = 'WHITEFLAG'
            _wi = 2
            while wf_name in table.colnames:
                wf_name = f'WHITEFLAG_{_wi}'; _wi += 1
            table[wf_name] = whiteflag
            click.echo(f"  Written {w_name}: Allan whiteness ratio R (NaN where not computable)")
            click.echo(f"  Written {wf_name}: whiteness flag (1=bad, 0=good, -1=non-science)")

        if baseline_auto or baseline_order is not None:
            col_name = 'RMS_BASELINE'
            while col_name in table.colnames:
                suffix = int(col_name.split('_')[-1]) + 1 if col_name != 'RMS_BASELINE' else 2
                col_name = f'RMS_BASELINE_{suffix}'
            table[col_name] = rms_baseline_values
            click.echo(f"  Written {col_name}: per-spectrum RMS outside line window after baseline subtraction")

            if baseline_auto and chosen_orders is not None:
                bo_name = 'BLORDER'
                _bi = 2
                while bo_name in table.colnames:
                    bo_name = f'BLORDER_{_bi}'
                    _bi += 1
                table[bo_name] = chosen_orders
                click.echo(f"  Written {bo_name}: baseline order chosen per spectrum by --baseline auto "
                           f"(-1 = non-science row)")

        # --- Associate each science spectrum with its TSYS / TAU_SIG calibration row ---
        # Build TSYS_INDEX / TAU_SIG_INDEX here so every downstream stage — PCA or
        # not — inherits the channel-based tau/T_sys association without needing
        # prepare_for_pca (which the no-PCA pipelines skip). The values are absolute
        # row positions in THIS output file's row order; filter_fits remaps them
        # when it drops rows, and combine_fits shifts them by the concat offset, so
        # they stay valid all the way to create_datacube. Built only when TSYS/
        # TAU_SIG rows are actually present (they must be, to point at) — if absent,
        # any pre-existing columns are left untouched rather than clobbered with -1.
        has_cal_rows = bool(np.isin(objects_col_all, ['TSYS', 'TAU_SIG']).any())
        if has_cal_rows and object_filter:
            try:
                from .pca_analysis.prepare_for_pca import build_tau_tsys_indices
                tsys_index, tau_sig_index = build_tau_tsys_indices(data, object_filter)
                table['TSYS_INDEX']    = np.asarray(tsys_index,    dtype=np.int32)
                table['TAU_SIG_INDEX'] = np.asarray(tau_sig_index, dtype=np.int32)
                click.echo(
                    f"  TSYS_INDEX / TAU_SIG_INDEX: linked "
                    f"{int(np.sum(tsys_index >= 0))} / {int(np.sum(tau_sig_index >= 0))} "
                    f"of {n_science} science rows to calibration spectra")
            except Exception as _e:
                click.echo(click.style(
                    f"  Warning: could not build TSYS_INDEX/TAU_SIG_INDEX ({_e}) — "
                    f"channel-based weighting downstream will fall back to uniform",
                    fg='yellow'))
        elif not has_cal_rows:
            click.echo("  TSYS_INDEX / TAU_SIG_INDEX: skipped "
                       "(no TSYS/TAU_SIG calibration rows in this dataset)")
        else:  # cal rows present but no object filter to identify science rows
            click.echo("  TSYS_INDEX / TAU_SIG_INDEX: skipped "
                       "([parameters].object not set — cannot identify science rows)")

        # Preserve header and write output
        primary = hdul[0].copy()
        new_hdu = fits.BinTableHDU(table)
        new_hdu.name = matrix_hdu.name
        for key in matrix_hdu.header:
            if key not in ['NAXIS1', 'NAXIS2', 'TFIELDS'] and key != '':
                try:
                    new_hdu.header[key] = matrix_hdu.header[key]
                except (ValueError, KeyError):
                    pass

        # Write median std of science spectra outside the line window as header
        # keyword. Guard against an empty/all-NaN science set (e.g. a flight with
        # no rows matching object_filter): nanmedian of an empty array is NaN, and
        # FITS headers reject NaN — so skip the keyword rather than crash the run.
        finite_sci = sci_rms[np.isfinite(sci_rms)]
        if finite_sci.size > 0:
            median_std = float(np.median(finite_sci))
            new_hdu.header['STD'] = (median_std, 'Median RMS outside line window [K]')
            click.echo(f"\n  STD header keyword: {median_std:.6f} K")
        else:
            click.echo("\n  STD header keyword: skipped (no science spectra with valid RMS)")

        import os
        os.makedirs(os.path.dirname(os.path.abspath(output)), exist_ok=True)
        fits.HDUList([primary, new_hdu]).writeto(output, overwrite=True)
        click.echo(f"\n✓ Written: {output}")

        if whiteness_check:
            if np.isfinite(white_R_crit):
                pct = 100.0 * n_white_bad / n_science if n_science > 0 else 0.0
                click.echo(sep)
                click.echo(click.style(
                    f"  WHITENESS: {n_white_bad}/{n_science} science spectra "
                    f"({pct:.1f}%) identified as BAD",
                    fg=('yellow' if n_white_bad else 'green'), bold=True))
                click.echo(f"    (R > R_crit = {white_R_crit:.4f}, "
                           f"bin k={max(2, int(whiteness_bin))}, "
                           f"z={whiteness_sigma:g}·MAD)")
                click.echo(sep)
            else:
                click.echo("  WHITENESS: no cut applied (see note above)")

        hdul.close()

    except FileNotFoundError as e:
        click.echo(click.style(f"Error: {e}", fg='red'), err=True)
        sys.exit(1)
    except ValueError as e:
        click.echo(click.style(f"Error: {e}", fg='red'), err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(click.style(f"Unexpected error: {e}", fg='red'), err=True)
        import traceback
        traceback.print_exc()
        sys.exit(1)


@click.command()
@click.option('--fits', 'cube_fits', required=True, type=click.Path(exists=True),
              help='Path to the 3D FITS datacube to collapse.')
@click.option('--velocity-range', type=(float, float), default=None, nargs=2,
              help='Integration range in km/s (e.g., --velocity-range 500 550).')
@click.option('--zoom', 'zoom_size_arcmin', type=float, default=None,
              help='Show a centred zoom panel of this size in arcmin, '
                   'with colour scale recalculated for the zoomed region.')
@click.option('--zoom-x', type=float, default=None,
              help='X pixel coordinate of the zoom region centre (default: map centre).')
@click.option('--zoom-y', type=float, default=None,
              help='Y pixel coordinate of the zoom region centre.')
@click.option('--zoom-ra', type=str, default=None,
              help='RA of the zoom region centre: decimal degrees or hh:mm:ss.s '
                   '(used when --zoom-x/y not given).')
@click.option('--zoom-dec', type=str, default=None,
              help='Dec of the zoom region centre: decimal degrees or dd:mm:ss.s')
@click.option('--region-x', type=float, default=None,
              help='X pixel coordinate of circular extraction aperture centre. '
                   'If omitted with --region-radius, uses map centre.')
@click.option('--region-y', type=float, default=None,
              help='Y pixel coordinate of circular extraction aperture centre.')
@click.option('--region-ra', type=str, default=None,
              help='RA of extraction aperture centre: decimal degrees or hh:mm:ss.s '
                   '(used when --region-x/y not given).')
@click.option('--region-dec', type=str, default=None,
              help='Dec of extraction aperture centre: decimal degrees or dd:mm:ss.s')
@click.option('--region-radius', 'region_radius_arcmin', type=float, default=None,
              help='Radius in arcmin of circular aperture for spectrum extraction. '
                   'Adds a separate spectrum panel and a circle on the map.')
@click.option('--use-wcs', is_flag=True, default=False,
              help='Display map axes in RA/Dec coordinates instead of pixel indices.')
@click.option('--plot', type=click.Path(), default=None,
              help='Save plot to this path.')
@click.option('--fits-output', type=click.Path(), default=None,
              help='Save collapsed 2D map as FITS to this path.')
@click.option('--no-show', is_flag=True, default=False,
              help='Do not open an interactive plot window.')
@click.option('--colormap', default='rainbow', show_default=True,
              help=(
                  'Matplotlib colormap for the intensity map. '
                  'Sequential: inferno, viridis, plasma, magma, cividis, hot, afmhot, gist_heat, YlOrRd, Blues, Greens. '
                  'Diverging: RdBu_r, seismic, bwr, coolwarm, PiYG. '
                  'Perceptual: cubehelix, turbo. '
                  'Classic: jet, rainbow, gray. '
                  'Append "_r" to any name to reverse it (e.g. viridis_r).'
              ))
@click.option('--coverage-threshold', type=float, default=0.0, show_default=True,
              help='Mask edge pixels whose gridding coverage (kernel weight sum) is below this '
                   'fraction of the peak coverage in the map. Requires a COVERAGE extension in '
                   'the cube FITS file (written automatically by create_datacube). '
                   'Default 0 = no coverage masking (show the entire map).')
@click.option('--mask-ra', type=str, default=None,
              help='RA centre of circular display mask (hh:mm:ss or degrees). '
                   'Defaults to map centre when --mask-radius is given without this option.')
@click.option('--mask-dec', type=str, default=None,
              help='Dec centre of circular display mask (dd:mm:ss or degrees). '
                   'Defaults to map centre when --mask-radius is given without this option.')
@click.option('--mask-radius', type=float, default=None,
              help='Radius in arcminutes of the circular display mask. '
                   'Pixels outside this circle are set to NaN in the plot.')
@click.option('--trim-edges', 'trim_edges', default=None,
              help='Peel a border off the irregular coverage footprint, following its '
                   'shape (erosion — NOT a rectangular crop), to shave the ragged '
                   'low-coverage edge while keeping the map outline. A bare number or "N%" '
                   'peels a border of N% of the smaller map axis; a unit suffix peels a '
                   'fixed angular border: "30arcsec"/"30\\"" or "0.5arcmin"/"0.5\'". '
                   'Interior holes are preserved (only the outer edge is trimmed).')
@click.option('--mode', default='moment0', show_default=True,
              type=click.Choice(['moment0', 'peak-intensity'], case_sensitive=False),
              help='Collapse mode. moment0: integrated intensity (K km/s). '
                   'peak-intensity: brightest channel within the velocity window (K). '
                   'Tip: for peak-intensity, use --percentile-clip 0 98 (or --suppress-high) '
                   'to prevent noisy edge pixels from dominating the colour scale.')
@click.option('--peak-range-int', 'peak_range_channels',
              is_flag=False, flag_value=5, default=None, type=int, metavar='N',
              help='Integrate N channels on each side of the per-pixel peak channel '
                   '(window = 2N+1 channels), producing a K km/s map. '
                   'If given without a number, N defaults to 5. '
                   'Overrides --mode.')
@click.option('--peak-smooth', 'peak_smooth_channels', type=float, default=None, metavar='SIGMA',
              help='For peak-intensity / --peak-range-int: locate each pixel\'s peak channel on a '
                   'spectrum Gaussian-smoothed by SIGMA channels along velocity, then read the RAW '
                   '(unsmoothed) cube there. Prevents a narrow residual spike — e.g. an un-refilled '
                   'telluric wing that survived --refill-telluric-noise — from being picked as the '
                   'peak, since smoothing discriminates by line width: the broad real line survives, '
                   'a 1–2 channel spike is crushed. Reported intensities stay in true K; only the '
                   'peak LOCATION uses the smoothed cube. Try SIGMA ~1–2.')
@click.option('--shuffle', is_flag=True, default=False,
              help='Velocity-field shuffle before moment-0: build a smooth per-pixel line-velocity '
                   'field from the cube (moment-1 on high-S/N pixels, smoothed+filled), shift every '
                   'spectrum so its line lands at a common reference velocity, then integrate a NARROW '
                   'window around it. Removes galaxy rotation from the velocity axis so a tight window '
                   'catches the line everywhere → far less integrated noise → faint extended/spiral '
                   'structure survives. Requires --velocity-range (the line window used to measure the '
                   'field). Writes the velocity field to <output>_vfield.fits for inspection. '
                   'Intended with --mode moment0.')
@click.option('--shuffle-snr', type=float, default=5.0, show_default=True, metavar='SNR',
              help='S/N floor for pixels contributing to the shuffle velocity field. Only pixels above '
                   'this get a measured line centroid; the rest inherit a smoothed value from neighbours.')
@click.option('--shuffle-field-smooth', type=float, default=6.0, show_default=True, metavar='PIX',
              help='Gaussian sigma (pixels) for smoothing/gap-filling the shuffle velocity field, so '
                   'faint pixels get a robust shift from their bright neighbours.')
@click.option('--shuffle-window', 'shuffle_window_kms', type=float, default=20.0, show_default=True,
              metavar='KMS',
              help='Half-width (km/s) of the narrow integration window applied AROUND the reference '
                   'velocity after shuffling. Slightly wider than the line to tolerate field error.')
@click.option('--shuffle-ref-velocity', type=float, default=None, metavar='KMS',
              help='Common velocity (km/s) lines are aligned to. Default: median of the measured field. '
                   'Pin to the systemic velocity for run-to-run comparability.')
@click.option('--shuffle-field-output', type=click.Path(), default=None,
              help='Explicit path for the shuffle velocity-field FITS (default: <output>_vfield.fits).')
@click.option('--suppress-negative', is_flag=True, default=False,
              help='Set vmin=0 in the colour scale, clipping negative (noise) values to the '
                   'bottom. Gives physically correct scaling for moment-0 maps where signal '
                   'is always positive.')
@click.option('--suppress-high', type=float, default=None, metavar='VALUE',
              help='Clip pixels with values above VALUE to VALUE. Useful for suppressing '
                   'bright artefacts that compress the colour scale.')
@click.option('--hex-plot', is_flag=True, default=False,
              help='Display the map as a hexagonal scatter plot sampled on a beam/2 hex grid, '
                   'similar to PyStructure.')
@click.option('--contour', is_flag=True, default=False,
              help='Display the map as contour lines instead of a filled image.')
@click.option('--stretch', default='linear', show_default=True,
              type=click.Choice(['linear', 'sqrt', 'power', 'asinh', 'symlog', 'log', 'histeq'],
                                 case_sensitive=False),
              help='Colour stretch for the intensity map. '
                   'sqrt/asinh/symlog compress bright regions and reveal faint detail '
                   '(asinh & symlog stay linear near zero; symlog also compresses negatives). '
                   'power is a tunable power law (see --gamma). log is the most aggressive '
                   'compression. histeq (histogram equalisation) maximises contrast but gives '
                   'a non-quantitative colourbar.')
@click.option('--gamma', type=float, default=1.0, show_default=True,
              help='Exponent for --stretch power: displayed ∝ value**gamma after vmin/vmax '
                   'scaling. gamma<1 lifts faint detail, gamma>1 emphasises bright peaks, '
                   'gamma=1 is linear, gamma=0.5 equals sqrt. Ignored unless --stretch power.')
@click.option('--output-noise-map', 'noise_map_output', type=click.Path(), default=None,
              help='Save the per-pixel noise map (RMS from line-free channels) as a FITS file. '
                   'Requires --velocity-range to identify line-free channels.')
@click.option('--smooth', 'smooth_sigma', type=float, default=None, metavar='SIGMA',
              help='Gaussian smooth the display map before plotting. '
                   'SIGMA is the kernel standard deviation in pixels.')
@click.option('--percentile-clip', type=(float, float), default=None, nargs=2,
              metavar='LO HI',
              help='Percentile range for the colour scale (default: 0 100, full data range). '
                   'Example: --percentile-clip 2 98 to clip outliers.')
@click.option('--snr-threshold', type=float, default=None, metavar='N',
              help='Mask pixels whose integrated intensity SNR is below N: they are '
                   'shown as 0 (the zero colour), like --suppress-negative does for '
                   'negatives, and are excluded from the mean-spectrum panel. '
                   'Requires --velocity-range to compute the noise map from line-free channels.')
@click.option('--drop-rms', type=float, default=None, metavar='LIMIT',
              help='Drop pixels whose per-pixel spectral RMS exceeds LIMIT (K): they are '
                   'shown as 0 (the zero colour), like --snr-threshold, and excluded from '
                   'the mean-spectrum panel. The RMS is measured on the already-gridded/'
                   'convolved cube over the velocity window given by --drop-window (default: '
                   'the line-free channels outside --velocity-range). Use it to remove noisy '
                   'edge/outlier spaxels that a fixed SNR cut leaves in. '
                   'Example: --drop-rms 1.5 --drop-window 0 400.')
@click.option('--drop-window', type=(float, float), default=None, nargs=2,
              metavar='VMIN VMAX',
              help='Velocity window (km/s) over which --drop-rms measures the per-pixel RMS, '
                   'e.g. --drop-window 0 400 to use a line-free stretch. If omitted, the '
                   'channels outside --velocity-range are used. Ignored without --drop-rms.')
@click.option('--align-peaks', is_flag=True, default=False,
              help='Visual aid: before averaging spectra into the zoom/region/full-map '
                   'spectrum panel(s), shift each one so its own peak (sub-channel, found '
                   'within --velocity-range if given) lands on the median peak velocity of '
                   'the stack. Use this to see line shape cleanly when the true centroid '
                   'velocity genuinely differs across spaxels (e.g. a rotation curve) — '
                   'otherwise that spread alone broadens a plain mean. Does not change the '
                   'moment-0 map itself, only these spectrum panels.')
@click.option('--collapse-weights', is_flag=True, default=False,
              help='Collapse the 3D WEIGHT_CUBE extension (written by create_datacube '
                   '--weight-channels --create-weights-datacube) instead of the data, '
                   'summing the per-channel gridding weights over --velocity-range into '
                   'a 2D weight/sensitivity map (units: weight, a plain sum — NOT x dv). '
                   'Writes to --fits-output, or <cube>_weights.fits next to the input '
                   'cube if omitted. All spatial/display options apply as usual (zoom, '
                   'region, --mask-radius, --trim-edges, --suppress-high, --stretch/'
                   '--gamma, --smooth, --percentile-clip, --hex-plot, --contour, --wcs, '
                   '--polygon, --colormap). Line/noise options are ignored with a warning '
                   '(--snr-threshold, --drop-rms, --shuffle, --align-peaks, '
                   '--output-noise-map, --mode/peak modes).')
@click.option('--polygon', is_flag=True, default=False,
              help='Draw a polygon on the map with the MOUSE to restrict the rendered '
                   'area: middle-click adds each vertex, left-click closes it (>=3 '
                   'vertices), right-click undoes the last one. Everything outside is '
                   'blanked. The chosen vertices are printed so the exact selection can '
                   'be reproduced later with --polygon-coords.')
@click.option('--polygon-coords', 'polygon_coords', default=None,
              metavar='"X1,Y1 X2,Y2 ..."',
              help='Reproduce a polygon selection headlessly from explicit pixel '
                   'vertices (space-separated "x,y" pairs), e.g. '
                   '--polygon-coords "120,110 180,115 170,175". No mouse window; '
                   'takes precedence over --polygon.')
def collapse_cube_cmd(cube_fits, velocity_range, zoom_size_arcmin,
                      zoom_x, zoom_y, zoom_ra, zoom_dec,
                      region_x, region_y, region_ra, region_dec,
                      region_radius_arcmin, use_wcs,
                      plot, fits_output, no_show, colormap,
                      coverage_threshold, mask_ra, mask_dec, mask_radius, trim_edges,
                      mode, peak_range_channels, peak_smooth_channels,
                      shuffle, shuffle_snr, shuffle_field_smooth, shuffle_window_kms,
                      shuffle_ref_velocity, shuffle_field_output,
                      suppress_negative, suppress_high,
                      hex_plot, contour, stretch, gamma,
                      noise_map_output, smooth_sigma, percentile_clip,
                      snr_threshold, drop_rms, drop_window,
                      align_peaks, collapse_weights, polygon, polygon_coords):
    """
    Collapse a 3D spectral datacube to a 2D integrated intensity map (moment-0).

    Sums the cube along the velocity axis (K km/s) within an optional velocity
    range.  Weights are already baked into the cube from the gridding step.

    Examples:

        collapse_cube --fits datacube.fits --velocity-range 500 550

        collapse_cube --fits datacube.fits --velocity-range 500 550 \\
            --zoom 5 --region-radius 2.5 --plot moment0.png

        collapse_cube --fits datacube.fits --velocity-range 500 550 \\
            --zoom 5 --region-x 128 --region-y 135 --region-radius 2.5 --use-wcs

        collapse_cube --fits datacube.fits --velocity-range 500 550 \\
            --zoom 5 --zoom-ra 202.4 --zoom-dec 47.2 \\
            --region-ra 202.4 --region-dec 47.2 --region-radius 2.5 --use-wcs
    """
    def _parse_angle(value, is_ra=False):
        """Parse decimal degrees or sexagesimal string → float degrees."""
        if value is None:
            return None
        try:
            return float(value)
        except ValueError:
            from astropy.coordinates import Angle
            import astropy.units as u
            unit = u.hourangle if is_ra else u.deg
            return float(Angle(value, unit=unit).deg)

    try:
        from oi_zeigt.mapping.gridding import collapse_cube

        # In weight-collapse mode, default the FITS output next to the cube so the
        # map lands in the pipeline's OUTDIR even without an explicit --fits-output.
        if collapse_weights and not fits_output:
            from pathlib import Path as _Path
            fits_output = str(_Path(cube_fits).with_suffix('')) + '_weights.fits'

        zoom_ra    = _parse_angle(zoom_ra,   is_ra=True)
        zoom_dec   = _parse_angle(zoom_dec,  is_ra=False)
        region_ra  = _parse_angle(region_ra,  is_ra=True)
        region_dec = _parse_angle(region_dec, is_ra=False)
        mask_ra_deg  = _parse_angle(mask_ra,  is_ra=True)
        mask_dec_deg = _parse_angle(mask_dec, is_ra=False)

        if peak_range_channels is not None:
            mode = 'peak-range-int'

        # --polygon        → interactive mouse selection.
        # --polygon-coords → explicit vertices (headless), takes precedence.
        select_polygon = bool(polygon)
        polygon_vertices = None
        if polygon_coords:
            select_polygon = False
            try:
                polygon_vertices = [
                    (float(p.split(',')[0]), float(p.split(',')[1]))
                    for p in polygon_coords.split()
                ]
            except (ValueError, IndexError):
                click.echo(click.style(
                    "Error: --polygon-coords must be space-separated \"x,y\" pairs, "
                    "e.g. \"120,110 180,115 170,175\" — or use --polygon to draw it "
                    "with the mouse.", fg='red'), err=True)
                sys.exit(1)
            if len(polygon_vertices) < 3:
                click.echo(click.style(
                    "Error: --polygon-coords needs at least 3 vertices", fg='red'), err=True)
                sys.exit(1)

        click.echo(f"Collapsing {'WEIGHT_CUBE' if collapse_weights else 'cube'}: {cube_fits}")
        if velocity_range:
            click.echo(f"  Velocity range: {velocity_range[0]:.1f} – {velocity_range[1]:.1f} km/s")
        if zoom_size_arcmin:
            if zoom_ra is not None and zoom_dec is not None:
                click.echo(f"  Zoom: {zoom_size_arcmin:.1f} arcmin centred on RA={zoom_ra}, Dec={zoom_dec}")
            elif zoom_x is not None and zoom_y is not None:
                click.echo(f"  Zoom: {zoom_size_arcmin:.1f} arcmin centred on pixel ({zoom_x}, {zoom_y})")
            else:
                click.echo(f"  Zoom: {zoom_size_arcmin:.1f} arcmin centred on map centre")
        if region_radius_arcmin:
            if region_ra is not None and region_dec is not None:
                ctr = f"RA={region_ra}, Dec={region_dec}"
            elif region_x is not None and region_y is not None:
                ctr = f"pixel ({region_x}, {region_y})"
            else:
                ctr = "map centre"
            click.echo(f"  Region: r = {region_radius_arcmin:.1f}′ at {ctr}")
        if align_peaks:
            click.echo("  Align peaks: spectrum panels will be peak-aligned before averaging")

        collapsed, header2d, fig = collapse_cube(
            cube_fits=cube_fits,
            velocity_range=velocity_range,
            zoom_size_arcmin=zoom_size_arcmin,
            zoom_x=zoom_x,
            zoom_y=zoom_y,
            zoom_ra=zoom_ra,
            zoom_dec=zoom_dec,
            region_x=region_x,
            region_y=region_y,
            region_ra=region_ra,
            region_dec=region_dec,
            region_radius_arcmin=region_radius_arcmin,
            use_wcs=use_wcs,
            fits_output=fits_output,
            noise_map_output=noise_map_output,
            plot_output=plot,
            colormap=colormap,
            coverage_threshold=coverage_threshold,
            mask_ra=mask_ra_deg,
            mask_dec=mask_dec_deg,
            mask_radius_arcmin=mask_radius,
            trim_edges=trim_edges,
            mode=mode,
            peak_range_channels=peak_range_channels if peak_range_channels is not None else 5,
            peak_smooth_channels=peak_smooth_channels,
            shuffle=shuffle,
            shuffle_snr=shuffle_snr,
            shuffle_field_smooth=shuffle_field_smooth,
            shuffle_window_kms=shuffle_window_kms,
            shuffle_ref_velocity=shuffle_ref_velocity,
            shuffle_field_output=shuffle_field_output,
            suppress_negative=suppress_negative,
            suppress_high=suppress_high,
            hex_plot=hex_plot,
            contour=contour,
            stretch=stretch,
            gamma=gamma,
            smooth_sigma=smooth_sigma,
            percentile_clip=percentile_clip,
            snr_threshold=snr_threshold,
            drop_rms=drop_rms,
            drop_window=drop_window,
            align_peaks=align_peaks,
            collapse_weights=collapse_weights,
            select_polygon=select_polygon,
            polygon=polygon_vertices,
        )

        click.echo(f"\n✓ Collapsed map: {collapsed.shape[1]} × {collapsed.shape[0]} pixels")

        if not no_show:
            plt.show()
        plt.close()

    except (FileNotFoundError, ValueError) as e:
        click.echo(click.style(f"Error: {e}", fg='red'), err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(click.style(f"Unexpected error: {e}", fg='red'), err=True)
        import traceback
        traceback.print_exc()
        sys.exit(1)


@click.command()
@click.option('--fits', 'cube_files', required=False, multiple=True,
              type=click.Path(),
              help='Gridded FITS datacube(s) to collapse and compare. Repeat the flag, '
                   'pass several paths, or use shell globs: --fits a.fits b.fits, '
                   '--fits a.fits --fits b.fits, or --fits *_center_*.fits '
                   '(quote a pattern to expand it here: --fits "*_center_*.fits").')
@click.argument('extra_cube_files', nargs=-1, type=click.Path())
@click.option('--velocity-range', type=(float, float), default=None, nargs=2,
              help='Integration range in km/s (e.g. --velocity-range 450 500). '
                   'Applied to every cube. Required with --shuffle.')
@click.option('--mode', default='moment0', show_default=True,
              type=click.Choice(['moment0', 'peak-intensity'], case_sensitive=False),
              help='Collapse mode for every cube. moment0: integrated intensity (K km/s). '
                   'peak-intensity: brightest channel in the window (K). Overridden by --peak-range-int.')
@click.option('--peak-range-int', 'peak_range_channels',
              is_flag=False, flag_value=5, default=None, type=int, metavar='N',
              help='Integrate N channels either side of each pixel peak (2N+1 window) -> K km/s. '
                   'Without a number, N=5. Overrides --mode.')
@click.option('--peak-smooth', 'peak_smooth_channels', type=float, default=None, metavar='SIGMA',
              help='For peak modes: locate the peak on a spectrum Gaussian-smoothed by SIGMA channels, '
                   'read the RAW cube there (crushes 1-2 channel telluric spikes). Try SIGMA ~1-2.')
@click.option('--shuffle', is_flag=True, default=False,
              help='Velocity-field shuffle before collapse: align each pixel line to a common '
                   'reference velocity, then integrate a narrow window (faint structure survives). '
                   'Requires --velocity-range. Applied per cube independently.')
@click.option('--shuffle-snr', type=float, default=5.0, show_default=True, metavar='SNR',
              help='S/N floor for pixels contributing to the shuffle velocity field.')
@click.option('--shuffle-field-smooth', type=float, default=6.0, show_default=True, metavar='PIX',
              help='Gaussian sigma (pixels) to smooth/gap-fill the shuffle velocity field.')
@click.option('--shuffle-window', 'shuffle_window_kms', type=float, default=20.0, show_default=True,
              metavar='KMS', help='Half-width (km/s) of the aligned integration window after shuffling.')
@click.option('--shuffle-ref-velocity', type=float, default=None, metavar='KMS',
              help='Common velocity (km/s) lines are aligned to. Default: per-cube median field. '
                   'Pin it (e.g. 476.5) for cube-to-cube comparability.')
@click.option('--coverage-threshold', type=float, default=0.0, show_default=True,
              help='Mask pixels below this fraction of peak gridding coverage (needs a COVERAGE '
                   'extension). Default 0 = no coverage masking (show the entire map).')
@click.option('--trim-edges', 'trim_edges', default=None,
              help='Peel a border off each cube coverage footprint (erosion, not a crop). '
                   'Bare number / "N%" = N% of the smaller axis; unit suffix (30arcsec, 0.5arcmin) = angular.')
@click.option('--suppress-negative', is_flag=True, default=False,
              help='Clip negative (noise) pixels to 0 and force vmin=0 (moment0 only).')
@click.option('--suppress-high', type=float, default=None, metavar='VALUE',
              help='Clip pixels above VALUE to VALUE (tame bright artefacts).')
@click.option('--smooth', 'smooth_sigma', type=float, default=None, metavar='SIGMA',
              help='Gaussian-smooth each display map by SIGMA pixels before plotting.')
@click.option('--compare-weights', is_flag=True, default=False,
              help='Compare the 3D WEIGHT_CUBE extensions instead of the science data: '
                   'each panel is the channel sum of gridded weights (Σ K·w). Requires cubes '
                   'built with --weight-channels --create-weights-datacube. Line/peak options '
                   '(--mode peak*, --peak-range-int, --peak-smooth, --shuffle) are ignored.')
@click.option('--use-wcs', is_flag=True, default=False,
              help='Show each panel in RA/Dec (WCS) instead of pixel indices.')
@click.option('--colormap', default='rainbow', show_default=True,
              help='Matplotlib colormap for all panels (append "_r" to reverse).')
@click.option('--stretch', default='sqrt', show_default=True,
              type=click.Choice(['linear', 'sqrt', 'power', 'asinh', 'symlog', 'log'],
                                 case_sensitive=False),
              help='Colour stretch shared by all panels. sqrt/asinh/symlog lift faint detail.')
@click.option('--gamma', type=float, default=1.0, show_default=True,
              help='Exponent for --stretch power (gamma<1 lifts faint, >1 emphasises peaks).')
@click.option('--percentile-clip', type=(float, float), default=None, nargs=2, metavar='LO HI',
              help='Percentile range for the colour scale (e.g. 2 98). Pooled across cubes for the shared scale.')
@click.option('--per-panel-scale', is_flag=True, default=False,
              help='Scale each panel to its OWN range (default: one shared scale + colorbar so '
                   'brightness is directly comparable between cubes).')
@click.option('--titles', default=None,
              help='Comma-separated panel titles, in file order (default: each file stem).')
@click.option('--output', '--plot', 'plot_output', type=click.Path(), default=None,
              help='Save the comparison figure to this path (PNG).')
@click.option('--no-show', is_flag=True, default=False,
              help='Do not open an interactive window (just save/return).')
def compare_maps_cmd(cube_files, extra_cube_files, velocity_range, mode, peak_range_channels,
                     peak_smooth_channels, shuffle, shuffle_snr, shuffle_field_smooth,
                     shuffle_window_kms, shuffle_ref_velocity, coverage_threshold,
                     trim_edges, suppress_negative, suppress_high, smooth_sigma,
                     compare_weights, use_wcs, colormap, stretch, gamma, percentile_clip,
                     per_panel_scale, titles, plot_output, no_show):
    """
    Collapse several gridded datacubes and draw their maps side by side.

    A lightweight, multi-cube companion to collapse_cube: it reuses the same
    collapse machinery (velocity range, peak / peak-range / moment0, --shuffle,
    coverage/edge-trim masks, --suppress-*, --smooth) but drops the interactive
    single-cube panels. The grid size is chosen from the number of cubes, and by
    default all panels share one colour scale + colorbar so map depth is directly
    comparable.

    Examples:

        compare_maps --fits pca.fits --fits nonpca.fits --velocity-range 350 540 \\
            --suppress-negative --stretch sqrt --peak-range-int 20 --output cmp.png

        compare_maps --fits a.fits b.fits c.fits --shuffle --velocity-range 400 560 \\
            --shuffle-ref-velocity 476.5 --shuffle-window 50 --trim-edges 3% --use-wcs

        compare_maps --fits *_center_*.fits --velocity-range 350 540 --peak-range-int 20
    """
    try:
        import os
        import glob as _glob
        from oi_zeigt.mapping.gridding import compare_maps

        # Collect cubes from --fits values, extra positional paths, and any
        # (quoted) shell globs, then expand/dedupe/validate.
        raw_paths = list(cube_files) + list(extra_cube_files)
        expanded = []
        for pat in raw_paths:
            hits = sorted(_glob.glob(pat))
            expanded.extend(hits if hits else [pat])   # keep literal if no match
        seen = set()
        cube_files = []
        for f in expanded:
            if f not in seen:
                seen.add(f)
                cube_files.append(f)

        if not cube_files:
            raise ValueError(
                "No FITS cubes given. Pass paths via --fits (e.g. --fits *_center_*.fits).")
        missing = [f for f in cube_files if not os.path.isfile(f)]
        if missing:
            raise FileNotFoundError("FITS file(s) not found: " + ", ".join(missing))

        title_list = [t.strip() for t in titles.split(',')] if titles else None
        eff_mode = 'weights' if compare_weights else (
            'peak-range-int' if peak_range_channels is not None else mode)

        click.echo(f"Comparing {len(cube_files)} cube(s) [{eff_mode}]:")
        for f in cube_files:
            click.echo(f"  - {f}")

        compare_maps(
            cube_files=cube_files,
            plot_output=plot_output,
            titles=title_list,
            use_wcs=use_wcs,
            colormap=colormap,
            stretch=stretch,
            gamma=gamma,
            percentile_clip=percentile_clip,
            per_panel_scale=per_panel_scale,
            show=not no_show,
            # forwarded to collapse_cube_to_map:
            velocity_range=velocity_range,
            mode=eff_mode,
            peak_range_channels=peak_range_channels if peak_range_channels is not None else 5,
            peak_smooth_channels=peak_smooth_channels,
            shuffle=shuffle,
            shuffle_snr=shuffle_snr,
            shuffle_field_smooth=shuffle_field_smooth,
            shuffle_window_kms=shuffle_window_kms,
            shuffle_ref_velocity=shuffle_ref_velocity,
            coverage_threshold=coverage_threshold,
            trim_edges=trim_edges,
            suppress_negative=suppress_negative,
            suppress_high=suppress_high,
            smooth_sigma=smooth_sigma,
            collapse_weights=compare_weights,
        )
        click.echo(f"\n✓ Comparison of {len(cube_files)} cube(s) done")
        plt.close('all')
    except (FileNotFoundError, ValueError) as e:
        click.echo(click.style(f"Error: {e}", fg='red'), err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(click.style(f"Unexpected error: {e}", fg='red'), err=True)
        import traceback
        traceback.print_exc()
        sys.exit(1)


@click.command()
@click.option(
    "--config",
    type=click.Path(exists=False),
    default=None,
    help="Path to config.toml file"
)
@click.option(
    "--fits",
    "fits_file",
    type=click.Path(exists=False),
    default=None,
    help="Path to FITS file to read directly"
)
@click.option(
    "--reduced", is_flag=True, default=False,
    help="Read from output.reduced_fits in config"
)
@click.option(
    "--clean", is_flag=True, default=False,
    help="Read from output.clean_fits in config"
)
@click.option(
    "--prepared", is_flag=True, default=False,
    help="Read from output.prepared_for_pca in config"
)
@click.option(
    "--pcad", is_flag=True, default=False,
    help="Read from output.pcad_fits in config"
)
@click.option(
    "--sky",
    type=click.Choice(["obs", "OBS", "fit", "FIT"], case_sensitive=False),
    default="fit",
    show_default=True,
    help="Which telluric object type to plot: 'fit' → S-H_FIT (default), 'obs' → S-H_OBS"
)
@click.option(
    "--output",
    type=click.Path(),
    default=None,
    help="Output file for the plot (PNG or PDF). If not given, displays interactively."
)
@click.option(
    "--no-show",
    is_flag=True,
    default=False,
    help="Do not open the interactive plot window (useful with --output)."
)
@click.option(
    "--mask-telluric",
    is_flag=True,
    default=False,
    help="Shade the telluric line region on each panel.  Reads the YAML path from "
         "config [pca][mission_parameters] or [input][mission_parameters], "
         "falling back to the bundled mission_id_parameters.yml."
)
@click.option(
    "--telluric-file",
    type=click.Path(),
    default=None,
    help="Path to mission_id_parameters YAML.  Implies --mask-telluric."
)
def examine_telluric(config, fits_file, reduced, clean, prepared, pcad, sky, output,
                     no_show, mask_telluric, telluric_file):
    """
    Plot the averaged telluric spectrum per mission_id/telescop combination.

    By default plots the average of OBJECT=S-H_FIT spectra.  Pass --sky obs
    to plot OBJECT=S-H_OBS instead.  One panel is produced per unique
    (MISSION_ID, TELESCOP) group found in the file.

    Examples:

        examine_telluric --config config.toml
        examine_telluric --fits data.fits --sky obs --output telluric_obs.png
    """
    import math as _math

    try:
        # Resolve FITS path
        fits_path = fits_file
        if fits_path is None:
            try:
                import tomllib
            except ModuleNotFoundError:
                import tomli as tomllib
            cfg_path = config or "config.toml"
            with open(cfg_path, "rb") as f:
                cfg = tomllib.load(f)
            out_cfg = cfg.get("output", {})
            if pcad:
                fits_path = out_cfg.get("pcad_fits")
            elif prepared:
                fits_path = out_cfg.get("prepared_for_pca")
            elif clean:
                fits_path = out_cfg.get("clean_fits")
            elif reduced:
                fits_path = out_cfg.get("reduced_fits")
            else:
                fits_path = out_cfg.get("prepared_for_pca") or out_cfg.get("reduced_fits")
            if not fits_path:
                click.echo(click.style("Error: could not determine FITS path from config.", fg="red"), err=True)
                sys.exit(1)

        hdul = read_fits(fits_path)
        matrix_hdu = next(
            (hdu for hdu in hdul
             if hasattr(hdu, "data") and hdu.data is not None
             and hdu.data.dtype.names is not None
             and "SPECTRUM" in hdu.data.dtype.names),
            None,
        )
        if matrix_hdu is None:
            click.echo(click.style("Error: no HDU with SPECTRUM column found.", fg="red"), err=True)
            sys.exit(1)

        data = matrix_hdu.data

        # Reconstruct velocity axis
        from .basic_io import reconstruct_velocity_axis
        try:
            vel_ms = reconstruct_velocity_axis(matrix_hdu)
            vel_kms = vel_ms / 1000.0
        except Exception as e:
            click.echo(f"Warning: could not reconstruct velocity axis ({e}); using channel numbers.")
            vel_kms = None

        # Determine target OBJECT string
        target_obj = "S-H_FIT" if sky.lower() == "fit" else "S-H_OBS"

        def _str(v):
            return v.decode().strip() if isinstance(v, bytes) else str(v).strip()

        # Group rows by MISSION_ID only (average all telescops within each mission)
        from collections import defaultdict
        groups = defaultdict(list)
        for i, row in enumerate(data):
            obj = _str(row["OBJECT"])
            if obj != target_obj:
                continue
            mission = _str(row["MISSION_ID"]) if "MISSION_ID" in data.dtype.names else "UNKNOWN"
            groups[mission].append(i)

        if not groups:
            click.echo(click.style(
                f"No spectra with OBJECT='{target_obj}' found in {fits_path}.", fg="red"), err=True)
            sys.exit(1)

        sorted_keys = sorted(groups.keys())
        n_groups = len(sorted_keys)
        click.echo(f"Found {n_groups} missions with OBJECT='{target_obj}'")

        # Load telluric mask ranges from YAML if requested
        telluric_ranges = {}   # mission_id → (v_min_km/s, v_max_km/s)
        if mask_telluric or telluric_file:
            import yaml as _yaml
            from pathlib import Path as _Path
            from oi_zeigt.pca_analysis.prepare_for_pca import _default_mission_params_file

            if telluric_file:
                yaml_path = _Path(telluric_file)
            else:
                # Read path from config: [pca] then [input], then bundled fallback.
                # Load cfg now if it wasn't already loaded (e.g. --fits was given directly).
                yaml_path = None
                _cfg_for_mask = cfg if 'cfg' in dir() else {}
                if not _cfg_for_mask and config:
                    try:
                        import tomllib as _tl
                    except ModuleNotFoundError:
                        import tomli as _tl
                    with open(config, 'rb') as _f:
                        _cfg_for_mask = _tl.load(_f)
                yaml_path = (_cfg_for_mask.get('pca', {}).get('mission_parameters') or
                             _cfg_for_mask.get('input', {}).get('mission_parameters'))
                if yaml_path:
                    yaml_path = _Path(yaml_path)
                else:
                    yaml_path = _default_mission_params_file()

            if yaml_path.exists():
                click.echo(f"Telluric mask: reading from {yaml_path}")
                with open(yaml_path, 'r') as _f:
                    _mission_params = _yaml.safe_load(_f) or {}
                for mission_key in sorted_keys:
                    for yaml_key, params in _mission_params.items():
                        if yaml_key and mission_key in yaml_key:
                            if params and 'telluric_line_center' in params:
                                center = float(params['telluric_line_center'])
                                width  = float(params.get('telluric_line_width', 30))
                                telluric_ranges[mission_key] = (center - width / 2,
                                                                center + width / 2)
                            break
                click.echo(f"  Telluric ranges found for "
                           f"{len(telluric_ranges)}/{n_groups} missions")
            else:
                click.echo(click.style(
                    f"Warning: telluric YAML not found: {yaml_path}", fg="yellow"))

        # Layout: roughly square grid
        n_cols = min(4, n_groups)
        n_rows = _math.ceil(n_groups / n_cols)
        fig, axes = plt.subplots(n_rows, n_cols,
                                 figsize=(5 * n_cols, 3 * n_rows),
                                 squeeze=False)
        fig.suptitle(f"Average {target_obj} per mission", fontsize=11)

        for idx, mission in enumerate(sorted_keys):
            row_idx = idx // n_cols
            col_idx = idx % n_cols
            ax = axes[row_idx][col_idx]

            indices = groups[mission]
            spectra = np.array([data["SPECTRUM"][i] for i in indices], dtype=np.float64)
            avg = np.nanmean(spectra, axis=0)

            x = vel_kms if vel_kms is not None else np.arange(len(avg))
            x_label = "Velocity (km/s)" if vel_kms is not None else "Channel"

            ax.plot(x, avg, linewidth=0.8, color="steelblue")
            ax.set_title(mission, fontsize=8)
            ax.set_xlabel(x_label, fontsize=7)
            ax.set_ylabel("T (K)", fontsize=7)
            ax.tick_params(labelsize=6)
            ax.grid(True, alpha=0.3)
            ax.text(0.97, 0.95, f"n={len(indices)}",
                    transform=ax.transAxes, ha="right", va="top", fontsize=6)

            if mission in telluric_ranges and vel_kms is not None:
                t_min, t_max = telluric_ranges[mission]
                ax.axvspan(t_min, t_max, alpha=0.25, color="orange", zorder=0,
                           label=f"telluric {t_min:.0f}–{t_max:.0f} km/s")
                ax.legend(fontsize=5, loc="upper left")

        # Hide unused axes
        for idx in range(n_groups, n_rows * n_cols):
            axes[idx // n_cols][idx % n_cols].set_visible(False)

        fig.tight_layout()

        if output:
            fig.savefig(output, dpi=150, bbox_inches="tight")
            click.echo(f"✓ Saved plot to {output}")
        if not no_show:
            plt.show()
        plt.close(fig)

    except (FileNotFoundError, ValueError) as e:
        click.echo(click.style(f"Error: {e}", fg="red"), err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(click.style(f"Unexpected error: {e}", fg="red"), err=True)
        import traceback
        traceback.print_exc()
        sys.exit(1)


@click.command()
@click.option('--fits', 'fits_path', type=click.Path(exists=True), default=None,
              help='Path to FITS file to inspect.')
@click.option('--config', type=click.Path(exists=True), default=None,
              help='Config file (uses output.pcad_fits if --fits not given).')
def print_pca_parameters(fits_path, config):
    """Print unique PCAPARAM values found in a FITS file."""
    if fits_path is None and config is None:
        click.echo("Error: provide --fits or --config.", err=True)
        sys.exit(1)

    if fits_path is None:
        cfg = get_config(config)
        fits_path = cfg.get('output', {}).get('pcad_fits')
        if not fits_path:
            click.echo("Error: config has no output.pcad_fits.", err=True)
            sys.exit(1)
        if not Path(fits_path).exists():
            click.echo(f"Error: file not found: {fits_path}", err=True)
            sys.exit(1)

    click.echo(f"Reading: {fits_path}")

    def _s(v):
        return v.decode('utf-8').strip() if isinstance(v, bytes) else str(v).strip()

    unique_vals = []
    counts = []
    source_label = ''

    with fits.open(fits_path, memmap=True) as hdul:
        # Case 1: binary table with PCAPARAM column (pca_correct output)
        data = None
        for hdu in hdul:
            if (hasattr(hdu, 'data') and hdu.data is not None
                    and hasattr(hdu.data, 'dtype')
                    and hdu.data.dtype.names is not None
                    and 'PCAPARAM' in hdu.data.dtype.names):
                data = hdu.data
                break

        if data is not None:
            raw = np.array([_s(v) for v in data['PCAPARAM']])
            unique_vals, counts = np.unique(raw, return_counts=True)
            source_label = f'{len(raw)} spectra'
        else:
            # Case 2: datacube with PCAPAR* header keywords
            hdr = hdul[0].header
            hdr_vals = [hdr[k] for k in hdr if k.startswith('PCAPAR')]
            if hdr_vals:
                unique_vals = hdr_vals
                counts = [None] * len(hdr_vals)
                source_label = 'cube header'
            else:
                click.echo("No PCAPARAM column or PCAPAR* header keywords found in this file.")
                return

    n_label = f"across {source_label}"
    click.echo(f"\nFound {len(unique_vals)} unique PCAPARAM value(s) {n_label}:\n")

    param_labels = {
        'nc':    'number of components',
        'vc':    'variance cutoff',
        'nr':    'noise ratio cutoff',
        'sk':    'smoothing kernel size',
        'cc':    'cut coefficients',
        'gnr':   'global noise ratio cutoff',
        'lk':    'line kernel size',
        'ls':    'line cutoff std (σ)',
        'lstsq': 'least squares',
    }

    for val, count in zip(unique_vals, counts):
        prefix = f"[{count} spectra]" if count is not None else "[cube header]"
        click.echo(f"  {prefix}  {val}")
        tokens = val.split()
        for token in tokens:
            if '=' in token:
                k, v = token.split('=', 1)
                label = param_labels.get(k, k)
                click.echo(f"    {k:5s}  ({label}) = {v}")
        click.echo()


@click.command()
@click.option('--image', 'image_cube', required=True, type=click.Path(exists=True),
              help='Gridded FITS datacube drawn as the filled colour map. Its celestial '
                   'grid is the reference: the contour cube is resampled onto it.')
@click.option('--contour', 'contour_cube', required=True, type=click.Path(exists=True),
              help='Gridded FITS datacube drawn as contours over the image cube.')
@click.option('--velocity-range', type=(float, float), default=None, nargs=2,
              help='Integration range in km/s (e.g. --velocity-range 0 300). Applied to BOTH '
                   'cubes on their own velocity axes, so cubes whose channel grids are offset '
                   'by a fractional channel still integrate the same physical window.')
@click.option('--mode', default='moment0', show_default=True,
              type=click.Choice(['moment0', 'peak-intensity'], case_sensitive=False),
              help='Collapse mode for both cubes. Overridden by --peak-range-int.')
@click.option('--peak-range-int', 'peak_range_channels',
              is_flag=False, flag_value=5, default=None, type=int, metavar='N',
              help='Integrate N channels either side of each pixel peak (2N+1 window) -> K km/s. '
                   'Without a number, N=5. Overrides --mode.')
@click.option('--peak-smooth', 'peak_smooth_channels', type=float, default=None, metavar='SIGMA',
              help='For peak modes: locate the peak on a spectrum smoothed by SIGMA channels, '
                   'read the RAW cube there.')
@click.option('--shuffle', is_flag=True, default=False,
              help='Velocity-field shuffle before collapse. Requires --velocity-range. '
                   'Applied to each cube independently.')
@click.option('--shuffle-snr', type=float, default=5.0, show_default=True, metavar='SNR',
              help='S/N floor for pixels contributing to the shuffle velocity field.')
@click.option('--shuffle-field-smooth', type=float, default=6.0, show_default=True, metavar='PIX',
              help='Gaussian sigma (pixels) to smooth/gap-fill the shuffle velocity field.')
@click.option('--shuffle-window', 'shuffle_window_kms', type=float, default=20.0, show_default=True,
              metavar='KMS', help='Half-width (km/s) of the aligned window after shuffling.')
@click.option('--shuffle-ref-velocity', type=float, default=None, metavar='KMS',
              help='Common velocity (km/s) lines are aligned to. Pin it for comparability.')
@click.option('--coverage-threshold', type=float, default=0.0, show_default=True,
              help='Mask pixels below this fraction of peak gridding coverage.')
@click.option('--trim-edges', 'trim_edges', default=None,
              help='Peel a border off each coverage footprint (erosion, not a crop).')
@click.option('--suppress-negative', is_flag=True, default=False,
              help='Clip negative (noise) pixels to 0 and force vmin=0 (moment0 only).')
@click.option('--suppress-high', type=float, default=None, metavar='VALUE',
              help='Clip pixels above VALUE to VALUE.')
@click.option('--smooth', 'smooth_sigma', type=float, default=None, metavar='SIGMA',
              help='Gaussian-smooth both display maps by SIGMA pixels before plotting.')
@click.option('--levels', default=None,
              help='Explicit contour levels in map units, comma-separated (e.g. "5,10,20,40"). '
                   'Overrides --level-fractions / --n-levels.')
@click.option('--level-fractions', default=None,
              help='Contour levels as comma-separated fractions of the contour map peak '
                   '(e.g. "0.2,0.4,0.6,0.8").')
@click.option('--n-levels', type=int, default=8, show_default=True,
              help='Number of contour levels, evenly spaced over 20-90%% of the contour peak.')
@click.option('--regrid-kernel', 'regrid_kernel_arcsec', type=float, default=None, metavar='ARCSEC',
              help='Gaussian FWHM for the cygrid resampling of the contour map. Default: one '
                   'target pixel. These maps are already beam-convolved, so re-gridding with '
                   'the full beam would smooth the contours relative to the image beneath them.')
@click.option('--contour-color', default='white', show_default=True,
              help='Matplotlib colour for the contour lines.')
@click.option('--contour-linewidth', type=float, default=0.9, show_default=True,
              help='Contour line width.')
@click.option('--label-contours', is_flag=True, default=False,
              help='Write the level value inline on each contour.')
@click.option('--colormap', default='rainbow', show_default=True,
              help='Matplotlib colormap for the image layer (append "_r" to reverse).')
@click.option('--stretch', default='linear', show_default=True,
              type=click.Choice(['linear', 'sqrt', 'power', 'asinh', 'symlog', 'log'],
                                case_sensitive=False),
              help='Colour stretch for the image layer.')
@click.option('--gamma', type=float, default=1.0, show_default=True,
              help='Exponent for --stretch power.')
@click.option('--percentile-clip', type=(float, float), default=None, nargs=2, metavar='LO HI',
              help='Percentile range for the colour scale (e.g. 2 98).')
@click.option('--titles', default=None,
              help='Comma-separated "image,contour" labels (default: each file name).')
@click.option('--output', '--plot', 'plot_output', type=click.Path(), default=None,
              help='Save the overlay figure to this path (PNG).')
@click.option('--no-show', is_flag=True, default=False,
              help='Do not open an interactive window (just save/return).')
def overlay_maps_cmd(image_cube, contour_cube, velocity_range, mode, peak_range_channels,
                     peak_smooth_channels, shuffle, shuffle_snr, shuffle_field_smooth,
                     shuffle_window_kms, shuffle_ref_velocity, coverage_threshold,
                     trim_edges, suppress_negative, suppress_high, smooth_sigma,
                     levels, level_fractions, n_levels, regrid_kernel_arcsec,
                     contour_color, contour_linewidth, label_contours,
                     colormap, stretch, gamma, percentile_clip, titles,
                     plot_output, no_show):
    """
    Overlay two gridded datacubes on one axes: one as colour, one as contours.

    A single-axes companion to compare_maps, which draws cubes side by side.
    Both cubes are collapsed with identical settings (so the layers show the
    same quantity over the same window), then the contour cube is resampled
    onto the image cube's celestial grid with cygrid — the same engine
    create_datacube grids with.

    That resampling is the point of the command. Two cubes can share a pixel
    scale and projection yet still differ in NAXIS and sit a fraction of a
    pixel apart (different CRPIX/CRVAL), in which case contouring one directly
    over the other draws them silently misaligned.

    Examples:

        overlay_maps --image nonpca_cube.fits --contour dr_cube.fits \\
            --velocity-range 0 300 --output overlay.png

        overlay_maps --image a.fits --contour b.fits --velocity-range 350 540 \\
            --level-fractions 0.3,0.5,0.7,0.9 --stretch sqrt --suppress-negative
    """
    try:
        from oi_zeigt.mapping.gridding import overlay_maps

        def _floats(raw):
            if not raw:
                return None
            try:
                return [float(v) for v in str(raw).replace(' ', '').split(',') if v]
            except ValueError:
                raise click.BadParameter(f"could not parse '{raw}' as comma-separated numbers")

        effective_mode = 'peak-range-int' if peak_range_channels is not None else mode.lower()
        title_list = [t.strip() for t in titles.split(',')] if titles else None

        overlay_maps(
            image_cube=image_cube,
            contour_cube=contour_cube,
            velocity_range=velocity_range,
            mode=effective_mode,
            peak_range_channels=peak_range_channels if peak_range_channels is not None else 5,
            peak_smooth_channels=peak_smooth_channels,
            shuffle=shuffle,
            shuffle_snr=shuffle_snr,
            shuffle_field_smooth=shuffle_field_smooth,
            shuffle_window_kms=shuffle_window_kms,
            shuffle_ref_velocity=shuffle_ref_velocity,
            coverage_threshold=coverage_threshold,
            trim_edges=trim_edges,
            suppress_negative=suppress_negative,
            suppress_high=suppress_high,
            smooth_sigma=smooth_sigma,
            levels=_floats(levels),
            level_fractions=_floats(level_fractions),
            n_levels=n_levels,
            regrid_kernel_arcsec=regrid_kernel_arcsec,
            contour_color=contour_color,
            contour_linewidth=contour_linewidth,
            label_contours=label_contours,
            colormap=colormap,
            stretch=stretch.lower(),
            gamma=gamma,
            percentile_clip=percentile_clip,
            titles=title_list,
            plot_output=plot_output,
            show=not no_show,
        )
    except Exception as e:
        click.echo(f"Error: {e}", err=True)
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    print_fits_info()
