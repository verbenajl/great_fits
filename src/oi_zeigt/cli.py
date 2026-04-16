"""
Command-line interface for OI ZEIGT.
"""

import sys
from pathlib import Path
from typing import Optional

import click
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

from .basic_io import read_fits_from_config, read_fits, get_config, combine_fits_from_list
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
            if hasattr(hdu, 'data') and hdu.data is not None:
                if 'SPECTRUM' in hdu.data.dtype.names:
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
            if hasattr(hdu, 'data') and hdu.data is not None:
                if 'SPECTRUM' in hdu.data.dtype.names:
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
            if hasattr(hdu, 'data') and hdu.data is not None:
                if 'SPECTRUM' in hdu.data.dtype.names:
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
    help="Keep rows where COLUMN >= VALUE (filter out below). "
         "Rows with NaN in the column always pass. Can be repeated."
)
@click.option(
    "--filter-above",
    type=(str, float),
    multiple=True,
    metavar="COLUMN VALUE",
    help="Keep rows where COLUMN <= VALUE (filter out above). "
         "Rows with NaN in the column always pass. Can be repeated."
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
def filter_fits(config: Optional[str], fits: Optional[str], object: Optional[str],
                nan_threshold: float, output_clean: Optional[str],
                output_rejected: Optional[str], remove: Optional[str],
                remove_values: tuple, apply_only_to_object: bool, filter_zero: bool,
                filter_below: tuple, filter_above: tuple,
                spectrum_peak_threshold: Optional[float], filter_tau: bool,
                filter_flight: tuple, filter_object_exact: Optional[str],
                filter_out_object_exact: Optional[str],
                exclude_obsmode: Optional[str]):
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

        # Keep only spectra with RMSRATIOB >= 2.0
        filter_fits --config config.toml --filter-below RMSRATIOB 2.0

        # Keep only spectra with 1.3 <= RMSRATIOB <= 3.0
        filter_fits --config config.toml --filter-below RMSRATIOB 1.3 --filter-above RMSRATIOB 3.0
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

        if spectrum_peak_threshold is not None:
            click.echo(f"Spectrum peak filter: reject science spectra with any |channel| > {spectrum_peak_threshold}")
        if filter_tau:
            click.echo("TAU filter: reject science spectra linked to TAU_SIG with channels outside [0.001, 1.0]")
        if filter_flight:
            click.echo(f"Flight filter: removing all rows with MISSION_ID containing: {', '.join(filter_flight)}")

        obj_exact_list = [s.strip() for s in filter_object_exact.split(',')] if filter_object_exact else None
        obj_out_exact_list = [s.strip() for s in filter_out_object_exact.split(',')] if filter_out_object_exact else None

        if obj_exact_list:
            click.echo(f"Object exact keep filter: keeping only OBJECT in: {', '.join(obj_exact_list)}")
        if obj_out_exact_list:
            click.echo(f"Object exact remove filter: removing rows with OBJECT in: {', '.join(obj_out_exact_list)}")

        exclude_obsmode_list = [s.strip() for s in exclude_obsmode.split(',')] if exclude_obsmode else None
        if exclude_obsmode_list:
            click.echo(f"OBSMODE exclusion filter: removing rows with OBSMODE in: {', '.join(exclude_obsmode_list)}")

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
            filter_object_exact=obj_exact_list,
            filter_out_object_exact=obj_out_exact_list,
            exclude_obsmode=exclude_obsmode_list,
        )
        
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
            if hasattr(hdu, "data") and hdu.data is not None and hasattr(hdu.data, "dtype"):
                if "MISSION_ID" in hdu.data.dtype.names:
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
def reduce_spectra_cmd(config, fits, clean, unblank, baseline, baseline_order, baseline_window,
                       smooth, decimate, extract, output):
    """
    Perform spectral reduction with selected methods.
    
    Applies reduction methods in sequence:
    1. Unblank (--unblank): Fill NaN values using interpolation
    2. Extract (--extract): Trim spectrum to velocity range from config [reduction].extract
    3. Baseline subtraction (--baseline): Remove polynomial baseline
    4. Smoothing (--smooth [N]): Apply boxcar smoothing with window N (default: 5).
    5. Decimation (--decimate): Keep every Nth channel using the --smooth window size as N; VELOCITY_AXIS is updated accordingly.

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
                baseline_order = baseline_from_config if baseline_from_config else 1
                try:
                    baseline_order = int(baseline_order)
                except Exception:
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
                   '(e.g. --weight-spectra RMSRATIOB or --weight-spectra RMS_BASELINE_2). '
                   'For RMSRATIO* columns a Gaussian transform exp(-(v-1)²/0.5²) is applied; '
                   'for all other columns the raw values are used directly as weights. '
                   'A WEIGHT_MAP extension is written to the output FITS file.')
@click.option('--weight-channels', is_flag=True, default=False,
              help='Weight each channel by exp(-tau)/T_sys from TSYS/TAU_SIG calibration spectra.')
@click.option('--kernel-fwhm', type=float, default=None,
              help='Gridding kernel FWHM in arcseconds. If not specified, uses the beam size. '
                   'Use a value smaller than the beam to minimize resolution degradation '
                   '(e.g., --kernel-fwhm 4.7 or --kernel-fwhm 7.0 for a 14″ beam).')
@click.option('--telescop', type=str, default=None,
              help='Telescope name written to the FITS header. '
                   'If not specified, reads from config [gridding].telescop (default: IRAM-30M).')
def create_datacube_cmd(config, fits_file, reduced, pcad, prepared, object, beamsize, pixsize, pixel_size_arcsec, output, plot, n_jobs, weight_spectra, weight_channels, kernel_fwhm, telescop):
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
        )

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
    type=click.Path(exists=True),
    required=True,
    help="Path to text file containing list of FITS files (one per line)"
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
def combine_fits(input, output, single_hdu):
    """
    Combine multiple FITS files into a single FITS file.
    
    Takes a text file with a list of FITS file paths (one per line)
    and combines them into a single FITS file.
    
    Examples:
    
        # Default: keep all HDUs as separate extensions
        combine_fits --input list.txt --output combined.fits
        
        # Merge all data into single HDU
        combine_fits --input list.txt --output combined.fits --single-hdu
    
    The input file can contain comments (lines starting with '#')
    and blank lines, which will be ignored.
    """
    try:
        click.echo("\n" + "="*70)
        click.echo("COMBINING FITS FILES")
        click.echo("="*70 + "\n")
        
        input_path = Path(input)
        output_path = Path(output)
        
        # Read the list file and display files to be combined
        click.echo(f"Reading FITS file list from: {input_path}\n")
        
        with open(input_path, 'r') as f:
            files = [line.strip() for line in f if line.strip() and not line.strip().startswith('#')]
        
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
        combined_hdul = combine_fits_from_list(input_path, output_path, single_hdu=single_hdu)
        
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
def prepare_for_pca(config: Optional[str], fits: Optional[str], output: Optional[str],
                   pca_source: Optional[str], object: Optional[str], mission_id: Optional[str],
                   scan: Optional[int], aor_id: Optional[str], fill_noise: bool,
                   filter_missions: bool, filter_flight: tuple):
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
@click.option('--baseline', 'baseline_order', type=int, default=None,
              help='Apply polynomial baseline subtraction of this order to science spectra before '
                   'computing RMS. The line window from config [reduction].line_window is excluded '
                   'from the fit. If not given, no baseline is applied.')
@click.option('--filter-flight', 'filter_flight', type=str, multiple=True,
              help='Remove all rows whose MISSION_ID contains this string (e.g. F528). Can be repeated.')
def post_process_data_cmd(config, pcad, clean, prepared, reduced, fits_input, output, refill_telluric, baseline_order, filter_flight):
    """
    Post-process spectra: optionally baseline-subtract, then compute per-spectrum RMS.

    Reads spectra from the PCA-corrected file (default) or an alternative dataset.
    Optionally applies a polynomial baseline subtraction (--baseline N) to science
    spectra, excluding the line window from the fit.  Then computes the noise RMS in
    channels outside [reduction].line_window (km/s) from config.toml, and writes the
    result as RMS / RMS_THEORETICAL / RMSRATIOB columns in the output FITS file.

    Input priority (first matching flag wins):
      --fits > --clean > --prepared > --reduced > --pcad (default)

    Examples:
        post_process_data --config config.toml
        post_process_data --config config.toml --pcad --baseline 3
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
            if hasattr(hdu, 'data') and hdu.data is not None:
                if 'SPECTRUM' in hdu.data.dtype.names:
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

        # --- Optional baseline subtraction (before RMS computation) ---
        rms_baseline_values = None
        if baseline_order is not None:
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
            _CHUNK = 1000
            for _start in range(0, len(sci_idx), _CHUNK):
                _idx = sci_idx[_start:_start + _CHUNK]
                spectra[_idx] = _reduce_baseline(spectra[_idx], order=baseline_order, window=bl_window)
            n_bl = int(np.sum(sci_bl_mask))
            win_str = f" (excluding channels {bl_window[0]}–{bl_window[1]})" if bl_window else ""
            click.echo(f"Baseline order {baseline_order} applied to {n_bl} science spectra{win_str}")
            # Write back so the output table contains the baselined spectra
            data['SPECTRUM'][:] = spectra.astype(data['SPECTRUM'].dtype)

            # Compute per-spectrum RMS outside the line window (same as reduce_spectra)
            rms_baseline_values = np.full(len(spectra), np.nan, dtype=np.float32)
            if bl_window is not None:
                bl_outside = np.ones(spectra.shape[1], dtype=bool)
                bl_outside[bl_window[0]:bl_window[1] + 1] = False
            else:
                bl_outside = np.ones(spectra.shape[1], dtype=bool)
            for _i in sci_idx:
                _ch = spectra[_i][bl_outside]
                _nv = int(np.sum(~np.isnan(_ch)))
                if _nv >= 2:
                    rms_baseline_values[_i] = np.nanstd(_ch.astype(np.float64))

        tsys_vals     = np.array(data['TSYS'],      dtype=np.float64)   # K
        deltav_vals   = np.array(data['DELTAV'],    dtype=np.float64)   # m/s (actual channel width after decimation)
        restfreq_vals = np.array(data['RESTFREQ'],  dtype=np.float64)   # Hz
        spectime_vals = np.array(data['SPECTIME'],  dtype=np.float64)   # s, on-source
        reftime_vals  = np.array(data['REFTIME'],   dtype=np.float64)   # s, off-source
        tau_vals      = np.array(data['TAU-ATM'],   dtype=np.float64)   # opacity
        elev_vals     = np.array(data['ELEVATIO'],  dtype=np.float64)   # degrees

        # --- Compute measured and theoretical RMS (fully vectorised) ---

        # Measured RMS: nanstd over outside-window channels for each spectrum
        outside_spectra = spectra[:, outside_mask]          # (n_spectra, n_outside)
        n_valid_per = np.sum(~np.isnan(outside_spectra), axis=1)
        with np.errstate(invalid='ignore'):
            rms_measured = np.where(
                n_valid_per >= 2,
                np.nanstd(outside_spectra, axis=1),
                np.nan
            ).astype(np.float32)

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

        with np.errstate(invalid='ignore', divide='ignore'):
            rms_ratio = np.where(
                np.isfinite(rms_measured) & np.isfinite(rms_theoretical) & (rms_theoretical > 0),
                rms_measured / rms_theoretical,
                np.nan
            ).astype(np.float32)

        n_valid_rms = int(np.sum(~np.isnan(rms_measured)))

        # --- Optionally refill telluric channels with post-PCA noise ---
        if refill_telluric:
            from .pca_analysis.prepare_for_pca import get_telluric_indices, load_mission_parameters
            click.echo("Refilling telluric channels with post-PCA Gaussian noise...")

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

                unique_missions = list(set(mission_ids))
                mission_params_cache = {mid: load_mission_parameters(mid) for mid in unique_missions}

                n_refilled = 0

                # Batch by mission_id: telluric mask is the same for all spectra
                # sharing the same mission_id and velocity axis.
                # Only science spectra are refilled — TSYS, TAU_SIG, SKYCHOPDIFF are left unchanged.
                for mid in unique_missions:
                    params = mission_params_cache.get(mid, {})
                    if not params or 'telluric_line_center' not in params:
                        continue

                    idx = np.where((mission_ids == mid) & science_mask)[0]

                    # Velocity axis reconstructed from VELOCITY/DELTAV/CRPIX1 (same for all spectra)
                    vel_kms = velocity_axis_kms

                    telluric_mask = get_telluric_indices(mid, vel_kms)
                    if telluric_mask is None or not np.any(telluric_mask):
                        continue

                    n_tel = int(np.sum(telluric_mask))
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

                click.echo(f"  Refilled telluric channels in {n_refilled} / {n_spectra} spectra")
                # Write modified spectra back into the data array
                data['SPECTRUM'][:] = spectra.astype(data['SPECTRUM'].dtype)

        # --- Restrict summary to science rows only ---
        # Science rows are those whose OBJECT matches object_filter.
        # If object_filter is not set, fall back to excluding known calibration
        # row types (TSYS, TAU_SIG, SKYCHOPDIFF) so calibration spectra do not
        # skew the statistics.
        CAL_TYPES = {'TSYS', 'TAU_SIG', 'SKYCHOPDIFF'}
        objects_col_all = np.array([
            s.decode().strip() if isinstance(s, bytes) else str(s).strip()
            for s in data['OBJECT']
        ])
        if object_filter:
            is_science_all = objects_col_all == object_filter
        else:
            is_science_all = np.array([o not in CAL_TYPES for o in objects_col_all])
        sci_rms    = rms_measured[is_science_all]
        sci_rmsth  = rms_theoretical[is_science_all]
        sci_ratio  = rms_ratio[is_science_all]
        n_science  = int(np.sum(is_science_all))

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

        click.echo(sep)

        # --- Build output table: same columns as input, add/replace columns ---
        table = Table(data)
        table['RMS']             = rms_measured
        table['RMS_THEORETICAL'] = rms_theoretical
        table['RMSRATIOB']       = rms_ratio    # ratio of measured RMS to theoretical radiometer RMS

        if baseline_order is not None:
            col_name = 'RMS_BASELINE'
            while col_name in table.colnames:
                suffix = int(col_name.split('_')[-1]) + 1 if col_name != 'RMS_BASELINE' else 2
                col_name = f'RMS_BASELINE_{suffix}'
            table[col_name] = rms_baseline_values
            click.echo(f"  Written {col_name}: per-spectrum RMS outside line window after baseline subtraction")

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

        # Write median std of science spectra outside the line window as header keyword
        median_std = float(np.nanmedian(sci_rms[np.isfinite(sci_rms)]))
        new_hdu.header['STD'] = (median_std, 'Median RMS outside line window [K]')
        click.echo(f"\n  STD header keyword: {median_std:.6f} K")

        import os
        os.makedirs(os.path.dirname(os.path.abspath(output)), exist_ok=True)
        fits.HDUList([primary, new_hdu]).writeto(output, overwrite=True)
        click.echo(f"\n✓ Written: {output}")

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
@click.option('--colormap', default='inferno', show_default=True,
              help=(
                  'Matplotlib colormap for the intensity map. '
                  'Sequential: inferno, viridis, plasma, magma, cividis, hot, afmhot, gist_heat, YlOrRd, Blues, Greens. '
                  'Diverging: RdBu_r, seismic, bwr, coolwarm, PiYG. '
                  'Perceptual: cubehelix, turbo. '
                  'Classic: jet, rainbow, gray. '
                  'Append "_r" to any name to reverse it (e.g. viridis_r).'
              ))
@click.option('--coverage-threshold', type=float, default=0.3, show_default=True,
              help='Mask edge pixels whose gridding coverage (kernel weight sum) is below this '
                   'fraction of the peak coverage in the map. Requires a COVERAGE extension in '
                   'the cube FITS file (written automatically by create_datacube). '
                   '0 = no coverage masking.')
@click.option('--mask-ra', type=str, default=None,
              help='RA centre of circular display mask (hh:mm:ss or degrees). '
                   'Defaults to map centre when --mask-radius is given without this option.')
@click.option('--mask-dec', type=str, default=None,
              help='Dec centre of circular display mask (dd:mm:ss or degrees). '
                   'Defaults to map centre when --mask-radius is given without this option.')
@click.option('--mask-radius', type=float, default=None,
              help='Radius in arcminutes of the circular display mask. '
                   'Pixels outside this circle are set to NaN in the plot.')
@click.option('--suppress-negative', is_flag=True, default=False,
              help='Set vmin=0 in the colour scale, clipping negative (noise) values to the '
                   'bottom. Gives physically correct scaling for moment-0 maps where signal '
                   'is always positive.')
@click.option('--hex-plot', is_flag=True, default=False,
              help='Display the map as a hexagonal scatter plot sampled on a beam/2 hex grid, '
                   'similar to PyStructure.')
@click.option('--contour', is_flag=True, default=False,
              help='Display the map as contour lines instead of a filled image.')
@click.option('--stretch', default='linear', show_default=True,
              type=click.Choice(['linear', 'sqrt', 'asinh', 'log'], case_sensitive=False),
              help='Colour stretch for the intensity map. '
                   'sqrt/asinh compress bright regions and reveal faint detail. '
                   'log is the most aggressive compression.')
def collapse_cube_cmd(cube_fits, velocity_range, zoom_size_arcmin,
                      zoom_x, zoom_y, zoom_ra, zoom_dec,
                      region_x, region_y, region_ra, region_dec,
                      region_radius_arcmin, use_wcs,
                      plot, fits_output, no_show, colormap,
                      coverage_threshold, mask_ra, mask_dec, mask_radius,
                      suppress_negative, hex_plot, contour, stretch):
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

        zoom_ra    = _parse_angle(zoom_ra,   is_ra=True)
        zoom_dec   = _parse_angle(zoom_dec,  is_ra=False)
        region_ra  = _parse_angle(region_ra,  is_ra=True)
        region_dec = _parse_angle(region_dec, is_ra=False)
        mask_ra_deg  = _parse_angle(mask_ra,  is_ra=True)
        mask_dec_deg = _parse_angle(mask_dec, is_ra=False)

        click.echo(f"Collapsing cube: {cube_fits}")
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
            plot_output=plot,
            colormap=colormap,
            coverage_threshold=coverage_threshold,
            mask_ra=mask_ra_deg,
            mask_dec=mask_dec_deg,
            mask_radius_arcmin=mask_radius,
            suppress_negative=suppress_negative,
            hex_plot=hex_plot,
            contour=contour,
            stretch=stretch,
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
def examine_telluric(config, fits_file, reduced, clean, prepared, pcad, sky, output, no_show):
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
    with fits.open(fits_path, memmap=True) as hdul:
        data = None
        for hdu in hdul:
            if hasattr(hdu, 'data') and hdu.data is not None and hasattr(hdu.data, 'dtype'):
                if 'PCAPARAM' in hdu.data.dtype.names:
                    data = hdu.data
                    break
        if data is None:
            click.echo("No PCAPARAM column found in this file.")
            return

    def _s(v):
        return v.decode('utf-8').strip() if isinstance(v, bytes) else str(v).strip()

    raw = np.array([_s(v) for v in data['PCAPARAM']])
    unique_vals, counts = np.unique(raw, return_counts=True)

    click.echo(f"\nFound {len(unique_vals)} unique PCAPARAM value(s) across {len(raw)} spectra:\n")

    param_labels = {
        'vc': 'variance cutoff',
        'nr': 'noise ratio cutoff',
        'sk': 'smoothing kernel size',
        'cc': 'cut coefficients',
        'gnr': 'global noise ratio cutoff',
        'lk': 'line kernel size',
        'lstsq': 'least squares',
    }

    for val, count in zip(unique_vals, counts):
        click.echo(f"  [{count} spectra]  {val}")
        tokens = val.split()
        for token in tokens:
            if '=' in token:
                k, v = token.split('=', 1)
                label = param_labels.get(k, k)
                click.echo(f"    {k:4s}  ({label}) = {v}")
        click.echo()


if __name__ == "__main__":
    print_fits_info()
