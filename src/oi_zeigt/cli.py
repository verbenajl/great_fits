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
                            average_spectra_from_config)
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


def plot_sample_spectra(config: Optional[str], fits: Optional[str], reduced: bool, clean: bool, prepared: bool,
                       object: Optional[str], num_spectra: int, output: Optional[str]):
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
            if config or (reduced or clean or prepared):
                try:
                    import tomllib
                except ModuleNotFoundError:
                    import tomli as tomllib
                
                config_path = config or "config.toml"
                with open(config_path, 'rb') as f:
                    cfg = tomllib.load(f)
                    config_data = cfg
                    
                    # Handle output file flags
                    if reduced or clean or prepared:
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
        
        # Sample spectra
        sample_size = min(num_spectra, len(matching_indices))
        sampled_indices = np.random.choice(matching_indices, size=sample_size, replace=False)
        sampled_indices = np.sort(sampled_indices)
        original_indices = sampled_indices.tolist()
        
        spectra_to_plot = [data[i] for i in sampled_indices]
        
        click.echo(f"Plotting {sample_size} spectra for object '{object}'")
        click.echo(f"Original FITS row indices (for reference): {original_indices}\n")
        
        # Check if VELOCITY_AXIS column exists, or create one on-the-fly
        velocity_axis_from_fits = None
        has_velocity_axis_column = 'VELOCITY_AXIS' in data.dtype.names
        
        if not has_velocity_axis_column:
            # Try to create velocity axis on-the-fly
            nchans = data['SPECTRUM'][0].shape[0]
            velocity_axis_from_fits = _create_velocity_axis_from_fits(matrix_hdu, nchans)
            if velocity_axis_from_fits is not None:
                click.echo("✓ Created velocity axis on-the-fly from FITS parameters")
            else:
                click.echo("  No velocity axis column; using channel indices for x-axis")
        else:
            click.echo("✓ Using velocity axis from FITS column")
        
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
            
            # Determine x-axis: priority is VELOCITY_AXIS column, then on-the-fly creation
            x_axis = None
            x_label = "Channel"
            
            if has_velocity_axis_column:
                try:
                    velocity_axis = spectrum_data['VELOCITY_AXIS']
                    if velocity_axis is not None and len(velocity_axis) == len(spectrum):
                        x_axis = velocity_axis / 1000.0  # Convert m/s to km/s
                        x_label = "Velocity (km/s)"
                except (IndexError, TypeError):
                    pass
            elif velocity_axis_from_fits is not None:
                # Use the on-the-fly created axis
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
            
            ax.set_title(
                f"{nan_indicator} Row {orig_idx}: {obj_name} ({nan_frac:.1%} NaN)",
                fontsize=10,
                color=title_color,
                weight='bold'
            )
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
def filter_fits(config: Optional[str], fits: Optional[str], object: Optional[str],
                nan_threshold: float, output_clean: Optional[str], 
                output_rejected: Optional[str], remove: Optional[str],
                remove_values: tuple, apply_only_to_object: bool, filter_zero: bool):
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
        
        if output_rejected is None and config:
            try:
                output_rejected = config_data.get("output", {}).get("rejected_fits")
            except (NameError, KeyError):
                pass
        
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
            filter_zero_spectra=filter_zero
        )
        
        # Get statistics before closing
        from astropy.io import fits as fits_lib
        hdul_clean = fits_lib.open(clean_path)
        hdul_rejected = fits_lib.open(rejected_path)
        
        n_clean = len(hdul_clean[1].data) if len(hdul_clean) > 1 else 0
        n_rejected = len(hdul_rejected[1].data) if len(hdul_rejected) > 1 and hdul_rejected[1].data is not None else 0
        
        hdul_clean.close()
        hdul_rejected.close()
        
        hdul.close()
        
        # Print results
        click.echo("="*70)
        click.echo(click.style("✓ Successfully created FITS files", fg="green"))
        click.echo("="*70)
        click.echo(f"\nClean FITS file: {clean_path}")
        click.echo(f"  Records: {n_clean}")
        click.echo(f"  (Contains all non-{object} objects + filtered {object} spectra)")
        click.echo(f"  (Only M51 spectra with < {nan_threshold:.1%} NaN channels)")
        click.echo(f"\nRejected FITS file: {rejected_path}")
        click.echo(f"  Records: {n_rejected}")
        
        # Print detailed rejection statistics
        click.echo(f"\n  Rejection breakdown:")
        click.echo(f"    - NaN threshold violations: {stats['rejected_nan']}")
        if filter_zero:
            click.echo(f"    - All-zero spectra: {stats['rejected_zero']}")
        if stats['rejected_removed'] > 0:
            click.echo(f"    - Removed by --remove criteria: {stats['rejected_removed']}")
        
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
@click.option('--smooth', is_flag=True, default=False,
              help='Apply smoothing.')
@click.option('--smooth-window', type=int, default=5,
              help='Smoothing window size (default: 5).')
@click.option('--decimate', type=int, default=None,
              help='Decimate spectra by taking every Nth channel after smoothing (e.g. --decimate 5). VELOCITY_AXIS is updated accordingly.')
@click.option('--output', type=click.Path(), default=None,
              help='Output FITS file path.')
def reduce_spectra_cmd(config, fits, clean, unblank, baseline, baseline_order, baseline_window,
                       smooth, smooth_window, decimate, output):
    """
    Perform spectral reduction with selected methods.
    
    Applies reduction methods in sequence:
    1. Unblank (--unblank): Fill NaN values using interpolation
    2. Baseline subtraction (--baseline): Remove polynomial baseline
    3. Smoothing (--smooth): Apply boxcar smoothing
    4. Decimation (--decimate N): Keep every Nth channel; VELOCITY_AXIS is updated accordingly.
       Typically used together with --smooth --smooth-window N for proper Nyquist sampling.

    More methods can be added in the future.
    
    Input file:
    - By default: Uses [input].fits_file from config
    - With --fits: Uses specified FITS file
    - With --clean: Uses [output].clean_fits from config
    
    If output path is not specified:
    - Uses config [output].reduced_fits if available
    - Falls back to 'reduced_data.fits' with a warning
    
    Examples:
        reduce_spectra --config config.toml --baseline
        reduce_spectra --config config.toml --clean --baseline
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
        
        # Add extraction from config if available
        reduction_cfg = cfg.get('reduction', {})
        extract_cfg = reduction_cfg.get('extract', None)
        if extract_cfg is not None:
            try:
                if isinstance(extract_cfg, (list, tuple)) and len(extract_cfg) == 2:
                    # Extract range from config is in km/s
                    extract_km_s = [float(extract_cfg[0]), float(extract_cfg[1])]
                    extract_m_s = [extract_km_s[0] * 1000.0, extract_km_s[1] * 1000.0]
                    methods['extract'] = extract_m_s
                    methods['extract_mode'] = 'velocity'
            except (ValueError, TypeError, IndexError):
                pass
        
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
        
        if smooth:
            methods['smooth'] = {
                'window_size': smooth_window
            }

        if decimate is not None and decimate > 1:
            methods['decimate'] = {'factor': decimate}

        if not methods:
            click.echo(click.style("No reduction methods selected. Use --unblank, --baseline, and/or --smooth.", 
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
        
        if smooth:
            click.echo(f"    - Smoothing (window={smooth_window})")
            methods_applied = True

        if decimate is not None and decimate > 1:
            click.echo(f"    - Decimation (factor={decimate})")
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
@click.option('--reduced', is_flag=True, default=False,
              help='Use reduced_data.fits instead of the input data.')
@click.option('--output', type=click.Path(), default=None,
              help='Output FITS file path for averaged spectra.')
@click.option('--object', type=str, multiple=True, default=None,
              help='Object(s) to average (e.g., "M51CENTER" or "S-H_FIT"). Can be specified multiple times. If specified, only these objects are averaged and shown together. Defaults to "object" from config.toml if not specified.')
@click.option('--group-by', type=str, default=None,
              help='Column to group by before averaging. If not specified, defaults to "OBJECT" from config or "OBJECT" column.')
@click.option('--no-group', is_flag=True, default=False,
              help='Disable grouping and compute a single global average.')
@click.option('--plot', type=click.Path(), default=None,
              help='Output path for plot (e.g., average_plot.png). If not specified, plot is shown but not saved.')
def average_cmd(config, reduced, output, object, group_by, no_group, plot):
    """
    Average spectra from a FITS file, grouped by object by default.
    
    Reads input file from config.toml [input].fits_file by default, or uses
    [output].reduced_fits if --reduced flag is specified.
    
    By default, computes per-object averages (grouped by OBJECT column).
    Use --object to average only specific object(s) - can be used multiple times.
    Use --no-group to compute a single global average instead.
    
    Object name defaults to [parameters].object from config.toml if not specified.
    
    Outputs:
    - Averaged FITS file with SPECTRUM, STD, COUNT, and RMS columns
    - Plot showing averaged spectra (always displayed, optionally saved)
    
    Examples:
        average --config config.toml
        average --config config.toml --output my_average.fits
        average --config config.toml --reduced
        average --config config.toml --reduced --plot avg.png
        average --config config.toml --object M51CENTER
        average --config config.toml --object M51CENTER --object S-H_FIT
        average --config config.toml --no-group
    """
    try:
        # Load config
        config_path = config if config else None
        cfg = get_config(config_path) if config_path else {}
        
        # Determine FITS file to process
        hdul = None
        if reduced:
            # Use reduced_data.fits from config
            output_cfg = cfg.get('output', {})
            reduced_fits_path = output_cfg.get('reduced_fits', None)
            if not reduced_fits_path:
                click.echo(click.style("Error: --reduced flag specified but [output].reduced_fits not defined in config", fg="red"), err=True)
                sys.exit(1)
            hdul = read_fits(reduced_fits_path)
        else:
            # Use input file from config (default behavior)
            pass  # Let average_spectra_from_config handle reading when hdul is None
        
        # Determine object(s) to process
        if not object:
            # Try to get from config
            parameters_cfg = cfg.get('parameters', {})
            obj_from_config = parameters_cfg.get('object', None)
            if obj_from_config:
                object = (obj_from_config,)  # Convert to tuple for consistency
        
        # If specific object(s) requested, override grouping settings
        if object:
            group_by_col = 'OBJECT' if len(object) > 1 else None
        else:
            # Determine grouping column
            if no_group:
                group_by_col = None
            elif group_by:
                # User explicitly specified
                group_by_col = group_by
            else:
                # Default to OBJECT, optionally override from config
                reduction_cfg = cfg.get('reduction', {})
                group_by_col = reduction_cfg.get('group_by', 'OBJECT')
        
        # If specific object(s) is requested, filter the FITS data first
        if object:
            if hdul is None:
                if reduced:
                    output_cfg = cfg.get('output', {})
                    reduced_fits_path = output_cfg.get('reduced_fits', None)
                    hdul = fits.open(reduced_fits_path)
                else:
                    hdul = read_fits_from_config(config_path)
            
            # Filter to only the specified objects
            matrix_hdu = None
            for idx, hdu in enumerate(hdul):
                if hasattr(hdu, 'data') and hdu.data is not None:
                    if 'SPECTRUM' in hdu.data.dtype.names and 'OBJECT' in hdu.data.dtype.names:
                        matrix_hdu = hdu
                        matrix_hdu_index = idx
                        break
            
            if matrix_hdu is not None:
                data = matrix_hdu.data
                obj_col = data['OBJECT']
                # Find matching objects for ALL requested object filters (case-insensitive substring match)
                # A spectrum matches if it matches ANY of the requested objects
                object_mask = np.zeros(len(obj_col), dtype=bool)
                for requested_obj in object:
                    obj_matches = np.array([
                        requested_obj.upper() in (obj.decode().strip() if isinstance(obj, bytes) else str(obj).strip()).upper()
                        for obj in obj_col
                    ])
                    object_mask |= obj_matches
                
                if not np.any(object_mask):
                    obj_str = ", ".join(object)
                    click.echo(click.style(f"Error: No objects containing '{obj_str}' found in FITS file", fg="red"), err=True)
                    sys.exit(1)
                
                # Create filtered FITS
                filtered_data = data[object_mask]
                new_table = fits.BinTableHDU(filtered_data)
                new_table.name = matrix_hdu.name
                # Copy header
                for key in matrix_hdu.header:
                    if key not in ['NAXIS1', 'NAXIS2', 'TFIELDS'] and key != '':
                        try:
                            new_table.header[key] = matrix_hdu.header[key]
                        except (ValueError, KeyError):
                            pass
                
                filtered_hdul = fits.HDUList([hdul[0].copy(), new_table])
                hdul = filtered_hdul
        
        # Apply averaging
        output_path = average_spectra_from_config(
            config_path=config_path,
            hdul=hdul,
            output_fits=output,
            group_by=group_by_col,
            overwrite=True
        )
        
        click.echo(f"\n✓ Averaging complete.")
        click.echo(f"  Output FITS: {output_path}")
        if group_by_col:
            click.echo(f"  Grouped by: {group_by_col}")
        else:
            click.echo(f"  Global average (no grouping)")
        click.echo()
        
        # Generate and display plot
        from astropy.io import fits as fits_module
        import matplotlib.pyplot as plt
        
        hdul_avg = fits_module.open(output_path)
        data = hdul_avg[1].data
        
        fig, ax = plt.subplots(figsize=(12, 6))
        
        if 'OBJECT' in data.dtype.names:
            # Multiple averaged spectra (grouped)
            for i, obj in enumerate(data['OBJECT']):
                obj_name = obj.decode().strip() if isinstance(obj, bytes) else str(obj).strip()
                spectrum = data['SPECTRUM'][i]
                ax.plot(spectrum, label=obj_name, alpha=0.7, linewidth=1)
            ax.legend(loc='best', fontsize=9)
            ax.set_title('Averaged Spectra (Grouped)')
        else:
            # Single averaged spectrum
            spectrum = data['SPECTRUM'][0]
            ax.plot(spectrum, label='Average', linewidth=1.5)
            ax.fill_between(np.arange(len(spectrum)), 
                           spectrum - data['STD'][0],
                           spectrum + data['STD'][0],
                           alpha=0.3, label='±1σ')
            ax.legend(loc='best')
            ax.set_title('Averaged Spectrum')
        
        ax.set_xlabel('Channel')
        ax.set_ylabel('Flux')
        ax.grid(True, alpha=0.3)
        plt.tight_layout()
        
        # Save plot if requested
        if plot:
            plt.savefig(plot, dpi=150)
            click.echo(f"  Plot saved: {plot}")
        
        # Always show the plot
        plt.show()
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
@click.option('--object', type=str, default=None,
              help='Filter by object name. If not specified, uses "object" from config.toml if available.')
@click.option('--all', 'all_metrics', is_flag=True, default=False,
              help='Include all available quality metrics.')
@click.option('--rmsratio', is_flag=True, default=False,
              help='Include RMSRATIO metric.')
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
def spechistogram_cmd(config, fits, reduced, clean, rejected, prepared, pcad, postfiltered, object, all_metrics, rmsratio, squality, roll_rms_n, mh2o, 
                      tsys, tau_atm, chi_sqr, err_pwv, rms_baseline, bins, plot):
    """
    Generate histograms of multiple spectral quality metrics.
    
    Creates a combined multi-panel figure showing histograms of selected quality metrics.
    Each metric shows mean and median lines.
    
    Input file priority: --fits option > output flags (--clean, --rejected, --prepared, --pcad, --postfiltered) > 
                        config [input].fits_file > --reduced flag > config [output].reduced_fits
    
    Available metrics:
    - --rmsratio: RMS ratio quality metric
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
            metrics = ['rmsratio', 'squality', 'roll_rms_n', 'mh2o', 'tsys', 'tau_atm', 'chi_sqr', 'err_pwv', 'rms_baseline']
        else:
            # Build list from individual flags
            if rmsratio:
                metrics.append('rmsratio')
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
            click.echo(click.style("Error: No metrics selected. Use --all or at least one of: --rmsratio, --squality, --roll-rms-n, --mh2o, --tsys, --tau-atm, --chi-sqr, --err-pwv, --rms-baseline", fg="red"), err=True)
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
def rmsratio_cmd(config, reduced, clean, rejected, prepared, pcad, object, bins, plot):
    """
    Analyze RMSRATIO quality metric and generate histogram.
    
    Reads input file from config.toml [input].fits_file by default, or uses one of:
    - [output].reduced_fits if --reduced flag is specified
    - [output].clean_fits if --clean flag is specified
    - [output].rejected_fits if --rejected flag is specified
    - [output].prepared_for_pca if --prepared flag is specified
    - [output].pcad_fits if --pcad flag is specified
    
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
        if reduced or clean or rejected or prepared or pcad:
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
              help='Column name for per-spectrum weighting (e.g., RMSRATIO). If specified, each integrated spectrum is weighted by 1/weight_value (lower values get higher weight).')
@click.option('--scatter', is_flag=True, default=False,
              help='Show observation points as scatter plot on map.')
@click.option('--plot', type=click.Path(), default=None,
              help='Output path for plot (e.g., integrated_map.png). If not specified, plot is shown but not saved.')
@click.option('--fits-output', type=click.Path(), default=None,
              help='Output path for FITS file (e.g., integrated_map.fits). If not specified, FITS is not saved.')
def map_integrated_cmd(config, fits, reduced, clean, pcad, rejected, postfiltered, prepared, object, beamsize, pixsize, velocity_range, weight_column, scatter, plot, fits_output):
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

        # Create integrated map
        grid_map, wcs_header, fig = create_integrated_map(
            hdul,
            object_filter=object_filter,
            beamsize_deg=beamsize_deg,
            pixsize=pixsize_deg,
            show_scatter=scatter,
            velocity_range=velocity_range,
            weight_column=weight_column
        )

        click.echo(f"\n✓ Integrated intensity map created")
        click.echo(f"  Beam size: {beamsize_deg*3600:.2f}\"")
        if pixsize_deg:
            click.echo(f"  Pixel size: {pixsize_deg*3600:.3f}\"")
        if object_filter:
            click.echo(f"  Object filter: {object_filter}")
        if weight_column:
            click.echo(f"  Per-spectrum weighting: {weight_column} (using 1/value)")
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
@click.option('--object', type=str, default=None,
              help='Filter by object name (partial match). If not specified, uses "object" from config.toml if available.')
@click.option('--beamsize', type=float, default=None,
              help='Beam size in degrees for gridding kernel. If not specified, reads from config [gridding].beamsize_arcsec.')
@click.option('--pixsize', type=float, default=None,
              help='Map pixel size in degrees. If not specified, uses beamsize/3 or config [gridding].pixel_size_arcsec.')
@click.option('--velocity-range', type=(float, float), default=None, nargs=2,
              help='Velocity range in km/s (e.g., --velocity-range 450 500). If not specified, integrates entire spectrum.')
@click.option('--weight-column', type=str, default=None,
              help='Column name for per-spectrum weighting (e.g., RMSRATIO).')
@click.option('--scatter', is_flag=True, default=False,
              help='Show observation points as scatter plot on map.')
@click.option('--plot', type=click.Path(), default=None,
              help='Output path for plot (e.g., compare_map.png). If not specified, plot is shown but not saved.')
def compare_map_integrated_cmd(config, object, beamsize, pixsize, velocity_range, weight_column, scatter, plot):
    """
    Compare integrated intensity maps from reduced, prepared, and PCA-corrected datasets.

    Creates a three-panel figure showing the integrated map from:
      - Panel 1: reduced data   (output.reduced_fits from config)
      - Panel 2: prepared data  (output.prepared_for_pca from config)
      - Panel 3: PCA-corrected  (output.pcad_fits from config)

    Examples:
        compare_map_integrated --config config.toml
        compare_map_integrated --config config.toml --object M51 --velocity-range 450 500
        compare_map_integrated --config config.toml --weight-column RMSRATIO --plot compare.png
    """
    try:
        from oi_zeigt.mapping.gridding import (
            get_gridding_params_from_config, format_ra_hms, format_dec_dms, HAS_CYGRID
        )
        from matplotlib.ticker import FuncFormatter, MaxNLocator

        config_path = config if config else None
        cfg = get_config(config_path) if config_path else {}

        beamsize_deg, pixsize_deg = get_gridding_params_from_config(
            config_path=config_path,
            beamsize_deg=beamsize,
            pixsize_deg=pixsize
        )

        output_cfg = cfg.get('output', {})
        reduced_path = output_cfg.get('reduced_fits')
        prepared_path = output_cfg.get('prepared_for_pca')
        pcad_path = output_cfg.get('pcad_fits')

        missing = []
        if not reduced_path:
            missing.append('output.reduced_fits')
        if not prepared_path:
            missing.append('output.prepared_for_pca')
        if not pcad_path:
            missing.append('output.pcad_fits')
        if missing:
            click.echo(click.style(f"Error: missing config keys: {', '.join(missing)}", fg='red'), err=True)
            sys.exit(1)

        if object is None:
            object_filter = cfg.get('parameters', {}).get('object', None)
        else:
            object_filter = object

        datasets = [
            ('Reduced', reduced_path),
            ('Prepared', prepared_path),
            ('PCA corrected', pcad_path),
        ]

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
                weight_column=weight_column,
            )
            plt.close(_fig)
            hdul.close()
            maps.append((label, grid_map, wcs_header))

        # Build 3-panel comparison figure
        fig, axes = plt.subplots(1, 3, figsize=(18, 6))

        for ax, (label, grid_map, wcs_header) in zip(axes, maps):
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
        weight_info = f', weighted by {weight_column}' if weight_column else ''
        suptitle = f'Integrated Intensity Comparison ({gridding_method}, beam={beamsize_deg*3600:.1f}″{weight_info})'
        if object_filter:
            suptitle += f' — {object_filter}'
        if velocity_range:
            suptitle += f' [{velocity_range[0]:.0f}–{velocity_range[1]:.0f} km/s]'
        fig.suptitle(suptitle, fontsize=13, fontweight='bold')
        fig.tight_layout()

        click.echo(f"\n✓ Comparison maps created")
        click.echo(f"  Beam size: {beamsize_deg*3600:.2f}\"")
        if object_filter:
            click.echo(f"  Object filter: {object_filter}")
        if weight_column:
            click.echo(f"  Per-spectrum weighting: {weight_column}")
        if velocity_range:
            click.echo(f"  Velocity range: {velocity_range[0]:.1f} - {velocity_range[1]:.1f} km/s")

        if plot:
            plt.savefig(plot, dpi=150)
            click.echo(f"  Plot saved: {plot}")

        click.echo()
        plt.show()
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
@click.option('--reduced', is_flag=True, default=False,
              help='Use reduced_data.fits (output.reduced_fits from config).')
@click.option('--pcad', is_flag=True, default=False,
              help='Use pca_corrected.fits (output.pcad_fits from config).')
@click.option('--prepared', is_flag=True, default=False,
              help='Use prepared_for_pca.fits (output.prepared_for_pca from config).')
@click.option('--object', type=str, default=None,
              help='Filter by object name. If not specified, uses "object" from config.toml if available.')
@click.option('--beamsize', type=float, default=None,
              help='Beam size in degrees for gridding kernel. If not specified, reads from config [gridding].beamsize_arcsec.')
@click.option('--pixsize', type=float, default=None,
              help='Map pixel size in degrees. If not specified, uses beamsize/3 or config [gridding].pixel_size_arcsec.')
@click.option('--output', type=click.Path(), default=None,
              help='Output FITS file path. If not specified, uses "datacube" from config.toml or ./datacube.fits.')
@click.option('--plot', type=click.Path(), default=None,
              help='Output path for diagnostic plot (e.g., datacube_slices.png). If not specified, plot is shown but not saved.')
@click.option('--n-jobs', type=int, default=-1,
              help='Number of parallel workers for channel gridding. -1 = all CPUs (default), 1 = sequential.')
def create_datacube_cmd(config, reduced, pcad, prepared, object, beamsize, pixsize, output, plot, n_jobs):
    """
    Create a full 3D spectral datacube by gridding spectra across spatial and spectral axes.
    
    Grids each velocity channel separately onto a spatial map, creating a proper 3D datacube
    (nvel, dec, ra) with WCS headers. Uses cygrid for optimal Gaussian kernel gridding
    with scipy fallback.
    
    Reads input file from (in order of priority):
    1. --pcad flag (uses output.pcad_fits from config)
    2. --prepared flag (uses output.prepared_for_pca from config)
    3. --reduced flag (uses output.reduced_fits from config)
    4. [input].fits_file from config.toml otherwise

    Output file path is determined by (in order of priority):
    1. --output command line option
    2. [output].datacube from config.toml
    3. Default: ./datacube.fits

    Gridding parameters read from config.toml [gridding] section if not specified on command line.

    Examples:
        create_datacube --config config.toml
        create_datacube --config config.toml --pcad --object M51
        create_datacube --config config.toml --prepared --object M51
        create_datacube --config config.toml --reduced --object M51
        create_datacube --config config.toml --output my_datacube.fits --plot slices.png
        create_datacube --config config.toml --beamsize 0.3 --pixsize 0.05
    """
    try:
        from oi_zeigt.mapping.gridding import get_gridding_params_from_config
        
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
        output_cfg = cfg.get('output', {})
        if pcad:
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
        telescop = gridding_cfg.get('telescop', '')

        datacube, wcs_header, fig = create_spectral_datacube(
            hdul,
            beamsize_deg=beamsize_deg,
            pixsize=pixsize_deg,
            object_filter=object_filter,
            output_file=output_file,
            telescop=telescop,
            n_jobs=n_jobs,
        )

        click.echo(f"\n✓ Spectral datacube created")
        click.echo(f"  Beam size: {beamsize_deg:.4f}°")
        if pixsize_deg:
            click.echo(f"  Pixel size: {pixsize_deg:.6f}°")
        if object_filter:
            click.echo(f"  Object filter: {object_filter}")
        click.echo(f"  Datacube shape: {datacube.shape[0]} channels × {datacube.shape[1]} × {datacube.shape[2]} pixels")
        click.echo(f"  Output file: {output_file}")

        # Save plot if requested
        if plot:
            plt.savefig(plot, dpi=150)
            click.echo(f"  Plot saved: {plot}")

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
    help="Input FITS file to prepare (default: output.reduced_fits from config.toml)"
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
def prepare_for_pca(config: Optional[str], fits: Optional[str], output: Optional[str],
                   pca_source: Optional[str], object: Optional[str], mission_id: Optional[str],
                   scan: Optional[int]):
    """
    Prepare FITS data for PCA analysis.
    
    Filters data to include only specified PCA source and object,
    then fills telluric line regions with Gaussian noise.
    
    Uses configuration from config.toml by default:
    - Input: output.reduced_fits
    - Output: output.prepared_for_pca
    - PCA source: pca.pca_source
    - Object filter: parameters.object
    
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
            scan=scan
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


if __name__ == "__main__":
    print_fits_info()
