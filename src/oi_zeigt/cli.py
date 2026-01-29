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
from .mapping.gridding import create_map_from_column, create_integrated_map


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
    hdul.info()
    click.echo("="*70 + "\n")
    
    # Print additional details
    click.echo(f"Number of HDUs: {len(hdul)}\n")
    
    for i, hdu in enumerate(hdul):
        click.echo(f"HDU {i}: {hdu.name} ({type(hdu).__name__})")
        if hdu.data is not None:
            click.echo(f"  Data shape: {hdu.data.shape}")
            click.echo(f"  Data type: {hdu.data.dtype}")
        if hdu.header:
            click.echo(f"  Header keywords: {len(hdu.header)}")
        click.echo()
    
    # Print object information
    _print_object_info(hdul)


def _print_object_info(hdul):
    """
    Print information about unique objects in the FITS file.
    
    Parameters
    ----------
    hdul : astropy.io.fits.HDUList
        The FITS HDU list to analyze.
    """
    # Look for OBJECT column in binary tables
    for hdu in hdul:
        if hasattr(hdu, 'data') and hdu.data is not None:
            # Check if 'OBJECT' column exists
            if 'OBJECT' in hdu.data.dtype.names:
                objects = hdu.data['OBJECT']
                
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
def print_fits_info(config: Optional[str], fits: Optional[str]):
    """
    Print basic information about a FITS file.
    
    Can read the FITS file path from a config.toml file or directly specify it.
    
    Examples:
    
        # Read FITS file from config.toml
        print_fits_info --config config.toml
        
        # Read FITS file directly
        print_fits_info --fits /path/to/file.fits
        
        # Use default config.toml from current/parent directory
        print_fits_info
    """
    try:
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
    help="Path to FITS file to read directly"
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


def plot_sample_spectra(config: Optional[str], fits: Optional[str], object: Optional[str],
                       num_spectra: int, output: Optional[str]):
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
        # Always try to load config (for object name, thresholds, etc.)
        config_data = {}
        try:
            if config:
                config_data = get_config(config)
                click.echo(f"Loaded config from: {config}")
            else:
                config_data = get_config()
                click.echo("Loaded config from default location")
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
    help="Path to FITS file to read directly"
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
def plot_skies(config: Optional[str], fits: Optional[str], num_spectra: int, 
              output: Optional[str]):
    """
    Plot a sample of sky/background spectra (SKYCHOPDIFF or SKY-DIFF observations).
    
    Examples:
    
        plot_skies --config config.toml
        plot_skies --config config.toml --num-spectra 30
        plot_skies --fits /path/to/file.fits --output sky_plot.pdf
    """
    try:
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
    help="Path to FITS file to read directly"
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
    help="Output path for clean FITS file (default: clean_data.fits)"
)
@click.option(
    "--output-rejected",
    type=click.Path(),
    default=None,
    help="Output path for rejected FITS file (default: rejected_data.fits)"
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
def filter_fits(config: Optional[str], fits: Optional[str], object: Optional[str],
                nan_threshold: float, output_clean: Optional[str], 
                output_rejected: Optional[str], remove: Optional[str],
                remove_values: tuple, apply_only_to_object: bool):
    """
    Filter FITS data by object, NaN content, and/or column values.
    
    Creates two FITS files:
    1. Clean file: Spectra passing NaN threshold filter
    2. Rejected file: Spectra failing NaN threshold filter or matching removal criteria
    
    By default, NaN filtering is applied to ALL spectra. Use --apply-only-to-object
    to filter only the target object and keep all other objects regardless of NaN content.
    
    Examples:
    
        # Filter ALL spectra by NaN threshold (default)
        filter_fits --config config.toml
        
        # Filter only M51 spectra, keep all other objects
        filter_fits --config config.toml --apply-only-to-object
        
        # Filter all with custom thresholds
        filter_fits --config config.toml --nan-threshold 0.75
        
        # Filter specific object and specify output files
        filter_fits --config config.toml --object "M51" \\
            --output-clean m51_clean.fits --output-rejected m51_rejected.fits
        
        # Remove specific AOR_ID values AND filter all spectra
        filter_fits --config config.toml --remove AOR_ID \\
            --remove-values 04_0116_0020609 --remove-values 04_0116_0020506
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
        
        # Check for removal criteria
        if remove and remove_values:
            click.echo(f"Also removing rows where {remove} = {', '.join(remove_values)}\n")
        
        # Filter and save
        clean_path, rejected_path = filter_and_save_fits(
            hdul,
            object_name=object,
            nan_threshold=nan_threshold,
            output_clean=output_clean,
            output_rejected=output_rejected,
            remove_column=remove,
            remove_values=list(remove_values) if remove_values else None,
            apply_to_all=not apply_only_to_object
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
        click.echo(f"  (Contains only {object} spectra with >= {nan_threshold:.1%} NaN channels)")
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
              help='Channel range [start end] to exclude from baseline fitting.')
@click.option('--smooth', is_flag=True, default=False,
              help='Apply smoothing.')
@click.option('--smooth-window', type=int, default=5,
              help='Smoothing window size (default: 5).')
@click.option('--output', type=click.Path(), default=None,
              help='Output FITS file path.')
def reduce_spectra_cmd(config, fits, clean, unblank, baseline, baseline_order, baseline_window, 
                       smooth, smooth_window, output):
    """
    Perform spectral reduction with selected methods.
    
    Applies reduction methods in sequence:
    1. Unblank (--unblank): Fill NaN values using interpolation
    2. Baseline subtraction (--baseline): Remove polynomial baseline
    3. Smoothing (--smooth): Apply boxcar smoothing
    
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
        
        if baseline:
            # Get baseline order
            if baseline_order is None:
                reduction_cfg = cfg.get('reduction', {})
                baseline_order = reduction_cfg.get('baseline_order', 
                                                   reduction_cfg.get('baseline', 1))
                try:
                    baseline_order = int(baseline_order)
                except Exception:
                    baseline_order = 1
            
            # Get baseline window
            if baseline_window is None:
                reduction_cfg = cfg.get('reduction', {})
                window_cfg = reduction_cfg.get('window', None)
                if window_cfg is not None:
                    try:
                        if isinstance(window_cfg, (list, tuple)) and len(window_cfg) == 2:
                            # Window from config is in km/s
                            baseline_window = window_cfg
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
        if unblank:
            click.echo(f"    - Unblank (fill NaN values with linear interpolation)")
        if baseline:
            if baseline_window:
                click.echo(f"    - Baseline subtraction (order={baseline_order}, window={baseline_window[0]}-{baseline_window[1]})")
            else:
                click.echo(f"    - Baseline subtraction (order={baseline_order})")
        if smooth:
            click.echo(f"    - Smoothing (window={smooth_window})")
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
@click.option('--reduced', is_flag=True, default=False,
              help='Use reduced_data.fits instead of the input data.')
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
@click.option('--bins', type=int, default=30,
              help='Number of histogram bins (default: 30).')
@click.option('--plot', type=click.Path(), default=None,
              help='Output path for plot. If not specified, plot is shown but not saved.')
def spechistogram_cmd(config, reduced, object, all_metrics, rmsratio, squality, roll_rms_n, mh2o, 
                      tsys, tau_atm, chi_sqr, err_pwv, bins, plot):
    """
    Generate histograms of multiple spectral quality metrics.
    
    Creates a combined multi-panel figure showing histograms of selected quality metrics.
    Each metric shows mean and median lines.
    
    Reads input file from config.toml [input].fits_file by default, or uses
    [output].reduced_fits if --reduced flag is specified.
    
    Available metrics:
    - --rmsratio: RMS ratio quality metric
    - --squality: Signal quality flag
    - --roll-rms-n: RMS rolloff (mean across channels)
    - --mh2o: Water vapor column density
    - --tsys: System temperature
    - --tau-atm: Atmospheric optical depth
    - --chi-sqr: Chi-square fit value
    - --err-pwv: PWV error
    
    Use --all to include all available metrics.
    
    Examples:
        spechistogram --config config.toml --rmsratio --tsys
        spechistogram --config config.toml --reduced --all
        spechistogram --config config.toml --all --object M51 --plot metrics.png
    """
    try:
        # Load config
        config_path = config if config else None
        cfg = get_config(config_path) if config_path else {}
        
        # Determine FITS file to process
        if reduced:
            # Use reduced_data.fits from config
            output_cfg = cfg.get('output', {})
            reduced_fits_path = output_cfg.get('reduced_fits', None)
            if not reduced_fits_path:
                click.echo(click.style("Error: --reduced flag specified but [output].reduced_fits not defined in config", fg="red"), err=True)
                sys.exit(1)
            hdul = fits.open(reduced_fits_path)
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
        
        # Build list of metrics to plot
        metrics = []
        
        # If --all is specified, enable all metrics
        if all_metrics:
            metrics = ['rmsratio', 'squality', 'roll_rms_n', 'mh2o', 'tsys', 'tau_atm', 'chi_sqr', 'err_pwv']
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
        
        if not metrics:
            click.echo(click.style("Error: No metrics selected. Use --all or at least one of: --rmsratio, --squality, --roll-rms-n, --mh2o, --tsys, --tau-atm, --chi-sqr, --err-pwv", fg="red"), err=True)
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
@click.option('--object', type=str, default=None,
              help='Filter by object name (e.g., "M51"). Defaults to "object" from config.toml if not specified. If not in config either, all spectra are analyzed.')
@click.option('--bins', type=int, default=30,
              help='Number of histogram bins (default: 30).')
@click.option('--plot', type=click.Path(), default=None,
              help='Output path for plot (e.g., rmsratio_hist.png). If not specified, plot is shown but not saved.')
def rmsratio_cmd(config, reduced, object, bins, plot):
    """
    Analyze RMSRATIO quality metric and generate histogram.
    
    Reads input file from config.toml [input].fits_file by default, or uses
    [output].reduced_fits if --reduced flag is specified.
    
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
        rmsratio --config config.toml --object M51
        rmsratio --config config.toml --bins 50 --plot hist.png
    """
    try:
        # Load config
        config_path = config if config else None
        cfg = get_config(config_path) if config_path else {}
        
        # Determine FITS file to process
        if reduced:
            # Use reduced_data.fits from config
            output_cfg = cfg.get('output', {})
            reduced_fits_path = output_cfg.get('reduced_fits', None)
            if not reduced_fits_path:
                click.echo(click.style("Error: --reduced flag specified but [output].reduced_fits not defined in config", fg="red"), err=True)
                sys.exit(1)
            hdul = fits.open(reduced_fits_path)
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
@click.option('--object', type=str, default=None,
              help='Filter by object name (partial match). If not specified, uses "object" from config.toml if available.')
@click.option('--beamsize', type=float, default=None,
              help='Beam size in degrees for gridding kernel. If not specified, reads from config [gridding].beamsize_arcsec.')
@click.option('--pixsize', type=float, default=None,
              help='Map pixel size in degrees. If not specified, uses beamsize/3 or config [gridding].pixel_size_arcsec.')
@click.option('--scatter', is_flag=True, default=False,
              help='Show observation points as scatter plot on map.')
@click.option('--plot', type=click.Path(), default=None,
              help='Output path for plot (e.g., integrated_map.png). If not specified, plot is shown but not saved.')
@click.option('--fits-output', type=click.Path(), default=None,
              help='Output path for FITS file (e.g., integrated_map.fits). If not specified, FITS is not saved.')
def map_integrated_cmd(config, fits, reduced, clean, object, beamsize, pixsize, scatter, plot, fits_output):
    """
    Create a spatial map of integrated spectral intensity.
    
    Integrates the spectrum across all frequency channels for each observation,
    then creates a WCS-based spatial map with cygrid gridding.
    
    Reads input file from:
    1. --fits parameter if specified (overrides config)
    2. [output].reduced_fits if --reduced flag is specified
    3. [input].fits_file from config.toml otherwise
    
    Gridding parameters read from config.toml [gridding] section if not specified on command line.

    Examples:
        map_integrated --config config.toml
        map_integrated --config config.toml --fits /path/to/data.fits
        map_integrated --config config.toml --reduced --object M51
        map_integrated --config config.toml --fits data.fits --beamsize 0.3 --plot integrated.png
        map_integrated --config config.toml --fits data.fits --pixsize 0.05 --object M51
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
        if fits:
            # Use directly specified FITS file (overrides everything)
            click.echo(f"Reading FITS file: {fits}")
            hdul = read_fits(fits)
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
            show_scatter=scatter
        )

        click.echo(f"\n✓ Integrated intensity map created")
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
@click.option('--reduced', is_flag=True, default=False,
              help='Use reduced_data.fits instead of the input data.')
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
def create_datacube_cmd(config, reduced, object, beamsize, pixsize, output, plot):
    """
    Create a full 3D spectral datacube by gridding spectra across spatial and spectral axes.
    
    Grids each velocity channel separately onto a spatial map, creating a proper 3D datacube
    (nvel, dec, ra) with WCS headers. Uses cygrid for optimal Gaussian kernel gridding
    with scipy fallback.
    
    Reads input file from config.toml [input].fits_file by default, or uses
    [output].reduced_fits if --reduced flag is specified.
    
    Output file path is determined by (in order of priority):
    1. --output command line option
    2. [output].datacube from config.toml
    3. Default: ./datacube.fits

    Gridding parameters read from config.toml [gridding] section if not specified on command line.
    
    Examples:
        create_datacube --config config.toml
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
        if reduced:
            # Use reduced_data.fits from config
            output_cfg = cfg.get('output', {})
            reduced_fits_path = output_cfg.get('reduced_fits', None)
            if not reduced_fits_path:
                click.echo(click.style("Error: --reduced flag specified but [output].reduced_fits not defined in config", fg="red"), err=True)
                sys.exit(1)
            hdul = fits.open(reduced_fits_path)
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
        
        datacube, wcs_header, fig = create_spectral_datacube(
            hdul,
            beamsize_deg=beamsize_deg,
            pixsize=pixsize_deg,
            object_filter=object_filter,
            output_file=output_file
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


if __name__ == "__main__":
    print_fits_info()
