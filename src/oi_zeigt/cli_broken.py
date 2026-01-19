"""
Command-line interface for OI ZEIGT.
"""

import sys
from pathlib import Path
from typing import Optional

import click
import numpy as np
import matplotlib.pyplot as plt

from .basic_io import read_fits_from_config, read_fits, get_config
from .reduction.core import (analyze_spectrum_values, detect_blank_channels, 
                            detect_nan_channels, filter_and_save_fits,
                            apply_baseline_from_config, reduce_spectra_from_config)


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
    c        click.style(f"Unexpected error: {e}", fg="red"), err=True)
        import traceback
        traceback.print_exc()
        sys.exit(1)


@click.command()
@click.option('--config', type=click.Path(exists=True), 
              help='Path to config.toml file.')
@click.option('--fits', type=click.Path(exists=True), 
              help='Path to FITS file (overrides config).')
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
def reduce_spectra(config, fits, baseline, baseline_order, baseline_window, 
                   smooth, smooth_window, output):
    """
    Perform spectral reduction with selected methods.
    
    Applies reduction methods in sequence:
    1. Baseline subtraction (--baseline)
    2. Smoothing (--smooth)
    
    More methods can be added in the future.
    
    If output path is not specified:
    - Uses config [output].reduced_fits if available
    - Falls back to 'reduced_data.fits' with a warning
    
    Examples:
        reduce_spectra --config config.toml --baseline
        reduce_spectra --config config.toml --baseline --baseline-order 2 --baseline-window 100 120
        reduce_spectra --fits clean_data.fits --baseline --smooth --output my_reduced.fits
    """
    try:
        # Load config
        config_path = config if config else None
        cfg = get_config(config_path) if config_path else {}
        
        # Determine FITS file to process
        if fits:
            hdul = read_fits(fits)
        else:
            hdul = None  # Let reduce_spectra_from_config handle reading
        
        # Build methods dictionary based on flags
        methods = {}
        
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
                            baseline_window = (int(window_cfg[0]), int(window_cfg[1]))
                    except (ValueError, TypeError, IndexError):
                        baseline_window = None
            
            methods['baseline'] = {
                'order': baseline_order,
                'window': baseline_window
            }
        
        if smooth:
            methods['smooth'] = {
                'window_size': smooth_window
            }
        
        if not methods:
            click.echo(click.style("No reduction methods selected. Use --baseline and/or --smooth.", 
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
        if baseline:
            click.echo(f"    - Baseline subtraction (order={baseline_order}", end='')
            if baseline_window:
                click.echo(f", window={baseline_window[0]}-{baseline_window[1]})")
            else:
                click.echo(")")
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
        sys.exit(1)o(f"Number of HDUs: {len(hdul)}\n")
    
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
        
        # Print details
        _print_fits_details(hdul)
        hdul.close()
        
        click.echo(click.style("✓ Done", fg="green"))
        
    except FileNotFoundError as e:
        click.echo(click.style(f"Error: {e}", fg="red"), err=True)
        sys.exit(1)
    except KeyError as e:
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
    "--object",
    type=str,
    default=None,
    help="Object name to filter (searches within OBJECT field)"
)
@click.option(
    "--output",
    type=click.Path(),
    default=None,
    help="Output file to save the plot (e.g., plot.png)"
)
@click.option(
    "--num-spectra",
    type=int,
    default=20,
    help="Number of spectra to plot (default: 20)"
)
def plot_sample_spectra(config: Optional[str], fits: Optional[str], 
                       object: Optional[str], output: Optional[str], 
                       num_spectra: int):
    """
    Plot a sample of spectra from the FITS file.
    
    Reads spectra from the FITS file and plots them in separate subplots.
    Can filter by object name if specified in config or as argument.
    
    The spectrum index shown in each title can be used to access the spectrum
    in the data array: data[spectrum_index]['SPECTRUM']
    
    Examples:
    
        # Plot 20 spectra using config file
        plot_sample_spectra --config config.toml
        
        # Plot 15 spectra of a specific object
        plot_sample_spectra --config config.toml --object "M51" --num-spectra 15
        
        # Save plot to file
        plot_sample_spectra --config config.toml --output spectra.png
    """
    try:
        # Read config and FITS file
        if fits:
            click.echo(f"Reading FITS file: {fits}")
            hdul = read_fits(fits)
        elif config:
            click.echo(f"Reading config from: {config}")
            hdul = read_fits_from_config(config)
            config_data = get_config(config)
        else:
            click.echo("Reading config from default location...")
            hdul = read_fits_from_config()
            config_data = get_config()
        
        # Get object name from config if not specified as argument
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
            click.echo(click.style("Warning: No object specified. Will plot any spectra found.", fg="yellow"))
        
        # Find the binary table HDU with spectra
        matrix_hdu = None
        for hdu in hdul:
            if hasattr(hdu, 'data') and hdu.data is not None:
                if 'SPECTRUM' in hdu.data.dtype.names:
                    matrix_hdu = hdu
                    break
        
        if matrix_hdu is None:
            raise ValueError("No HDU with SPECTRUM column found")
        
        # Filter by object if specified
        data = matrix_hdu.data
        if object:
            mask = np.array([object.lower() in str(obj).lower() 
                           for obj in data['OBJECT']])
            filtered_indices = np.where(mask)[0]
            filtered_data = data[mask]
            click.echo(f"Found {len(filtered_data)} entries matching object '{object}'")
        else:
            filtered_indices = np.arange(len(data))
            filtered_data = data
            click.echo(f"Using all {len(filtered_data)} entries")
        
        if len(filtered_data) == 0:
            raise ValueError(f"No spectra found matching object '{object}'")
        
        # Limit to num_spectra
        num_to_plot = min(num_spectra, len(filtered_data))
        sample_positions = np.random.choice(len(filtered_data), size=num_to_plot, replace=False)
        
        # Get both the original indices and the data
        original_indices = filtered_indices[sample_positions]
        spectra_to_plot = filtered_data[sample_positions]
        
        click.echo(f"Plotting {num_to_plot} sample spectra...")
        click.echo(f"Original indices in FITS data: {sorted(original_indices)}\n")
        
        # Create figure with subplots (4 rows x 5 columns = 20 plots)
        rows = int(np.ceil(num_to_plot / 5))
        cols = min(5, num_to_plot)
        
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
    It samples spectra and shows the most common values, which often reveals
    the blank value marker (commonly 0.0, NaN, or a specific negative number).
    
    Examples:
    
        # Analyze using config file
        analyze_blanks --config config.toml
        
        # Analyze direct file
        analyze_blanks --fits /path/to/file.fits
        
        # Sample more spectra for better statistics
        analyze_blanks --config config.toml --sample-size 500
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
        
        click.echo(f"\nAnalyzing {sample_size} spectra...\n")
        
        # Analyze
        analysis = analyze_spectrum_values(hdul, sample_size=sample_size)
        
        # Print results
        click.echo("="*70)
        click.echo("SPECTRUM VALUE ANALYSIS")
        click.echo("="*70 + "\n")
        
        click.echo(f"Total values analyzed: {analysis['total_values_analyzed']:,}")
        click.echo(f"Value range: [{analysis['min_value']:.2e}, {analysis['max_value']:.2e}]\n")
        
        click.echo("Statistics:")
        click.echo(f"  Mean:     {analysis['mean_value']:12.4e}")
        click.echo(f"  Median:   {analysis['median_value']:12.4e}")
        click.echo(f"  Std Dev:  {analysis['std_value']:12.4e}\n")
        
        click.echo("Special value counts:")
        click.echo(f"  Zeros (0.0):           {analysis['zero_count']:10,}")
        click.echo(f"  Negative values:       {analysis['negative_count']:10,}")
        click.echo(f"  Very small (< 1e-10):  {analysis['very_small_count']:10,}\n")
        
        click.echo("="*70)
        click.echo("TOP 20 MOST COMMON VALUES")
        click.echo("="*70)
        click.echo(f"{'Value':<20} {'Count':<15} {'Percentage':<10}")
        click.echo("-"*70)
        
        for value, count in analysis['common_values']:
            percentage = (count / analysis['total_values_analyzed']) * 100
            click.echo(f"{value:<20.4e} {count:<15,} {percentage:>8.2f}%")
        
        click.echo("="*70 + "\n")
        
        click.echo(click.style(
            "💡 Tip: The most common value is likely your blank/missing marker.",
            fg="blue"
        ))
        click.echo(click.style(
            "   Use this value with split_spectra_by_blanks(blank_value=...)",
            fg="blue"
        ))
        click.echo()
        
        hdul.close()
        click.echo(click.style("✓ Done", fg="green"))
        
    except FileNotFoundError as e:
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
    "--skydiff",
    is_flag=True,
    default=False,
    help="Plot SKY-DIFF spectra instead of SKYCHOPDIFF (default: SKYCHOPDIFF)"
)
@click.option(
    "--output",
    type=click.Path(),
    default=None,
    help="Output file to save the plot (e.g., plot.png)"
)
@click.option(
    "--num-spectra",
    type=int,
    default=20,
    help="Number of spectra to plot (default: 20)"
)
def plot_skies(config: Optional[str], fits: Optional[str], skydiff: bool,
               output: Optional[str], num_spectra: int):
    """
    Plot sample sky observation spectra (SKYCHOPDIFF or SKY-DIFF).
    
    By default plots SKYCHOPDIFF spectra. Use --skydiff flag to plot SKY-DIFF instead.
    
    Examples:
    
        # Plot 20 SKYCHOPDIFF spectra
        plot_skies --config config.toml
        
        # Plot 20 SKY-DIFF spectra
        plot_skies --config config.toml --skydiff
        
        # Plot 15 SKY-DIFF spectra and save
        plot_skies --config config.toml --skydiff --num-spectra 15 --output skies.png
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
        
        # Determine sky type to plot
        sky_type = "SKY-DIFF" if skydiff else "SKYCHOPDIFF"
        click.echo(f"Looking for {sky_type} spectra...\n")
        
        # Find the binary table HDU with spectra
        matrix_hdu = None
        for hdu in hdul:
            if hasattr(hdu, 'data') and hdu.data is not None:
                if 'SPECTRUM' in hdu.data.dtype.names:
                    matrix_hdu = hdu
                    break
        
        if matrix_hdu is None:
            raise ValueError("No HDU with SPECTRUM column found")
        
        # Filter by sky type
        data = matrix_hdu.data
        
        # Check if OBJECT column exists and filter accordingly
        if 'OBJECT' in data.dtype.names:
            mask = np.array([sky_type.lower() in str(obj).lower() 
                           for obj in data['OBJECT']])
            filtered_indices = np.where(mask)[0]
            filtered_data = data[mask]
        else:
            # If no OBJECT column, use all data
            click.echo(click.style("Warning: No OBJECT column found. Using all spectra.", fg="yellow"))
            filtered_indices = np.arange(len(data))
            filtered_data = data
        
        if len(filtered_data) == 0:
            raise ValueError(f"No spectra found with OBSMODE matching '{sky_type}'")
        
        click.echo(f"Found {len(filtered_data)} {sky_type} entries")
        
        # Limit to num_spectra
        num_to_plot = min(num_spectra, len(filtered_data))
        sample_positions = np.random.choice(len(filtered_data), size=num_to_plot, replace=False)
        
        # Get both the original indices and the data
        original_indices = filtered_indices[sample_positions]
        spectra_to_plot = filtered_data[sample_positions]
        
        click.echo(f"Plotting {num_to_plot} sample {sky_type} spectra...")
        click.echo(f"Original indices in FITS data: {sorted(original_indices)}\n")
        
        # Create figure with subplots (4 rows x 5 columns = 20 plots)
        rows = int(np.ceil(num_to_plot / 5))
        cols = min(5, num_to_plot)
        
        fig, axes = plt.subplots(rows, cols, figsize=(15, 3*rows))
        
        # Flatten axes array for easier iteration
        if rows == 1 and cols == 1:
            axes = np.array([[axes]])
        elif rows == 1 or cols == 1:
            axes = axes.reshape(rows, cols)
        
        # Plot each spectrum
        for plot_num, (ax, spectrum_data, orig_idx) in enumerate(zip(axes.flat, spectra_to_plot, original_indices)):
            spectrum = spectrum_data['SPECTRUM']
            
            # Try to get object name if available
            obj_name = "Unknown"
            if 'OBJECT' in data.dtype.names:
                obj_data = spectrum_data['OBJECT']
                obj_name = obj_data.decode().strip() if isinstance(obj_data, bytes) else str(obj_data).strip()
            
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
        
        fig.suptitle(f"{sky_type} Spectra", fontsize=14, weight='bold', y=0.98)
        
        plt.tight_layout(rect=[0, 0.04, 1, 0.97])
        
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
    "--object",
    type=str,
    default=None,
    help="Object name to filter (e.g., M51). If not specified, reads from config."
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
def filter_fits(config: Optional[str], fits: Optional[str], object: Optional[str],
                nan_threshold: float, output_clean: Optional[str], 
                output_rejected: Optional[str]):
    """
    Filter FITS data by object and NaN content.
    
    Creates two FITS files:
    1. Clean file: All non-target objects + target object spectra with < nan_threshold NaNs
    2. Rejected file: Target object spectra with >= nan_threshold NaNs
    
    Examples:
    
        # Filter using config file (reads object name from config)
        filter_fits --config config.toml
        
        # Filter with custom thresholds
        filter_fits --config config.toml --nan-threshold 0.75
        
        # Filter specific object and specify output files
        filter_fits --config config.toml --object "M51" \\
            --output-clean m51_clean.fits --output-rejected m51_rejected.fits
    """
    try:
        # Read FITS file
        if fits:
            click.echo(f"Reading FITS file: {fits}")
            hdul = read_fits(fits)
        elif config:
            click.echo(f"Reading config from: {config}")
            hdul = read_fits_from_config(config)
            config_data = get_config(config)
        else:
            click.echo("Reading config from default location...")
            hdul = read_fits_from_config()
            config_data = get_config()
        
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
        click.echo(f"Keeping spectra with < {nan_threshold:.1%} NaN channels")
        click.echo(f"Rejecting spectra with >= {nan_threshold:.1%} NaN channels\n")
        
        # Filter and save
        clean_path, rejected_path = filter_and_save_fits(
            hdul,
            object_name=object,
            nan_threshold=nan_threshold,
            output_clean=output_clean,
            output_rejected=output_rejected
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
              help='Channel range [start end] to ignore during baseline fitting. '
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
def reduce_spectra_cmd(config, fits, baseline, baseline_order, baseline_window, 
                       smooth, smooth_window, output):
    """
    Perform spectral reduction with selected methods.
    
    Applies reduction methods in sequence:
    1. Baseline subtraction (--baseline)
    2. Smoothing (--smooth)
    
    More methods can be added in the future.
    
    If output path is not specified:
    - Uses config [output].reduced_fits if available
    - Falls back to 'reduced_data.fits' with a warning
    
    Examples:
        reduce_spectra --config config.toml --baseline
        reduce_spectra --config config.toml --baseline --baseline-order 2 --baseline-window 100 120
        reduce_spectra --fits clean_data.fits --baseline --smooth --output my_reduced.fits
    """
    try:
        # Load config
        config_path = config if config else None
        cfg = get_config(config_path) if config_path else {}
        
        # Determine FITS file to process
        if fits:
            hdul = read_fits(fits)
        else:
            hdul = None  # Let reduce_spectra_from_config handle reading
        
        # Build methods dictionary based on flags
        methods = {}
        
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
                            baseline_window = (int(window_cfg[0]), int(window_cfg[1]))
                    except (ValueError, TypeError, IndexError):
                        baseline_window = None
            
            methods['baseline'] = {
                'order': baseline_order,
                'window': baseline_window
            }
        
        if smooth:
            methods['smooth'] = {
                'window_size': smooth_window
            }
        
        if not methods:
            click.echo(click.style("No reduction methods selected. Use --baseline and/or --smooth.", 
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
        if baseline:
            click.echo(f"    - Baseline subtraction (order={baseline_order}", end='')
            if baseline_window:
                click.echo(f", window={baseline_window[0]}-{baseline_window[1]})")
            else:
                click.echo(")")
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

