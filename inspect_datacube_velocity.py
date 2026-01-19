#!/usr/bin/env python3
"""
Practical example: Verify velocity axis in the output datacube

This script demonstrates how to:
1. Load a datacube created with the modified create_spectral_datacube()
2. Extract and inspect the velocity axis parameters
3. Convert pixel coordinates to velocity
4. Select velocity ranges
"""

from pathlib import Path
from astropy.io import fits
from astropy.wcs import WCS
import numpy as np

def inspect_datacube_velocity_axis(fits_file):
    """
    Inspect the velocity axis of a spectral datacube FITS file.
    
    Parameters
    ----------
    fits_file : str or Path
        Path to the datacube FITS file
    """
    fits_file = Path(fits_file)
    
    if not fits_file.exists():
        print(f"Error: File not found: {fits_file}")
        return
    
    print(f"Inspecting datacube: {fits_file}")
    print()
    
    # Open the FITS file
    with fits.open(fits_file) as hdul:
        header = hdul[0].header
        data = hdul[0].data
        
        print("=" * 80)
        print("DATACUBE STRUCTURE")
        print("=" * 80)
        print(f"Data shape: {data.shape}")
        print(f"  Channels (NAXIS3): {data.shape[0]}")
        print(f"  Dec pixels (NAXIS2): {data.shape[1]}")
        print(f"  RA pixels (NAXIS1): {data.shape[2]}")
        print()
        
        # Extract velocity axis parameters
        print("=" * 80)
        print("VELOCITY AXIS PARAMETERS")
        print("=" * 80)
        
        # WCS information
        crval3 = header.get('CRVAL3', None)
        crpix3 = header.get('CRPIX3', None)
        cdelt3 = header.get('CDELT3', None)
        ctype3 = header.get('CTYPE3', None)
        cunit3 = header.get('CUNIT3', None)
        
        print(f"CTYPE3 (axis type): {ctype3}")
        print(f"CUNIT3 (units): {cunit3}")
        print(f"CRVAL3 (reference value): {crval3:.2f}")
        print(f"CRPIX3 (reference pixel): {crpix3:.1f}")
        print(f"CDELT3 (pixel scale): {cdelt3:.2f}")
        print()
        
        # Calculate velocity for each channel
        nchans = header.get('NAXIS3', 0)
        print("=" * 80)
        print("VELOCITY RANGE")
        print("=" * 80)
        
        # Using FITS WCS convention (1-indexed pixels)
        pixels = np.arange(nchans)
        velocities = crval3 + (pixels - (crpix3 - 1)) * cdelt3
        
        print(f"Number of channels: {nchans}")
        print(f"First channel (pixel 0): {velocities[0]:.2f} m/s = {velocities[0]/1000:.2f} km/s")
        print(f"Last channel (pixel {nchans-1}): {velocities[-1]:.2f} m/s = {velocities[-1]/1000:.2f} km/s")
        print(f"Total velocity span: {(velocities[-1] - velocities[0]):.2f} m/s = {(velocities[-1] - velocities[0])/1000:.2f} km/s")
        print()
        
        # Find channels within a velocity range
        print("=" * 80)
        print("EXAMPLE: SELECT VELOCITY RANGE")
        print("=" * 80)
        
        v_min = 400000.0  # 400 km/s
        v_max = 500000.0  # 500 km/s
        
        mask = (velocities >= v_min) & (velocities <= v_max)
        channels = np.where(mask)[0]
        
        print(f"Velocity range: {v_min/1000:.2f} to {v_max/1000:.2f} km/s")
        print(f"Channels in this range: {len(channels)} channels")
        if len(channels) > 0:
            print(f"  First channel: {channels[0]} (velocity {velocities[channels[0]]/1000:.2f} km/s)")
            print(f"  Last channel: {channels[-1]} (velocity {velocities[channels[-1]]/1000:.2f} km/s)")
        print()
        
        # Beam information
        print("=" * 80)
        print("BEAM INFORMATION")
        print("=" * 80)
        bmaj = header.get('BMAJ', None)
        bmin = header.get('BMIN', None)
        bpa = header.get('BPA', None)
        
        if bmaj is not None:
            print(f"BMAJ (beam major axis): {bmaj:.4f}°")
        if bmin is not None:
            print(f"BMIN (beam minor axis): {bmin:.4f}°")
        if bpa is not None:
            print(f"BPA (beam position angle): {bpa:.1f}°")
        print()
        
        # Data statistics
        print("=" * 80)
        print("DATA STATISTICS")
        print("=" * 80)
        print(f"Data type: {data.dtype}")
        print(f"Min value: {np.nanmin(data):.6f}")
        print(f"Max value: {np.nanmax(data):.6f}")
        print(f"Mean value: {np.nanmean(data):.6f}")
        print(f"Median value: {np.nanmedian(data):.6f}")
        print(f"Std deviation: {np.nanstd(data):.6f}")
        print(f"NaN pixels: {np.sum(~np.isfinite(data))}")
        print(f"Valid pixels: {np.sum(np.isfinite(data))}")
        print()
        
        # Channel statistics
        print("=" * 80)
        print("CHANNEL ANALYSIS")
        print("=" * 80)
        
        # Find channel with most data
        valid_per_channel = np.sum(np.isfinite(data), axis=(1, 2))
        max_channel = np.argmax(valid_per_channel)
        mean_per_channel = np.nanmean(data, axis=(1, 2))
        
        print(f"Channel with most valid pixels: {max_channel} ({valid_per_channel[max_channel]} pixels)")
        print(f"  Velocity: {velocities[max_channel]/1000:.2f} km/s")
        print()
        
        # Cube statistics
        peak_channel = np.nanargmax(mean_per_channel)
        print(f"Channel with highest mean intensity: {peak_channel}")
        print(f"  Velocity: {velocities[peak_channel]/1000:.2f} km/s")
        print(f"  Mean intensity: {mean_per_channel[peak_channel]:.6f}")


if __name__ == "__main__":
    # Example usage
    print("DATACUBE VELOCITY AXIS INSPECTION EXAMPLE")
    print()
    
    # You would replace this with your actual datacube file
    datacube_file = Path("/home/diskB/sofiaobsdata/m51/m51_central_fits/fits_files_processed/datacube.fits")
    
    if datacube_file.exists():
        inspect_datacube_velocity_axis(datacube_file)
    else:
        print(f"Note: Example datacube not found at {datacube_file}")
        print()
        print("To use this script:")
        print("  1. Create a datacube using: create_datacube --config config.toml")
        print("  2. Run this script with the path to your datacube FITS file")
        print()
        print("Example:")
        print("  python3 inspect_datacube_velocity.py /path/to/datacube.fits")
