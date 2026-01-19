#!/usr/bin/env python3
"""
Test script to verify that mapping functions work with multi-HDU FITS files.
"""

from pathlib import Path
from astropy.io import fits
import numpy as np
import sys

# Add src to path for imports
sys.path.insert(0, str(Path(__file__).parent / 'src'))

from oi_zeigt.mapping.gridding import _get_celestial_coords, create_integrated_map

# Test on the combined multi-HDU file
fits_file = Path('/home/diskB/sofiaobsdata/m51/m51_central_region/kalibrated_all/m51_central_tile.fits')

if not fits_file.exists():
    print(f"Error: Test file not found: {fits_file}")
    sys.exit(1)

print(f"Testing mapping on multi-HDU combined file: {fits_file}")
print(f"File size: {fits_file.stat().st_size / (1024**3):.2f} GB\n")

# Open the FITS file
with fits.open(fits_file) as hdul:
    # Print file structure
    print("FITS file structure:")
    hdul.info()
    print()
    
    # Count rows in each HDU
    total_rows = 0
    for i, hdu in enumerate(hdul):
        if hasattr(hdu, 'data') and hdu.data is not None and 'SPECTRUM' in hdu.data.dtype.names:
            print(f"HDU {i} ({hdu.name}): {len(hdu.data)} rows")
            total_rows += len(hdu.data)
    
    print(f"Total rows across all HDUs: {total_rows}\n")
    
    # Test _get_celestial_coords with multi-HDU file
    print("Testing _get_celestial_coords()...")
    try:
        ras, decs = _get_celestial_coords(hdul)
        print(f"✓ Successfully retrieved coordinates from all HDUs")
        print(f"  Total coordinate points: {len(ras)}")
        print(f"  RA range: {ras.min():.4f}° to {ras.max():.4f}°")
        print(f"  Dec range: {decs.min():.4f}° to {decs.max():.4f}°")
        
        # Verify we got coordinates from ALL HDUs (should be ~48,188 total)
        if len(ras) == total_rows:
            print(f"✓ PASSED: Got all {total_rows} coordinates (all HDUs processed)")
        else:
            print(f"✗ FAILED: Got {len(ras)} coordinates but expected {total_rows} (all HDUs)")
            sys.exit(1)
    except Exception as e:
        print(f"✗ FAILED: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
    
    print()
    
    # Test create_integrated_map with multi-HDU file
    print("Testing create_integrated_map()...")
    try:
        # This will create a map from integrated spectra across all HDUs
        map_data, wcs_header, fig = create_integrated_map(
            hdul,
            spectrum_column='SPECTRUM',
            beamsize_deg=0.25,
            pixsize=None,
            object_filter=None
        )
        
        print(f"✓ Successfully created integrated map from all HDUs")
        print(f"  Map shape (Dec, RA): {map_data.shape}")
        print(f"  Map data range: {np.nanmin(map_data):.2e} to {np.nanmax(map_data):.2e}")
        
        # Close the figure to avoid display issues
        import matplotlib.pyplot as plt
        plt.close(fig)
        
    except Exception as e:
        print(f"✗ FAILED: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

print("\n✓ All tests PASSED!")
print("Mapping functions now correctly process all HDUs in multi-HDU combined files.")
