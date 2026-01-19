#!/usr/bin/env python3
"""
Test script for map_integrated functionality.

Tests the create_spectral_cube function with the available test data.
"""

from pathlib import Path
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt

# Add src to path for imports
import sys
sys.path.insert(0, str(Path(__file__).parent / 'src'))

from oi_zeigt.mapping.gridding import create_spectral_cube

# Test data
fits_file = Path(__file__).parent / 'averaged_data.fits'

if not fits_file.exists():
    print(f"Error: Test file not found: {fits_file}")
    sys.exit(1)

print(f"Loading test FITS file: {fits_file}")
print(f"File size: {fits_file.stat().st_size / 1024:.1f} KB\n")

# Open the FITS file
hdul = fits.open(fits_file)

# Print file info
print("FITS file info:")
hdul.info()
print()

# Try to create spectral cube
try:
    print("Attempting to create spectral cube...")
    cube_data, wcs_header, header = create_spectral_cube(
        hdul,
        object_filter=None,
        grid_method='linear',
        pixsize=0.01,
        beamsize_deg=0.05,
    )
    
    print(f"✓ Spectral cube created successfully!")
    print(f"  Cube shape (channels, DEC, RA): {cube_data.shape}")
    print(f"  Data type: {cube_data.dtype}")
    print(f"  Data range: [{np.nanmin(cube_data):.2f}, {np.nanmax(cube_data):.2f}]")
    print(f"  Number of valid pixels: {np.sum(np.isfinite(cube_data))}")
    print()
    
    # Try to plot integrated intensity
    print("Creating integrated intensity map...")
    integrated_map = np.nansum(cube_data, axis=0)
    
    plt.figure(figsize=(10, 8))
    plt.imshow(integrated_map, origin='lower', cmap='viridis')
    plt.colorbar(label='Integrated Intensity (K)')
    plt.xlabel('RA (pixels)')
    plt.ylabel('DEC (pixels)')
    plt.title('Integrated Spectral Intensity Map')
    plt.tight_layout()
    
    output_file = Path(__file__).parent / 'test_integrated_map.png'
    plt.savefig(output_file, dpi=100)
    print(f"✓ Plot saved to: {output_file}")
    
    plt.show()
    
except Exception as e:
    print(f"✗ Error: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)
finally:
    hdul.close()

print("\n✓ Test completed successfully!")
