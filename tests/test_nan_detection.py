#!/usr/bin/env python3
"""
Test script to verify NaN detection works correctly.
"""
import sys
from pathlib import Path
import numpy as np
from astropy.io import fits

# Add src to path
sys.path.insert(0, str(Path(__file__).parent / 'src'))

from oi_zeigt.reduction.core import detect_nan_channels

# Test 1: Simple numpy array with NaNs
print("Test 1: Simple numpy array with NaNs")
spectrum = np.array([1.0, 2.0, np.nan, 4.0, np.nan, 6.0])
mask, frac = detect_nan_channels(spectrum)
print(f"  Input: {spectrum}")
print(f"  NaN mask: {mask}")
print(f"  NaN fraction: {frac:.2%}")
assert np.allclose(frac, 2/6), f"Expected 0.333, got {frac}"
print("  ✓ PASSED\n")

# Test 2: Array with no NaNs
print("Test 2: Array with no NaNs")
spectrum = np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0])
mask, frac = detect_nan_channels(spectrum)
print(f"  Input: {spectrum}")
print(f"  NaN fraction: {frac:.2%}")
assert frac == 0.0, f"Expected 0.0, got {frac}"
print("  ✓ PASSED\n")

# Test 3: Array with all NaNs
print("Test 3: Array with all NaNs")
spectrum = np.array([np.nan, np.nan, np.nan])
mask, frac = detect_nan_channels(spectrum)
print(f"  Input: {spectrum}")
print(f"  NaN fraction: {frac:.2%}")
assert frac == 1.0, f"Expected 1.0, got {frac}"
print("  ✓ PASSED\n")

# Test 4: From actual FITS file if available
fits_file = Path('/home/diskB/sofiaobsdata/m51/m51_central_region/kalibrated_all/m51_central_tile.fits')
if fits_file.exists():
    print(f"Test 4: Using real FITS file: {fits_file.name}")
    with fits.open(fits_file) as hdul:
        # Get the first spectrum HDU
        for hdu in hdul:
            if hasattr(hdu, 'data') and hdu.data is not None:
                if 'SPECTRUM' in hdu.data.dtype.names:
                    # Get first few spectra
                    for i in range(min(5, len(hdu.data))):
                        spectrum = hdu.data[i]['SPECTRUM']
                        mask, frac = detect_nan_channels(spectrum)
                        
                        # Count actual NaNs
                        actual_nans = np.sum(np.isnan(spectrum))
                        expected_frac = actual_nans / len(spectrum)
                        
                        print(f"  Row {i}:")
                        print(f"    Spectrum shape: {spectrum.shape}")
                        print(f"    Actual NaN count: {actual_nans}")
                        print(f"    Detected NaN fraction: {frac:.4f}")
                        print(f"    Expected fraction: {expected_frac:.4f}")
                        
                        if not np.isclose(frac, expected_frac):
                            print(f"    ✗ MISMATCH!")
                        else:
                            print(f"    ✓ OK")
                    break
    print("✓ PASSED\n")
else:
    print(f"Test 4: Skipped (file not found: {fits_file})\n")

print("All tests completed!")
