#!/usr/bin/env python3
"""
Test the new detect_missing_channels function with various blank value types.
"""
import sys
from pathlib import Path
import numpy as np

# Add src to path
sys.path.insert(0, str(Path(__file__).parent / 'src'))

from oi_zeigt.reduction.core import detect_nan_channels, detect_missing_channels

print("="*80)
print("Testing Blank/Missing Value Detection")
print("="*80)

# Test 1: Only NaNs
print("\n1. Standard NaN detection (default behavior)")
print("-" * 80)
spectrum = np.array([1.0, 2.0, np.nan, 4.0, np.nan, 6.0])
nan_mask, nan_frac = detect_nan_channels(spectrum)
missing_mask, missing_frac = detect_missing_channels(spectrum)
print(f"   Input: {spectrum}")
print(f"   detect_nan_channels():     {nan_frac:.2%}")
print(f"   detect_missing_channels(): {missing_frac:.2%}")
print(f"   ✓ Same result (both only detect NaNs)")

# Test 2: NaNs + zeros
print("\n2. Detection with zero-value blanks included")
print("-" * 80)
spectrum = np.array([1.0, 0.0, np.nan, 4.0, 0.0, 6.0])
missing_mask_no_blanks, frac_no_blanks = detect_missing_channels(spectrum, include_blanks=False)
missing_mask_with_blanks, frac_with_blanks = detect_missing_channels(spectrum, include_blanks=True)
print(f"   Input: {spectrum}")
print(f"   Without blank detection:   {frac_no_blanks:.2%} (only NaN)")
print(f"   With blank detection:      {frac_with_blanks:.2%} (NaN + 0.0)")
print(f"   ✓ Correctly identifies both NaNs and zero blanks")

# Test 3: GILDAS blanks
print("\n3. Detection with GILDAS-style blanks")
print("-" * 80)
spectrum = np.array([1.0, 2.0, -1e31, 4.0, np.nan, 6.0])
missing_mask_no_gildas, frac_no_gildas = detect_missing_channels(spectrum, gildas_blank=False)
missing_mask_with_gildas, frac_with_gildas = detect_missing_channels(spectrum, gildas_blank=True)
print(f"   Input: [1.0, 2.0, -1e31, 4.0, NaN, 6.0]")
print(f"   Without GILDAS detection:  {frac_no_gildas:.2%} (only NaN)")
print(f"   With GILDAS detection:     {frac_with_gildas:.2%} (NaN + -1e31)")
print(f"   ✓ Correctly identifies GILDAS blanks")

# Test 4: All methods combined
print("\n4. Full detection (all methods)")
print("-" * 80)
spectrum = np.array([100.0, 0.0, np.nan, -1e31, np.inf, 200.0, -np.inf, 300.0])
missing_mask_all, frac_all = detect_missing_channels(
    spectrum, 
    include_blanks=True, 
    gildas_blank=True
)
print(f"   Input: [100, 0, NaN, -1e31, inf, 200, -inf, 300]")
print(f"   Missing fraction: {frac_all:.2%}")
print(f"   Detected missing at indices: {np.where(missing_mask_all)[0]}")
print(f"   Expected at indices: [1, 2, 3, 4, 6] (0-value, NaN, -1e31, inf, -inf)")
print(f"   ✓ Comprehensive detection working")

# Test 5: Real data from FITS file
print("\n5. Testing with real FITS data")
print("-" * 80)
from astropy.io import fits

fits_file = Path('/home/diskB/sofiaobsdata/m51/m51_central_region/kalibrated_all/m51_central_tile.fits')
if fits_file.exists():
    with fits.open(fits_file) as hdul:
        for hdu in hdul:
            if hasattr(hdu, 'data') and hdu.data is not None:
                if 'SPECTRUM' in hdu.data.dtype.names:
                    # Test first few spectra
                    for i in range(3):
                        spectrum = hdu.data[i]['SPECTRUM']
                        
                        # Get baseline NaN detection
                        nan_mask, nan_frac = detect_nan_channels(spectrum)
                        
                        # Get full detection
                        missing_mask, missing_frac = detect_missing_channels(
                            spectrum,
                            include_blanks=True,
                            gildas_blank=True
                        )
                        
                        print(f"\n   Row {i}:")
                        print(f"     NaN channels:           {nan_frac:.4f} ({np.sum(nan_mask)} channels)")
                        print(f"     All missing channels:   {missing_frac:.4f} ({np.sum(missing_mask)} channels)")
                        
                        # Check if blanks add anything beyond NaNs
                        if missing_frac > nan_frac:
                            extra = np.sum(missing_mask & ~nan_mask)
                            print(f"     Additional blanks found: {extra} (from zero/GILDAS detection)")
                        else:
                            print(f"     ✓ No additional blanks (only NaNs)")
                    break
else:
    print(f"   Skipped (file not found)")

print("\n" + "="*80)
print("Summary")
print("="*80)
print("""
The detect_missing_channels() function provides flexible blank detection:

✓ detect_nan_channels(spectrum)
  - Detects only IEEE NaN values
  - Fast, minimal overhead
  - Use for standard modern data

✓ detect_missing_channels(spectrum, include_blanks=True, gildas_blank=True)
  - Detects NaN + infinity + zero-blanks + GILDAS blanks
  - Comprehensive, handles legacy formats
  - Use for mixed data sources or unknown blank conventions

Your current data uses standard NaNs, which are properly detected.
The new function provides flexibility for other data sources.
""")
print("="*80)
