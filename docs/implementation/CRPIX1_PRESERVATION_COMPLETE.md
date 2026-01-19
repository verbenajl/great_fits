# CRPIX1 Preservation in Single-HDU FITS Files - Implementation Complete

## ✅ PROBLEM SOLVED

When combining multi-HDU FITS files into single-HDU files, the critical **CRPIX1** parameter (reference pixel for spectral axis) was being **lost**.

## What Was Fixed

### File: `src/oi_zeigt/basic_io.py`

**Function: `combine_fits_files()`** (lines ~219-233)

Changed the header keyword copying from a hardcoded list:
```python
# OLD CODE (missing CRPIX1 and other important keywords)
for key in ['CRVAL2', 'CRVAL3', 'CRPIX2', 'CRPIX3', 'CTYPE2', 'CTYPE3', 'CUNIT2', 'CUNIT3']:
    if key in first_spectrum_hdu.header:
        combined_bintable_hdu.header[key] = first_spectrum_hdu.header[key]
```

To a comprehensive organized list:
```python
# NEW CODE (includes CRPIX1 and all WCS parameters)
spatial_wcs_keys = ['CRVAL2', 'CRVAL3', 'CRPIX2', 'CRPIX3', 'CTYPE2', 'CTYPE3', 'CUNIT2', 'CUNIT3', 'CDELT2', 'CDELT3']
spectral_wcs_keys = ['CRVAL1', 'CRPIX1', 'CTYPE1', 'CUNIT1', 'CDELT1']  # ← ADDED CRPIX1!

all_keys_to_copy = spatial_wcs_keys + spectral_wcs_keys + spectral_param_keys

for key in all_keys_to_copy:
    if key in first_spectrum_hdu.header:
        combined_bintable_hdu.header[key] = first_spectrum_hdu.header[key]
```

## What Gets Preserved Now

### Header Keywords (WCS Parameters)
✅ **Spatial axes** (RA/Dec):
- CRVAL2, CRVAL3 (reference coordinates)
- CRPIX2, CRPIX3 (reference pixels)
- CTYPE2, CTYPE3 (axis types)
- CUNIT2, CUNIT3 (axis units)
- CDELT2, CDELT3 (pixel scales)

✅ **Spectral axis** (frequency):
- **CRVAL1** (reference frequency)
- **CRPIX1** (reference pixel) ← **CRITICAL - NOW PRESERVED!**
- CTYPE1 (axis type: 'FREQ')
- CUNIT1 (axis unit: 'Hz')
- CDELT1 (frequency spacing per pixel)

### Table Columns (Observational Parameters)
✅ **Automatically preserved** when tables are concatenated:
- **VELOCITY** (LSR radial velocity, m/s)
- **DELTAV** (velocity spacing per channel, m/s)
- **RESTFREQ** (rest frequency, Hz)
- **VELDEF** (velocity definition, e.g., 'RADI-LSR')
- Plus 45 other observational columns

## Test Results

### Original Multi-HDU File
```
File: m51_central_tile.fits
CRPIX1 = 505.850218558238  (in each MATRIX HDU)
VELOCITY = 470000.0  (m/s, in table column)
DELTAV = 500.0  (m/s, in table column)
RESTFREQ = 1.9005369e+12  (Hz, in table column)
VELDEF = 'RADI-LSR'  (in table column)
```

### After Single-HDU Combination
```
✓ CRPIX1 = 505.850218558238  [PRESERVED]
✓ VELOCITY = 470000.0  [PRESERVED via table]
✓ DELTAV = 500.0  [PRESERVED via table]
✓ RESTFREQ = 1.9005369e+12  [PRESERVED via table]
✓ VELDEF = 'RADI-LSR'  [PRESERVED via table]
✓ Successfully persisted to disk and read back
```

## How It Works

### Multi-HDU File Structure
```
HDU 0: PRIMARY (minimal header)
HDU 1: MATRIX (14,658 observations, 49 columns)
  Header:  CRPIX1=505.85, CRVAL1=0.0, CTYPE1=FREQ, etc.
  Columns: VELOCITY, DELTAV, RESTFREQ, VELDEF, SPECTRUM, ...
HDU 2: MATRIX (20,510 observations)
  ...
HDU 3: MATRIX (3,878 observations)
  ...
HDU 4: MATRIX (31,822 observations)
  ...
Total: 70,868 observations
```

### Single-HDU Output After Combination
```
HDU 0: PRIMARY (from first file)
HDU 1: SPECTRA (combined from all MATRIX HDUs)
  Header:  CRPIX1=505.85 ✓ (copied from first HDU)
  Columns: All 70,868 observations concatenated
           VELOCITY, DELTAV, RESTFREQ, VELDEF, SPECTRUM, ...
```

## Impact on Velocity Axis Reconstruction

With CRPIX1 now preserved in single-HDU files, the velocity axis can be properly reconstructed:

$$v(i) = \text{CRVAL3} + (i - (\text{CRPIX3} - 1)) \times \text{CDELT3}$$

For M51 data:
- CRVAL3 = 470,000 m/s (from VELOCITY column)
- CRPIX3 = 505.85 (from CRPIX1 header keyword, now preserved!)
- CDELT3 = 500 m/s (from DELTAV column)

**Channel 505**: v = 470,000 + (505 - 505) × 500 = **470,000 m/s** ✓

## Files Modified
✅ `src/oi_zeigt/basic_io.py` - Line ~219-233 in `combine_fits_files()` function

## Files Created (for testing/documentation)
✅ `test_crpix1_preservation.py` - Comprehensive test script
✅ `check_spectral_params_location.py` - Parameter location analysis

## Usage

The fix is **automatic** - no changes needed to existing code:

```python
from src.oi_zeigt.basic_io import combine_fits_from_list

# Combine with single-HDU option (now preserves CRPIX1!)
combined = combine_fits_from_list("input_list.txt", "output.fits", single_hdu=True)

# CRPIX1 is now in the output file!
with fits.open("output.fits") as hdul:
    crpix1 = hdul[1].header.get('CRPIX1')  # Will be ~505.85 for M51
```

## Verification

Run the test script to verify CRPIX1 preservation:
```bash
python3 test_crpix1_preservation.py
```

Expected output:
```
✓ Original CRPIX1 found: 505.850218558238
✓ Combined CRPIX1 found: 505.850218558238
✓ CRPIX1 CORRECTLY PRESERVED!
```

## Summary

✅ **CRPIX1 is now preserved** when combining FITS files into single-HDU format
✅ **All WCS keywords are now copied** systematically (not just RA/Dec)
✅ **All table columns are preserved** (including spectral parameters)
✅ **Backward compatible** - existing code works without changes
✅ **Tested and verified** with real M51 data

Single-HDU FITS files created with the updated code will have all the information needed for proper velocity axis reconstruction!
