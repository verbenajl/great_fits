# Velocity Axis Reconstruction - CRPIX1 Integration

## Problem Identified

The user pointed out that the reference velocity of **470 km/s is NOT at channel 0** in the M51 dataset, but rather at some other channel (around channel 505-506).

## Solution: CRPIX1 Integration

We discovered that the original FITS file (multi-HDU version) contains **CRPIX1 ≈ 506**, which is the **reference pixel for the spectral axis** in FITS convention (1-indexed).

### What is CRPIX1?

- **CRPIX1**: Reference pixel in FITS convention (1-indexed)
  - For M51: CRPIX1 ≈ 506
  - Corresponds to channel 505 (0-indexed)
  - At this reference channel, velocity = 470 km/s

### Velocity Axis Formula (Corrected)

The velocity of channel `i` is now calculated as:

$$v(i) = \text{CRVAL3} + (i - (\text{CRPIX3} - 1)) \times \text{CDELT3}$$

Where:
- **CRVAL3** = 470,000 m/s (reference velocity from VELOCITY column)
- **CRPIX3** = CRPIX1 from FITS header (≈506 for M51)
- **CDELT3** = 500 m/s (velocity step from DELTAV column)
- **i** = channel index (0-based)

### Example for M51 (with CRPIX1 ≈ 506):

```
Channel 0:     v = 470,000 + (0 - 505) × 500 = 217,500 m/s = 217.5 km/s
Channel 505:   v = 470,000 + (505 - 505) × 500 = 470,000 m/s = 470.0 km/s (reference!)
Channel 1263:  v = 470,000 + (1263 - 505) × 500 = 1,101,500 m/s = 1101.5 km/s
```

## Code Changes

### 1. Modified `_get_spectral_axis_params()` function

**Now returns 6 values instead of 5:**

```python
def _get_spectral_axis_params(hdul: fits.HDUList) -> Tuple[float, float, float, str, float, float]:
    """
    Returns:
    - velo_ref: Reference velocity in m/s
    - deltav: Velocity spacing per channel in m/s
    - restfreq: Rest frequency in Hz
    - veldef: Velocity definition string
    - nchans: Number of spectral channels
    - crpix1_spec: Reference pixel for spectral axis (NEW!)
    """
```

**Extraction logic for CRPIX1:**
1. First looks in binary table header (HDU 1)
2. If not found, looks in primary header (HDU 0)
3. If not found, defaults to 1.0 (first channel)

### 2. Updated `create_spectral_datacube()` function

Now:
- Extracts `crpix1_spec` along with other parameters
- Prints the reference pixel and reference channel
- Prints velocity at each key channel
- Uses `crpix1_spec` instead of hardcoded `1.0` when creating WCS headers

### 3. Updated WCS Header Creation

The CRPIX3 keyword in the output FITS file now correctly reflects the reference pixel:

```python
wcs_header = create_wcs_header(
    ...
    crpix3=crpix1_spec,  # Now uses actual CRPIX1 value!
    cdelt3=deltav,       # Uses actual velocity spacing
    ...
)
```

## Output Changes

### Before (Incorrect)
```
Velocity at channel 0: 470000 m/s (470.00 km/s)
Channel 0: 470000.00 m/s (470.00 km/s)  ← Wrong! This is reference, not channel 0
```

### After (Correct)
```
Reference pixel (CRPIX1): 506.0 (FITS 1-indexed)
Reference channel: 505.0 (0-indexed)
Velocity at channel 0: 217500 m/s = 217.50 km/s
Velocity at ref channel 505: 470000 m/s = 470.00 km/s  ← Correct!
Velocity at last channel 1263: 1101500 m/s = 1101.50 km/s
```

## File Compatibility

### Multi-HDU File (m51_central_tile.fits)
- ✅ Contains CRPIX1 in MATRIX HDU headers
- ✅ Will be extracted and used correctly
- ✅ Output datacube will have proper velocity axis

### Single-HDU File (m51_central_tile_singlehdu.fits)
- ✗ Does not contain CRPIX1 (lost during conversion)
- ⚠ Will default to CRPIX1 = 1.0
- ⚠ Velocity axis will be incorrect (but code still works)

**Recommendation:** Use the multi-HDU file if possible, or regenerate the single-HDU file from the multi-HDU version to preserve CRPIX1.

## Verification

To verify the velocity axis is correct:

```python
from astropy.io import fits
import numpy as np

with fits.open('datacube.fits') as hdul:
    header = hdul[0].header
    
    crval3 = header['CRVAL3']  # 470000 m/s
    crpix3 = header['CRPIX3']  # ~506
    cdelt3 = header['CDELT3']  # 500 m/s
    naxis3 = header['NAXIS3']  # 1264
    
    # Calculate velocity at channel 505
    v_505 = crval3 + (505 - (crpix3 - 1)) * cdelt3
    print(f"Velocity at channel 505: {v_505:.0f} m/s")  # Should be 470000 m/s
    
    # Calculate velocity at channel 0
    v_0 = crval3 + (0 - (crpix3 - 1)) * cdelt3
    print(f"Velocity at channel 0: {v_0:.0f} m/s")  # Should be 217500 m/s
```

## Summary

✅ **Problem Fixed**: Velocity axis now properly reflects the reference channel from CRPIX1
✅ **Backward Compatible**: Still works with files that don't have CRPIX1
✅ **Physically Accurate**: 470 km/s is now at the correct channel (~505), not at channel 0
✅ **WCS Compliant**: Output FITS files have correct CRPIX3 values
