# Summary: Velocity Axis Fix Implementation

## The Problem You Identified

You noticed: **"The 470 km/s velocity cannot be at channel 1, it must be at some other channel"**

This was **exactly right**. Investigation revealed:
- ✅ CRPIX1 ≈ 506 in the multi-HDU FITS file
- ✅ This means channel 505 (0-indexed) is the reference
- ✅ At channel 505, velocity = 470 km/s (verified!)

## What We Fixed

### The Old Code Had:
```python
cdelt3 = 1.0              # WRONG: treats channels as unitless
crpix3 = 1.0              # WRONG: assumes reference at channel 0
```

### The New Code Has:
```python
cdelt3 = deltav            # ✅ 500 m/s per channel
crpix3 = crpix1_spec       # ✅ ~506 from FITS header
```

## Implementation Details

### What Changed

**File: `src/oi_zeigt/mapping/gridding.py`**

1. **New Function: `_get_spectral_axis_params()`**
   - Now extracts 6 values (was 5)
   - **NEW**: Extracts CRPIX1 from FITS header
   - Searches table header first, then primary header
   - Defaults to 1.0 if missing (with warning)

2. **Updated: `create_spectral_datacube()`**
   - Gets all 6 parameters including crpix1_spec
   - Uses actual CDELT3 and CRPIX3 in WCS headers
   - Prints reference channel information

### Velocity Calculation

For any channel i:
$$v(i) = 470,000 + (i - 505) \times 500 \text{ m/s}$$

Example:
- Channel 0: v = 470,000 + (0 - 505) × 500 = 217,500 m/s (217.5 km/s)
- Channel 505: v = 470,000 + (505 - 505) × 500 = 470,000 m/s (470.0 km/s) ✅
- Channel 1263: v = 470,000 + (1263 - 505) × 500 = 1,101,500 m/s (1101.5 km/s)

## Current Status

✅ **Code Modified and Tested**
- Function `_get_spectral_axis_params()` extracts CRPIX1
- Function `create_spectral_datacube()` uses it in WCS headers
- Test runs successfully with M51 data
- Syntax verified (no compilation errors)

⚠️ **Single-HDU File Limitation**
- The single-HDU FITS file doesn't have CRPIX1
- Code defaults to CRPIX1=1.0 (with warning)
- Multi-HDU file has CRPIX1≈506 and will be used correctly

## Key Files

| File | Purpose |
|------|---------|
| `src/oi_zeigt/mapping/gridding.py` | Main implementation |
| `CRPIX1_FIX.md` | Detailed technical documentation |
| `VELOCITY_AXIS_RECONSTRUCTION_STATUS.md` | Complete implementation status |
| `visualize_velocity_axis.py` | Visualization of velocity-channel mapping |
| `velocity_axis_visualization.png` | Generated plots |

## Next: Regenerate Your Datacubes

To use the fixed velocity axis:

1. Use the **multi-HDU FITS file** if available (has CRPIX1)
2. Run the datacube creation script
3. Check that output has:
   - CRPIX3 ≈ 506 (or whatever CRPIX1 was in input)
   - CDELT3 = 500 m/s (not 1.0)
   - CRVAL3 = 470000 m/s

Your analysis will now work with proper velocity coordinates!

## Questions for You

1. **Should we recover CRPIX1** from the multi-HDU file when creating single-HDU versions?
2. **Do you want to regenerate** the datacubes with this fix?
3. **Should we add** CRPIX1 preservation to the FITS preparation step?

Let me know if you'd like to test with your actual data next!
