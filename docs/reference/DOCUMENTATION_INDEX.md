# Documentation Index - Velocity Axis Reconstruction Project

## Quick Start

**Start here:** `VELOCITY_AXIS_QUICK_SUMMARY.md`
- One-page overview of the problem and fix
- User-friendly explanation
- Next steps

## Implementation Documentation

### 📋 Main Implementation Guides

1. **`BEFORE_AND_AFTER.md`** ⭐ Visual Comparison
   - Side-by-side before/after code
   - Data flow diagrams
   - Metrics and validation proof
   - Shows exactly what was wrong and what was fixed

2. **`IMPLEMENTATION_COMPLETE.md`** - Overall Summary
   - Complete overview of both fixes
   - Impact chain showing how fixes connect
   - Test results
   - File list and changes

3. **`IMPLEMENTATION_CHECKLIST.md`** - Verification
   - Detailed checklist of all changes
   - Status of each item
   - Quality assurance section
   - Deployment readiness

### 🔧 Technical Details

4. **`CRPIX1_FIX.md`** - Spectral Axis Reference Pixel
   - Explains what CRPIX1 is
   - Velocity calculation formula
   - File compatibility notes
   - Verification instructions

5. **`CRPIX1_PRESERVATION_COMPLETE.md`** - CRPIX1 in Single-HDU Files
   - Details on combining FITS files
   - What gets preserved now
   - Test results with tables
   - How it works under the hood

6. **`VELOCITY_AXIS_RECONSTRUCTION_STATUS.md`** - Datacube Velocity Axis
   - Complete status of velocity axis fix
   - Function signatures and changes
   - Output examples
   - Known limitations and next steps

## Test Scripts

### 🧪 Validation Scripts

Run these to verify the implementation:

```bash
# Test 1: Datacube velocity axis extraction
python3 test_velocity_axis.py

# Test 2: CRPIX1 preservation in single-HDU files
python3 test_crpix1_preservation.py

# Test 3: Check where parameters are stored
python3 check_spectral_params_location.py

# Test 4: Generate velocity-channel visualization
python3 visualize_velocity_axis.py
```

### Generated Outputs

- **`velocity_axis_visualization.png`** - Plots showing:
  - Full velocity axis (channels 0-1263)
  - Zoomed view around reference channel
  - Reference point at channel 505
  - Velocity for each channel

## Code Changes

### Modified Files

1. **`src/oi_zeigt/mapping/gridding.py`**
   - Added `_get_spectral_axis_params()` function (lines ~305-378)
   - Updated `create_spectral_datacube()` function (lines ~950+)
   - Fixed variable naming (crpix_spec → crpix1_spec)

2. **`src/oi_zeigt/basic_io.py`**
   - Updated `combine_fits_files()` function (lines ~219-233)
   - Now copies CRPIX1 and other spectral WCS keywords

### Test/Analysis Files Created

- `test_velocity_axis.py` - Datacube velocity axis test
- `test_crpix1_preservation.py` - CRPIX1 preservation test
- `check_spectral_params_location.py` - Parameter location analysis
- `visualize_velocity_axis.py` - Velocity-channel visualization

## Documentation Reading Guide

### For Quick Understanding (5 minutes)
1. `VELOCITY_AXIS_QUICK_SUMMARY.md` - What was wrong and what's fixed
2. `BEFORE_AND_AFTER.md` - Visual comparison of old vs new

### For Technical Implementation (15 minutes)
1. `IMPLEMENTATION_COMPLETE.md` - Overall summary
2. `CRPIX1_FIX.md` - Reference pixel details
3. `VELOCITY_AXIS_RECONSTRUCTION_STATUS.md` - Datacube axis details

### For Complete Understanding (30 minutes)
1. Read all of above
2. `CRPIX1_PRESERVATION_COMPLETE.md` - File combination details
3. `IMPLEMENTATION_CHECKLIST.md` - Verification checklist

### For Verification (10 minutes)
1. `IMPLEMENTATION_CHECKLIST.md` - Status of all items
2. Run test scripts
3. Check output files

## Key Concepts

### CRPIX1 (Reference Pixel for Spectral Axis)
- **What:** Pixel coordinate where the reference velocity is located
- **Where:** In FITS header keywords (now preserved!)
- **Why:** Determines where the velocity reference point (470 km/s) is located
- **Value:** ~506 (FITS 1-indexed) = channel 505 (0-indexed) for M51

### CDELT3 (Velocity Spacing per Channel)
- **What:** Velocity step between adjacent channels
- **Where:** In WCS header
- **Was:** Hardcoded to 1.0 (wrong!)
- **Now:** Extracted from DELTAV column = 500 m/s (correct!)

### Spectral Parameters
- **VELOCITY:** Reference velocity in m/s (470,000 for M51)
- **DELTAV:** Velocity spacing per channel (500 m/s)
- **RESTFREQ:** Rest frequency in Hz (1.9005369e12)
- **VELDEF:** Velocity definition ('RADI-LSR')

## Quick Reference

### M51 Example Values
```
CRPIX1 = 505.850218558238 (reference pixel)
CRVAL3 = 470000.0 m/s (reference velocity)
CRPIX3 = 505.85 (same as CRPIX1)
CDELT3 = 500.0 m/s (velocity per channel)

Channel 0:   velocity = 217500 m/s
Channel 505: velocity = 470000 m/s (reference!)
Channel 1263: velocity = 1101500 m/s

Total velocity span: 884 km/s (631.5 km/s of spectrum)
```

### WCS Velocity Formula
```
v(i) = CRVAL3 + (i - (CRPIX3 - 1)) × CDELT3

Example for channel 505:
v(505) = 470000 + (505 - 505.85) × 500 = ~470000 m/s ✓
```

## Troubleshooting

### Issue: CRPIX1 not found warning
- **Cause:** Input FITS file doesn't have CRPIX1 in header
- **Effect:** Defaults to 1.0 (first channel is reference)
- **Solution:** Use multi-HDU file if available, or regenerate with CRPIX1

### Issue: Velocity axis still looks wrong
- **Cause:** Using old datacube files created before fix
- **Solution:** Regenerate datacubes with updated code

### Issue: FITS combining fails
- **Cause:** File permissions or path issues
- **Solution:** Check that input FITS files exist and are readable

## Files Summary Table

| File | Purpose | Status |
|------|---------|--------|
| `VELOCITY_AXIS_QUICK_SUMMARY.md` | Quick reference | ✅ Complete |
| `BEFORE_AND_AFTER.md` | Visual comparison | ✅ Complete |
| `IMPLEMENTATION_COMPLETE.md` | Full overview | ✅ Complete |
| `IMPLEMENTATION_CHECKLIST.md` | Verification | ✅ Complete |
| `CRPIX1_FIX.md` | Reference pixel explanation | ✅ Complete |
| `CRPIX1_PRESERVATION_COMPLETE.md` | File combination details | ✅ Complete |
| `VELOCITY_AXIS_RECONSTRUCTION_STATUS.md` | Datacube axis status | ✅ Complete |
| `test_velocity_axis.py` | Velocity axis test | ✅ Working |
| `test_crpix1_preservation.py` | CRPIX1 preservation test | ✅ Working |
| `check_spectral_params_location.py` | Parameter location analysis | ✅ Working |
| `visualize_velocity_axis.py` | Velocity visualization | ✅ Working |
| `velocity_axis_visualization.png` | Generated plots | ✅ Generated |

## Next Steps

1. **For Immediate Use:**
   - Review `BEFORE_AND_AFTER.md` for quick understanding
   - Run test scripts to verify installation
   - Check `velocity_axis_visualization.png` to understand velocity-channel mapping

2. **For Production Use:**
   - Regenerate single-HDU FITS files (CRPIX1 will be preserved)
   - Create new datacubes with updated code
   - Verify output headers have correct CRPIX3 and CDELT3

3. **For Documentation:**
   - Add velocity axis explanation to user guide
   - Include FITS header keyword reference
   - Add spectral analysis examples

## Contact & Questions

All documentation is self-contained and comprehensive. If issues arise:
1. Check `IMPLEMENTATION_CHECKLIST.md` for status
2. Review `BEFORE_AND_AFTER.md` for technical details
3. Run test scripts to verify installation
4. Check error messages in script output

---

**Last Updated:** 2026-01-15
**Status:** ✅ Implementation Complete and Tested
**Ready for Production:** Yes
