# Implementation Checklist - Velocity Axis Reconstruction

## Part 1: Datacube Velocity Axis Scaling

### Code Changes
- [x] Create `_get_spectral_axis_params()` function in gridding.py
  - [x] Extract VELOCITY from FITS table
  - [x] Extract DELTAV from FITS table
  - [x] Extract RESTFREQ from FITS table
  - [x] Extract VELDEF from FITS table
  - [x] Extract CRPIX1 from FITS header (new)
  - [x] Return 6 values including crpix1_spec

- [x] Update `create_spectral_datacube()` function
  - [x] Extract 6 spectral parameters (including crpix1_spec)
  - [x] Use deltav for CDELT3 (not hardcoded 1.0)
  - [x] Use crpix1_spec for CRPIX3 (not hardcoded 1.0)
  - [x] Print reference pixel and channel info
  - [x] Print velocity at key channels

- [x] Fix variable naming
  - [x] Change crpix_spec to crpix1_spec (location 1)
  - [x] Change crpix_spec to crpix1_spec (location 2)

### Testing
- [x] Verify syntax compiles (`python3 -m py_compile gridding.py`)
- [x] Create test script (`test_velocity_axis.py`)
- [x] Run test with real M51 data
- [x] Verify output shows correct velocity axis
- [x] Verify CRPIX1 extraction from header
- [x] Verify default to 1.0 when CRPIX1 missing

### Documentation
- [x] Create `CRPIX1_FIX.md` - Technical details
- [x] Create `VELOCITY_AXIS_RECONSTRUCTION_STATUS.md` - Complete status
- [x] Create `visualize_velocity_axis.py` - Visualization script
- [x] Create `velocity_axis_visualization.png` - Generated plots

## Part 2: CRPIX1 Preservation in Single-HDU Files

### Code Changes
- [x] Modify `combine_fits_files()` in basic_io.py
  - [x] Create spatial_wcs_keys list
  - [x] Create spectral_wcs_keys list (with CRPIX1)
  - [x] Create spectral_param_keys list
  - [x] Combine into all_keys_to_copy
  - [x] Copy all keys from first spectrum HDU to combined HDU

### Testing
- [x] Verify syntax compiles (`python3 -m py_compile basic_io.py`)
- [x] Create test script (`test_crpix1_preservation.py`)
- [x] Run test with real multi-HDU M51 file
- [x] Verify CRPIX1 is extracted from original file
- [x] Verify CRPIX1 is copied to combined file
- [x] Verify other spectral parameters preserved
- [x] Verify file can be written and read back
- [x] Confirm table columns with spectral data are preserved

### Investigation
- [x] Check where spectral parameters are stored
  - [x] CRPIX1 is in header keywords ✓
  - [x] VELOCITY, DELTAV, RESTFREQ, VELDEF are in table columns ✓
- [x] Verify multi-HDU file has CRPIX1 (~505.85)
- [x] Verify single-HDU file originally didn't have CRPIX1
- [x] Confirm concatenation preserves table columns

### Documentation
- [x] Create `CRPIX1_PRESERVATION_COMPLETE.md` - CRPIX1 fix details
- [x] Create `IMPLEMENTATION_COMPLETE.md` - Overall summary
- [x] Create `IMPLEMENTATION_CHECKLIST.md` - This file

## Integration & Workflow

### How They Work Together
```
Input: Multi-HDU FITS file
  ↓
[STEP 1] Combine FITS files (basic_io.py)
  - Concatenates all observations
  - Preserves CRPIX1 header keyword ✓
  - Preserves VELOCITY/DELTAV/RESTFREQ/VELDEF columns ✓
  - Output: Single-HDU FITS file with all parameters
  ↓
[STEP 2] Create datacube (gridding.py)
  - Reads single-HDU FITS file
  - Extracts spectral parameters including CRPIX1 ✓
  - Uses proper velocity spacing (DELTAV) ✓
  - Uses proper reference pixel (CRPIX1) ✓
  - Grids to proper velocity axis
  - Output: 3D datacube with correct velocity axis
```

## User Workflow

### Old Way (Broken)
```bash
$ combine_fits --input list.txt --output m51.fits --single-hdu
  # CRPIX1 was lost! ✗

$ oi_zeigt create-datacube --fits m51.fits
  # Velocity axis was wrong! ✗
  # CRPIX3=1.0, CDELT3=1.0
```

### New Way (Fixed)
```bash
$ combine_fits --input list.txt --output m51.fits --single-hdu
  # CRPIX1 now preserved! ✓

$ oi_zeigt create-datacube --fits m51.fits
  # Velocity axis is correct! ✓
  # CRPIX3~506 (from FITS file)
  # CDELT3=500 m/s (from FITS file)
```

## Files Status

### Code Files (Modified)
- ✅ `src/oi_zeigt/mapping/gridding.py`
  - Lines ~305-378: New `_get_spectral_axis_params()` function
  - Lines ~950+: Updated `create_spectral_datacube()` function

- ✅ `src/oi_zeigt/basic_io.py`
  - Lines ~219-233: Updated `combine_fits_files()` function

### Documentation Files (Created)
- ✅ `VELOCITY_AXIS_QUICK_SUMMARY.md`
- ✅ `CRPIX1_FIX.md`
- ✅ `VELOCITY_AXIS_RECONSTRUCTION_STATUS.md`
- ✅ `CRPIX1_PRESERVATION_COMPLETE.md`
- ✅ `IMPLEMENTATION_COMPLETE.md`
- ✅ `IMPLEMENTATION_CHECKLIST.md` (this file)

### Test Files (Created)
- ✅ `test_velocity_axis.py` - Tests datacube velocity axis
- ✅ `test_crpix1_preservation.py` - Tests CRPIX1 preservation
- ✅ `check_spectral_params_location.py` - Analyzes parameter storage
- ✅ `visualize_velocity_axis.py` - Generates velocity plots
- ✅ `velocity_axis_visualization.png` - Generated visualization

## Quality Assurance

### Code Quality
- [x] All syntax checked and validated
- [x] No compilation errors
- [x] Functions have proper docstrings
- [x] Changes are backward compatible
- [x] Error handling for missing parameters

### Testing
- [x] Unit tests created for both changes
- [x] Tests pass with real M51 data
- [x] File I/O tested (write and read back)
- [x] Parameter preservation verified
- [x] Edge cases handled (missing CRPIX1, etc.)

### Documentation
- [x] Changes documented clearly
- [x] Technical details explained
- [x] User workflow documented
- [x] Examples provided
- [x] Visualization generated

## Validation Results

### Test 1: Velocity Axis Extraction
```
✓ VELOCITY parameter: 470000 m/s (470.00 km/s)
✓ DELTAV parameter: 500 m/s (0.5000 km/s)
✓ RESTFREQ parameter: 1.9005369e+12 Hz
✓ VELDEF parameter: RADI-LSR
✓ CRPIX1 parameter: ~506 (when present in file)
✓ Datacube created: 1264×12×11 voxels
✓ WCS headers set correctly
```

### Test 2: CRPIX1 Preservation
```
✓ Original file has CRPIX1: 505.850218558238
✓ Combined file has CRPIX1: 505.850218558238
✓ Value correctly preserved
✓ File written and read back successfully
✓ No data loss in concatenation
```

## Deployment Readiness

- [x] Code changes complete and tested
- [x] No breaking changes to existing API
- [x] Backward compatible (handles missing CRPIX1)
- [x] Documentation complete
- [x] Tests included
- [x] User instructions provided

## Next Steps for User

1. **Regenerate single-HDU FITS** (optional, to get CRPIX1)
   ```bash
   combine_fits --input filelist.txt --output m51.fits --single-hdu
   ```

2. **Create datacubes** (will now use proper velocity axis)
   ```bash
   oi_zeigt create-datacube --fits m51.fits
   ```

3. **Verify output** (check CRPIX3 and CDELT3 in header)
   ```python
   from astropy.io import fits
   with fits.open('datacube.fits') as hdul:
       print(hdul[0].header['CRPIX3'])  # Should be ~506
       print(hdul[0].header['CDELT3'])  # Should be 500.0
   ```

4. **Use velocity coordinates** (now physically accurate!)
   ```python
   # Velocity at any channel i
   v = 470000 + (i - 505) * 500  # m/s
   ```

## Sign-Off

✅ **Implementation Complete**
✅ **All Tests Passing**
✅ **Documentation Complete**
✅ **Ready for Production Use**

Both the velocity axis scaling problem and the CRPIX1 preservation issue have been successfully resolved!
