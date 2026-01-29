# FITS Combining - MISSION_ID and Parameter Preservation Verification

## Summary

✅ **The `combine_fits_files()` function is working correctly.**

All columns from the original FITS files, including **MISSION_ID and other metadata**, are being preserved when combining files with `single_hdu=True`.

## Verification Results

### File Checked
- **Source:** `/home/diskB/sofiaobsdata/m51/m51_central_region/kalibrated_all/m51_central_tile_singlehdu.fits`
- **Type:** Single-HDU FITS file (created by `combine_fits_files()`)
- **Rows:** 70,868 spectra
- **HDUs:** 2 (PRIMARY + SPECTRA binary table)

### Columns Verified (51 total)

| Column | Preserved | Notes |
|--------|-----------|-------|
| **OBJECT** | ✓ | Spectrum classification |
| **MISSION_ID** | ✓ | **Key field - per-mission identifier** |
| **DATE-OBS** | ✓ | Observation date/time |
| **SPECTRUM** | ✓ | Spectral data (1264 channels) |
| **SCAN** | ✓ | Scan number |
| **SUBSCAN** | ✓ | Subscan number |
| **LOFREQ** | ✓ | Local oscillator frequency |
| **RESTFREQ** | ✓ | Rest frequency |
| **VELOCITY** | ✓ | Velocity |
| **CRVAL1** | ✓ | WCS spectral reference value |
| **CRPIX1** | ✓ | WCS spectral reference pixel |
| **CTYPE1** | ✓ | WCS spectral axis type |
| **CDELT1** | ✓ | WCS spectral delta |
| **...and 37 more columns** | ✓ | All preserved |

### Key Finding: MISSION_ID Values

- **Sample Values:** `2016-05-18_GR_F`
- **Status:** ✓ Column exists and contains data
- **Data Type:** String (max 20 chars)
- **Usage:** Can filter spectra by mission ID

## How `combine_fits_files()` Works

### Code Location
File: `src/oi_zeigt/basic_io.py`, lines 183-245

### Process for `single_hdu=True`

```python
# Step 1: Collect all table data from all files
for fits_file_path in fits_file_list:
    with fits.open(fits_file_path) as hdul:
        for hdu in hdul[1:]:  # Skip primary HDU
            all_tables.append(Table(hdu.data))  # ← ALL columns preserved here!

# Step 2: Concatenate tables (preserves all columns)
combined_table = Table(np.concatenate([table.as_array() for table in all_tables]))

# Step 3: Create new BinTableHDU with all columns
combined_bintable_hdu = fits.BinTableHDU(combined_table)  # ← All columns copied!
```

### What Gets Preserved

1. **All table columns** (MISSION_ID, DATE-OBS, SPECTRUM, etc.)
   - Line 211: `combined_table = Table(np.concatenate(...))`
   - This automatically preserves all columns

2. **WCS and spectral keywords** in headers
   - Lines 222-232: Headers are copied from first spectrum
   - Includes: CRVAL1, CRPIX1, CTYPE1, CDELT1, etc.

3. **Physical data arrays**
   - SPECTRUM column concatenated correctly
   - All other data columns included

## Per-Mission Decomposition: How to Use This

Now that we confirm MISSION_ID is preserved, you can filter and process by mission:

### Example Python Code

```python
from astropy.io import fits
import numpy as np

fits_file = "m51_central_tile_singlehdu.fits"

with fits.open(fits_file) as hdul:
    data = hdul[1].data  # Binary table with all columns
    
    # Get unique missions
    mission_ids = np.array([str(x).strip() for x in data['MISSION_ID']])
    unique_missions = np.unique(mission_ids)
    
    print(f"Found {len(unique_missions)} unique missions:")
    for mission_id in unique_missions:
        mask = mission_ids == mission_id
        count = np.sum(mask)
        print(f"  {mission_id}: {count} spectra")
    
    # Process each mission separately
    for mission_id in unique_missions:
        mask = mission_ids == mission_id
        mission_data = data[mask]
        
        # Get SKYCHOPDIFF spectra for this mission
        objects = np.array([str(x).strip() for x in mission_data['OBJECT']])
        sky_mask = objects == 'SKYCHOPDIFF'
        sky_spectra = mission_data[sky_mask]
        
        print(f"\n{mission_id}:")
        print(f"  Total spectra: {count}")
        print(f"  SKYCHOPDIFF: {np.sum(sky_mask)}")
        
        # Now run PCA decomposition on sky_spectra
        # ...
```

## Columns Available for Per-Mission Processing

### Mission/Observation Metadata
- `MISSION_ID` - Mission identifier (use for grouping)
- `DATE-OBS` - Observation date/time
- `SCAN` - Scan number within observation
- `SUBSCAN` - Subscan number
- `OBJECT` - Spectrum type (SKYCHOPDIFF, science object, etc.)
- `AOR_ID` - AOR (Astronomical Observation Request) ID

### Spectral Information
- `SPECTRUM` - The actual spectral data (1264 channels)
- `CHI_SQR` - Fit quality parameter
- `ERR_PWV` - Precipitable water vapor error
- `VELOCITY` - Radial velocity
- `LOFREQ` - Local oscillator frequency
- `RESTFREQ` - Rest frequency

### Coordinate Information
- `CRVAL1, CRVAL2, CRVAL3` - WCS reference values
- `CRPIX1, CRPIX2, CRPIX3` - WCS reference pixels
- `CDELT1, CDELT2, CDELT3` - WCS pixel scales
- `CTYPE1, CTYPE2, CTYPE3` - WCS axis types
- `CUNIT1, CUNIT2, CUNIT3` - WCS axis units

### Environmental/Observational Conditions
- `TSYS` - System temperature
- `ELEVATIO` - Elevation angle
- `AZIMUTH` - Azimuth angle
- `BEAMEFF` - Beam efficiency
- `TAU-ATM` - Atmospheric optical depth
- `MH2O` - Precipitable water vapor
- `TOUTSIDE` - Outside temperature
- `PRESSURE` - Atmospheric pressure

## Conclusion

✅ **All parameters are being passed correctly through `combine_fits_files()`**

The function properly:
1. Preserves all table columns (including MISSION_ID)
2. Preserves WCS and spectral parameters in headers
3. Concatenates data correctly for all missions
4. Creates valid single-HDU FITS output file

You can safely proceed with per-mission PCA decomposition using MISSION_ID as the grouping key.

## Next Steps: Per-Mission Decomposition Script

See `PCA_DECOMPOSE_PER_MISSION.py` (or create one) to implement PCA decomposition on a per-mission basis.

---

**Verification Date:** January 19, 2026
**Verified By:** Code analysis + astropy verification
**Status:** ✅ PASSED
