# VELOCITY_AXIS Implementation Summary

## Overview

The `reduce_spectra()` function in `src/oi_zeigt/reduction/core.py` has been enhanced to automatically add a **VELOCITY_AXIS** column to the output FITS file. This allows for direct correlation between spectrum values and their corresponding velocities.

## What Was Added

### 1. Helper Functions (core.py)

**`_extract_spectral_params(hdul)`**
- Extracts spectral parameters from FITS file headers and table columns
- Returns: dictionary with `velo_ref`, `deltav`, `crpix1_spec`, and `nchans`
- Parameters extracted from:
  - VELOCITY column → reference velocity in m/s
  - DELTAV column → velocity spacing per channel in m/s
  - CRPIX1 header → reference pixel (1-indexed, FITS convention)
  - SPECTRUM column → number of spectral channels

**`_create_velocity_axis(velo_ref, deltav, crpix1_spec, nchans)`**
- Creates velocity axis array using WCS formula
- Calculation: `v[i] = velo_ref + (i - (crpix1_spec - 1)) * deltav`
- Returns: numpy array of shape (nchans,) with velocities in m/s

### 2. Modified `reduce_spectra()` Function

After applying reduction methods (baseline, smoothing, etc.):
1. Extracts spectral parameters from original FITS file
2. Creates velocity axis using extracted parameters
3. Tiles velocity axis for all spectra (creates [n_spectra, n_channels] array)
4. Adds VELOCITY_AXIS column to output FITS table

**Error Handling**: If spectral parameters cannot be extracted, a warning is issued but processing continues (graceful degradation)

## Usage

### Basic Usage
```python
from oi_zeigt.reduction.core import reduce_spectra
from astropy.io import fits

# Load FITS file
hdul = fits.open('input.fits')

# Apply reduction with velocity axis
methods = {
    'baseline': {'order': 2, 'window': (100, 120)},
    'smooth': {'window_size': 5}
}
reduced_hdul = reduce_spectra(hdul, methods=methods)

# Save output
reduced_hdul.writeto('reduced.fits', overwrite=True)
```

### Accessing Velocity Axis
```python
from astropy.io import fits
import matplotlib.pyplot as plt

# Load reduced FITS file
with fits.open('reduced.fits') as hdul:
    data = hdul[1].data
    spectrum = data['SPECTRUM'][0]
    velocity = data['VELOCITY_AXIS'][0]

# Plot spectrum vs velocity
plt.plot(velocity, spectrum)
plt.xlabel('Velocity (m/s)')
plt.ylabel('Intensity')
plt.show()
```

## Output FITS Structure

### Original Columns (Preserved)
- All original columns from input FITS
- SPECTRUM: reduced spectra [n_spectra, n_channels]

### New Columns (Added)
- VELOCITY_AXIS: velocity array [n_spectra, n_channels]
  - Data type: float64
  - Unit: m/s
  - Same velocity axis for all spectra (unless velocities vary per observation)

### Header Keywords (Copied)
- VELOCITY, DELTAV, VELDEF, RESTFREQ (spectral parameters)
- CRPIX1 (reference pixel)
- All other original header keywords

## Example Output

For the test case:
- VELOCITY = 470000 m/s (reference)
- DELTAV = 500 m/s/channel
- CRPIX1 = 506 (reference pixel, 1-indexed)
- n_channels = 1264

The velocity axis covers:
```
Channel 0:    217500 m/s
Channel 505:  470000 m/s (reference channel)
Channel 1263: 849000 m/s
```

## Testing

Run the test suite:
```bash
cd /home/verbena/software/oi_zeigt
python3 test_reduce_spectra_velocity.py
```

Expected output:
```
✓ VELOCITY_AXIS column found in output
✓ Velocity axis values are correct!
```

## Integration with Pipeline

The velocity axis is **automatically added** whenever `reduce_spectra()` is called:

1. Via `reduce_spectra_from_config()` - automatic
2. Via direct `reduce_spectra()` call - automatic
3. Via CLI commands - automatic (if implemented)

No additional configuration needed. The velocity axis is extracted from the input FITS file and propagated to the output.

## Benefits

1. **No Additional Metadata Required**: Velocity axis is computed from existing FITS parameters
2. **Direct Spectrum-Velocity Correlation**: Each spectrum value directly maps to its velocity
3. **Plotting Ready**: Can immediately plot spectrum vs velocity without recalculating axis
4. **Gridding Compatible**: Velocity axis structure matches SPECTRUM column structure
5. **Graceful Degradation**: If spectral parameters missing, processing continues without velocity axis

## File Size Impact

- Adds ~10 MB per 70,000 spectra at 1264 channels (float64)
- Can be reduced to ~5 MB if stored as float32
- Alternative: Store only header keywords and reconstruct on-the-fly (more computation)

## Future Enhancements

1. Add option to store as float32 instead of float64 to reduce file size
2. Add option to skip VELOCITY_AXIS column if not needed
3. Support per-spectrum velocity variations (if observations have different VELOCITY values)
4. Add velocity-based masking (e.g., extract specific velocity ranges)

## Files Modified

- `src/oi_zeigt/reduction/core.py`: Added helper functions and VELOCITY_AXIS column creation
- `VELOCITY_AXIS_IMPLEMENTATION.md`: Implementation design document
- `test_reduce_spectra_velocity.py`: Test script

## Status

✅ **IMPLEMENTED AND TESTED**
- Helper functions working correctly
- VELOCITY_AXIS column being added to output
- Velocity values computed accurately
- Test suite passing
