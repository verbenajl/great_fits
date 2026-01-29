# Spectral Extraction Implementation - Summary

## What Was Implemented

Complete velocity-based spectral extraction feature for the reduction pipeline. Users can now extract spectra to a specific velocity range and apply baseline subtraction, with all parameters specified as **absolute velocities in km/s**.

## Configuration

Add these parameters to `config.toml` under `[reduction]`:

```toml
[reduction]
baseline = 3           # Polynomial order
window = [450, 500]   # Baseline window in km/s (absolute velocities)
extract = [350, 700]  # Extraction range in km/s (absolute velocities)
```

## How It Works

### Example: M51 Data

**Input:**
- 31,822 spectra × 1,264 channels
- Full velocity range: 217.6 - 849.1 km/s
- Channel spacing: 500 m/s per channel

**Configuration:**
- `extract=[350, 700]` km/s
- `window=[450, 500]` km/s

**Processing Steps:**

1. **Extraction**: Map velocities to channels
   - 350 km/s → channel 265
   - 700 km/s → channel 965
   - Extract 700 channels from each spectrum

2. **Baseline Subtraction**: Convert baseline window to extracted coordinates
   - Window [450, 500] km/s → channels [465, 565] in full spectrum
   - In extracted spectrum: channels [200, 300] (relative to start)
   - Fit 3rd-order polynomial
   - Subtract from all 700 extracted channels

3. **Output**: Write FITS with
   - Spectra: 31,822 × 700 channels
   - VELOCITY_AXIS column
   - All 52 original columns preserved

**Result:**
- 700 channels per spectrum (velocity range 350.1 - 699.6 km/s)
- Baseline subtracted (mean ≈ -4.77)
- File size: ~300 MB (compared to ~505 MB original)

## Key Features

✅ **Absolute velocity parameters** - Specify extraction and baseline window in km/s, independent of data  
✅ **Automatic channel mapping** - Velocities automatically converted to channel indices using FITS header parameters  
✅ **Flexible baseline window** - Baseline window can be anywhere in velocity space, even outside extraction range  
✅ **Smart coordinate adjustment** - Baseline window automatically adjusted to extracted spectrum coordinate system  
✅ **Velocity axis preserved** - VELOCITY_AXIS column reflects extracted channel range  
✅ **Column preservation** - All original FITS columns maintained  

## Test Files

### `test_extraction_from_config.py`
Full end-to-end test with real M51 data (31,822 spectra)
- Loads config from config.toml
- Applies extraction and baseline reduction
- Verifies output FITS file
- **Result: ✅ PASSED**

### `demo_extraction_parameters.py`
Educational demonstration showing parameter mapping
- Shows how km/s velocities map to channel indices
- Illustrates baseline window adjustment
- Explains processing order

### `test_reduced_extracted.fits`
Test output file (300 MB)
- 31,822 spectra × 700 channels
- VELOCITY_AXIS column present
- Baseline subtracted
- Can be used for downstream analysis

## Code Changes

**Modified file:** `src/oi_zeigt/reduction/core.py`

1. **`reduce_spectra_from_config()`** - Parses km/s parameters from config
2. **`reduce_spectra()`** - Handles extraction and baseline window conversion
3. **Helper functions** - Velocity-to-channel conversion utilities

**New documentation:**
- `EXTRACTION_VELOCITY_IMPLEMENTATION.md` - Technical reference
- `demo_extraction_parameters.py` - Educational walkthrough

## Ready for Production

The feature is fully implemented and tested:
- ✅ Configuration parsing works
- ✅ Velocity-to-channel conversion accurate
- ✅ Extraction produces correct output
- ✅ Baseline window adjustment correct
- ✅ VELOCITY_AXIS column generated correctly
- ✅ File integrity maintained
- ✅ All 31,822 spectra processed successfully

## Next Steps

1. **Use with real science data:**
   ```python
   from oi_zeigt.reduction.core import reduce_spectra_from_config
   output_path = reduce_spectra_from_config(
       config_path='config.toml',
       output_path='reduced_data.fits'
   )
   ```

2. **Adjust extraction range as needed:**
   - Modify `extract` values in config.toml
   - Re-run reduction
   - Parameters are flexible for different science requirements

3. **Combine with other reduction steps:**
   - Works alongside baseline subtraction
   - Compatible with smoothing if needed
   - Integrates with unblank filtering

4. **Apply PCA correction:**
   - Use mission-specific PCA bases (from earlier implementation)
   - Extract spectra to same velocity range
   - Compute/subtract PCA components
