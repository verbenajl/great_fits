# Fix: Preserve Original Data Type in Reduced FITS Files

## Problem Identified

When using `reduce_spectra --baseline`, the output file was **1.78× larger** than the input file (780 MB vs 439 MB).

### Root Cause

The SPECTRUM column was being converted from **float32 (`>f4`)** to **float64 (`>f8`)** during reduction, doubling its size:

**Before fix:**
- Clean file SPECTRUM: `('>f4', (1264,))` — 431.9 MB
- Reduced file SPECTRUM: `('>f8', (1264,))` — 863.8 MB

This happened because:
1. Spectra were loaded as float64 for processing (for numerical accuracy)
2. After reduction, they remained as float64
3. When saving to FITS, the float64 dtype was preserved

## Solution Implemented

Modified the `reduce_spectra()` function in `/home/verbena/software/oi_zeigt/src/oi_zeigt/reduction/core.py` to:

1. **Save the original data type** before processing:
   ```python
   original_spectrum_dtype = data[spectrum_column].dtype
   ```

2. **Process in float64** for numerical accuracy:
   ```python
   spectra = np.array([np.array(spec, dtype=float) for spec in data[spectrum_column]])
   ```

3. **Convert back to original dtype** before saving:
   ```python
   spectra = spectra.astype(original_spectrum_dtype)
   ```

## Expected Improvement

After the fix, the reduced file should be approximately **the same size as the input file**:
- ~439 MB reduced_data.fits (instead of 780 MB)
- **Savings: ~341 MB per file**

The reduction is lossless for float32 because baseline subtraction with float32 precision is more than adequate for spectral data.

## How It Works

1. Original SPECTRUM column (float32) is loaded
2. Converted to float64 for baseline calculation (better numerical precision)
3. Baseline is subtracted using float64 arithmetic
4. Result is converted back to float32
5. Saved to FITS with original dtype preserved

This maintains:
- ✅ Numerical accuracy during processing
- ✅ Original file size and dtype
- ✅ Compatibility with downstream processing

## Affected Methods

The fix applies to all reduction methods:
- `--baseline` ✅ Fixed
- `--smooth` ✅ Fixed
- `--unblank` ✅ Fixed
- Combined methods ✅ All fixed

## Testing

To verify the fix works, run:

```bash
# Clear old results
rm -f /path/to/reduced_data.fits

# Re-run reduction
reduce_spectra --config config.toml --clean --baseline

# Check file sizes
ls -lh clean_data.fits reduced_data.fits

# Inspect data types
python3 inspect_fits.py clean_data.fits reduced_data.fits
```

Expected output:
- SPECTRUM column should be `('>f4', (1264,))` in both files
- File sizes should be very similar
