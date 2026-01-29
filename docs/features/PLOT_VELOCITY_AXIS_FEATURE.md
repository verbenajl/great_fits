# Plot Sample Spectra - Velocity Axis Feature

## Overview

The `plot_sample_spectra` function has been updated to automatically use the **VELOCITY_AXIS** column for the x-axis when available, with automatic fallback to channel numbering for backward compatibility.

## Features

### Automatic Velocity Axis Detection
- **If VELOCITY_AXIS exists**: X-axis is plotted as **Velocity (km/s)**
- **If VELOCITY_AXIS missing**: X-axis falls back to **Channel** numbering
- **Zero configuration required**: Works automatically with both reduced and original data

### Velocity Axis Conversion
- Velocity values are stored in FITS as m/s
- Automatically converted to km/s for human-readable plots
- Extracted spectra show the velocity range (e.g., 350-700 km/s)

## Usage

### Plot with extracted+reduced data (shows velocity axis):
```bash
plot_sample_spectra --fits test_reduced_extracted.fits --object M51CENTER --num-spectra 8 --output plot.png
```

**Result**: X-axis shows "Velocity (km/s)" from 350-700 km/s

### Plot with original data (shows channel axis):
```bash
plot_sample_spectra --fits original_data.fits --object M51CENTER --num-spectra 8 --output plot.png
```

**Result**: X-axis shows "Channel" (0 to 1263)

### Using config file:
```bash
plot_sample_spectra --config config.toml --object M51CENTER --num-spectra 20
```

## Implementation Details

### Code Changes
File: `src/oi_zeigt/cli.py` in `plot_sample_spectra()` function

**Added logic:**
1. Check if `VELOCITY_AXIS` column exists in FITS data
2. For each spectrum:
   - Extract the corresponding velocity axis
   - Verify length matches spectrum length
   - Convert from m/s to km/s
   - Use as x-axis for plotting
3. Falls back to channel numbering if VELOCITY_AXIS unavailable

### Backward Compatibility
✓ Fully compatible with existing FITS files without VELOCITY_AXIS
✓ Works with original unreduced data
✓ Works with newly reduced data with extracted spectra
✓ No breaking changes to CLI interface

## Example Output

### With Velocity Axis (Reduced Data)
```
Plotting 8 spectra for object 'M51CENTER'
Original FITS row indices (for reference): [2991, 4988, 5337, 6390, 11264, 14003, 20207, 26895]

✓ Plot saved to test_plot_velocity.png
```

Each subplot shows:
- **X-axis**: Velocity (km/s) from 350-700 km/s
- **Y-axis**: Intensity (flux)
- **Title**: Row index, object name, NaN fraction

### With Channel Axis (Original Data)
```
✓ Plot saved to test_plot_original.png
```

Each subplot shows:
- **X-axis**: Channel (0-1263)
- **Y-axis**: Intensity (flux)
- **Title**: Row index, object name, NaN fraction

## Quality Indicators

The title of each spectrum subplot includes:
- ✓ (green) = No NaNs in spectrum
- ~ (blue) = <10% NaNs
- ! (orange) = 10-30% NaNs
- ✗ (red) = >30% NaNs

This color-coding is independent of the axis type used.

## Testing

Verified with:
- Original M51 data (1264 channels, no VELOCITY_AXIS) → Uses channel axis ✓
- Reduced extracted data (700 channels, has VELOCITY_AXIS) → Uses velocity axis ✓
- 8 sample M51CENTER spectra plotted successfully ✓

## Related Files

- `src/oi_zeigt/cli.py` - Updated `plot_sample_spectra()` function
- `src/oi_zeigt/reduction/core.py` - Creates VELOCITY_AXIS column during reduction
- `config.toml` - Configuration for extraction/baseline reduction
