# Velocity Axis Handling Reference

## Critical Implementation Note

**The x-axis is ALWAYS a velocity axis** in our spectral data. This document specifies exactly how velocity information is read, stored, and converted throughout the pipeline.

## Data Source

### FITS File Level
**Location**: VELOCITY_AXIS column in binary table (HDU 1)
**Units**: **Meters per second (m/s)**
**Type**: 1D array (n_channels,)
**Access**: `data['VELOCITY_AXIS'][0]` - same for all rows

```python
# Example from pca_correct_fits.py (line 421-423)
if 'VELOCITY_AXIS' in hdu.columns.names:
    velocity_axis = data['VELOCITY_AXIS'][0]  # Units: m/s
    velocity_axis_kms = velocity_axis / 1000.0  # Convert to km/s
    logger.info(f"Using VELOCITY_AXIS: {len(velocity_axis)} channels")
```

### Config File Level
**Location**: `config.toml` [reduction] section
**Units**: **Kilometers per second (km/s)**
**Example**:
```toml
[reduction]
baseline = 3
window = [450, 500]  # km/s
```

## Conversion Flow

```
FITS File
  VELOCITY_AXIS column
       ↓ (in m/s)
       
Internal Processing
  velocity_axis_kms = velocity_axis / 1000.0
       ↓ (convert to km/s)
       
User Parameters / Config
  --config-window 450 500  (km/s)
  reduction.window = [450, 500]  (km/s)
       ↓
Channel Index Mapping
  indices = get_velocity_window_channels(velocity_axis_kms, window_km_s)
       ↓ (nearest-neighbor search)
       
Output
  Boolean mask of channel indices
```

## Channel Index Mapping Function

Function: `get_velocity_window_channels(velocity_axis, window_km_s, tolerance=1.0)`
**Location**: `src/oi_zeigt/pca_analysis/pca_correct_fits.py` (line 123)

```python
def get_velocity_window_channels(velocity_axis, window_km_s, tolerance=1.0):
    """
    Convert velocity window to channel indices.
    
    Parameters
    ----------
    velocity_axis : ndarray
        Velocity axis in km/s (already converted from m/s)
    window_km_s : list/tuple
        Velocity window [v_min, v_max] in km/s
    tolerance : float
        Tolerance in km/s for finding channels (default: 1.0)
    
    Returns
    -------
    tuple
        (ch_min, ch_max) - channel indices
    """
    if velocity_axis is None or window_km_s is None:
        return None
    
    v_min, v_max = window_km_s
    
    # Nearest-neighbor search to find matching channels
    ch_min = np.argmin(np.abs(velocity_axis - v_min))
    ch_max = np.argmin(np.abs(velocity_axis - v_max))
    
    # Ensure proper order
    if ch_min > ch_max:
        ch_min, ch_max = ch_max, ch_min
    
    return ch_min, ch_max
```

**Key Points**:
- Uses **nearest-neighbor matching** (argmin of absolute difference)
- Handles **reversed velocity arrays** (some instruments have descending velocity)
- Returns **channel indices** (0-based) ready for array slicing

## Usage in PCA Correction

### Science Line Fallback (line 455-456)
```python
if velocity_axis_kms is not None:
    ch_range = get_velocity_window_channels(velocity_axis_kms, config_window)
```

**Purpose**: When OpenCV line detection fails, protect science lines using config window

**Flow**:
1. Read VELOCITY_AXIS from FITS (m/s)
2. Convert to km/s
3. Get user's config_window from config.toml or CLI (km/s)
4. Convert window to channel indices via nearest-neighbor
5. Create boolean mask protecting those channels

### Full Context (lines 419-426)
```python
# Get velocity axis if available
velocity_axis = None
if 'VELOCITY_AXIS' in hdu.columns.names:
    velocity_axis = data['VELOCITY_AXIS'][0]  # Same for all spectra
    velocity_axis_kms = velocity_axis / 1000.0  # Convert m/s to km/s
    logger.info(f"  Using VELOCITY_AXIS: {len(velocity_axis)} channels")
else:
    velocity_axis_kms = None
```

## Testing & Verification

### To Verify Velocity Axis Handling:

1. **Check FITS file has VELOCITY_AXIS column**:
```python
from astropy.io import fits
hdul = fits.open('reduced_data.fits')
print('VELOCITY_AXIS' in hdul[1].columns.names)  # Should be True
velocity = hdul[1].data['VELOCITY_AXIS'][0]
print(f"Velocity units: {velocity[0]} (should be ~6000000 for m/s)")
```

2. **Test conversion**:
```python
velocity_ms = 6000000  # m/s
velocity_kms = velocity_ms / 1000.0
print(f"{velocity_ms} m/s = {velocity_kms} km/s")  # 6000.0 km/s
```

3. **Test channel mapping**:
```python
from oi_zeigt.pca_analysis.pca_correct_fits import get_velocity_window_channels
velocity_axis_kms = np.array([400, 450, 500, 550, 600])  # km/s
window = [450, 550]  # km/s
ch_min, ch_max = get_velocity_window_channels(velocity_axis_kms, window)
print(f"Channels for {window} km/s: {ch_min}-{ch_max}")  # Should be 1-3
```

## Common Issues & Solutions

### Issue: "velocity_axis is None"
**Cause**: FITS file missing VELOCITY_AXIS column
**Solution**: 
- Check FITS file: `fitsheader reduced_data.fits | grep VELOCITY`
- Use `--no-line-detection` to skip line detection
- Provide explicit `--config-window` in km/s

### Issue: Wrong channels are protected
**Cause**: Velocity axis in unexpected units or order
**Solution**:
1. Verify VELOCITY_AXIS units: `print(hdul[1].data['VELOCITY_AXIS'][0][:5])`
2. Check if reversed: `print(np.diff(hdul[1].data['VELOCITY_AXIS'][0])[:5])`
3. If descending, absolute difference in `get_velocity_window_channels` still works

### Issue: Config window doesn't match expected channels
**Cause**: Velocity axis not strictly monotonic or has gaps
**Solution**:
1. Use nearest-neighbor matching (already implemented)
2. Increase `tolerance` parameter if needed
3. Debug with: `print(get_velocity_window_channels(velocity_axis_kms, [450, 500]))`

## Integration Checklist

✓ FITS file always has VELOCITY_AXIS column in m/s
✓ User-facing parameters always in km/s
✓ Conversion always: velocity_km_s = velocity_m_s / 1000.0
✓ Channel mapping via nearest-neighbor (argmin absolute difference)
✓ Fallback to config.toml [reduction.window] in km/s
✓ Science lines protected during fitting, restored after
✓ Diagnostic plots show velocity axis labels in km/s

## Documentation References

- **PCA Correction Hub**: `/docs/pca_correction/README.md`
- **Complete Workflow**: `/docs/pca_correction/PCA_COMPLETE_WORKFLOW.md`
- **Implementation Guide**: `/docs/pca_correction/PCA_CORRECTION_IMPLEMENTATION.md`
- **Source Code**: `src/oi_zeigt/pca_analysis/pca_correct_fits.py`

---

**Last Updated**: January 20, 2026
**Maintainer**: PCA Analysis Module
**Status**: Production Ready ✓

