# The Velocity Axis: Definitive Reference

**Date**: January 20, 2026  
**Status**: Critical Implementation Reference  
**Emphasis**: X-axis is ALWAYS velocity axis - read from VELOCITY_AXIS header

---

## 🚨 The Fundamental Rule

### **The X-axis is ALWAYS a velocity axis**

This is not wavelength. This is not channel index. **This is velocity.**

Every spectral element in our FITS data has an associated velocity value that determines its position on the x-axis. This velocity is:
- **Read from** the `VELOCITY_AXIS` column in the FITS binary table
- **Used to determine** which channel has which velocity
- **Essential for** science line protection and velocity windows
- **Non-negotiable** for correct data interpretation

---

## 📋 The Velocity Axis: Complete Specification

### Source of Truth

**Location in FITS File**:
```
HDU 1 (Binary Table)
  Column: VELOCITY_AXIS
  Rows: N_SPECTRA (same value for all rows)
  Data Type: Float64 (or similar)
```

**How to Access**:
```python
from astropy.io import fits
hdul = fits.open('reduced_data.fits')
hdu = hdul[1]
data = hdu.data

# The velocity axis
velocity_axis_ms = data['VELOCITY_AXIS'][0]  # m/s (same for all rows)
```

### Units: METERS per SECOND

**Why meters per second?**
- FITS standard convention for spectral axes
- Compatibility with astronomical software
- High precision for small velocity differences

**Example Values**:
```
velocity_axis_ms ≈ 6,000,000 m/s (not 6,000 km/s)
velocity_axis_ms[0] ≈ 5,700,000 m/s
velocity_axis_ms[-1] ≈ 6,300,000 m/s
```

**Actual Array**:
- 1D numpy array with shape `(n_channels,)`
- One velocity value per channel
- Typically monotonic (increasing or decreasing)

---

## 🔄 The Conversion Flow

### Step 1: Read from FITS (m/s)
```python
velocity_axis_ms = data['VELOCITY_AXIS'][0]
# Units: m/s
# Example: [5700000, 5702000, 5704000, ..., 6300000]
```

### Step 2: Convert to User Units (km/s)
```python
velocity_axis_kms = velocity_axis_ms / 1000.0
# Units: km/s
# Example: [5700, 5702, 5704, ..., 6300]
```

### Step 3: Match to User Parameters (km/s)
```python
config_window = [450, 500]  # User-specified in km/s
# This comes from config.toml [reduction.window]
# NO further conversion needed - already in km/s!
```

### Step 4: Convert Velocity Window to Channel Indices
```python
ch_min, ch_max = get_velocity_window_channels(
    velocity_axis_kms,     # Already in km/s
    config_window          # Already in km/s
)
# Returns: (channel_index_min, channel_index_max)
```

---

## 🔍 Channel Mapping: The Critical Function

### Purpose
Convert a velocity range (in km/s) to a range of channel indices

### Algorithm
```python
def get_velocity_window_channels(velocity_axis, window_km_s, tolerance=1.0):
    """
    Convert velocity window to channel indices.
    
    Uses NEAREST-NEIGHBOR matching to find channels.
    """
    if velocity_axis is None or window_km_s is None:
        return None
    
    v_min, v_max = window_km_s
    
    # Find nearest channel to minimum velocity
    ch_min = np.argmin(np.abs(velocity_axis - v_min))
    
    # Find nearest channel to maximum velocity
    ch_max = np.argmin(np.abs(velocity_axis - v_max))
    
    # Ensure proper order (handles reversed velocity arrays)
    if ch_min > ch_max:
        ch_min, ch_max = ch_max, ch_min
    
    return ch_min, ch_max
```

### Key Points
1. **Nearest-neighbor matching**: Uses `argmin(abs_difference)` not interpolation
2. **Velocity range input**: Both in km/s
3. **Channel index output**: 0-based indices for array slicing
4. **Handles reversed arrays**: Works whether velocity increases or decreases
5. **Inclusive range**: `[ch_min:ch_max+1]` includes both endpoints

### Example
```python
# Velocity axis (km/s): [5700, 5702, 5704, 5706, 5708, 5710]
# User window: [5703, 5707] km/s
# 
# Find nearest to 5703: channel 1 (5702 is closer than 5704)
# Find nearest to 5707: channel 3 (5706 is closer than 5708)
#
# Result: ch_min=1, ch_max=3
# Protected channels: [1, 2, 3] ← includes both endpoints
```

---

## 📊 Real-World Flow in pca_correct_fits.py

### Load and Prepare (Lines 419-426)
```python
# Get velocity axis if available
velocity_axis = None
if 'VELOCITY_AXIS' in hdu.columns.names:
    velocity_axis = data['VELOCITY_AXIS'][0]      # Read m/s from FITS
    velocity_axis_kms = velocity_axis / 1000.0    # Convert to km/s
    logger.info(f"Using VELOCITY_AXIS: {len(velocity_axis)} channels")
else:
    velocity_axis_kms = None
```

**What this does**:
- Reads 1D velocity array from FITS column (m/s)
- Converts to km/s for user-facing operations
- Handles missing VELOCITY_AXIS gracefully

### Detect Lines or Use Config Window (Lines 450-465)
```python
# If no lines detected and config_window provided, use config window
if science_line_mask is None and config_window is not None:
    logger.info(f"Using config window: {config_window[0]:.1f}-{config_window[1]:.1f} km/s")
    if velocity_axis_kms is not None:
        ch_range = get_velocity_window_channels(velocity_axis_kms, config_window)
        if ch_range:
            science_line_mask = np.zeros((len(indices), len(velocity_axis_kms)), dtype=bool)
            ch_min, ch_max = ch_range
            science_line_mask[:, ch_min:ch_max+1] = True
```

**What this does**:
1. Check if line detection succeeded
2. If not, fall back to config window (km/s)
3. Use velocity axis to convert window to channels
4. Create boolean mask for all spectra

### Per-Spectrum Fitting (Lines 470+)
```python
# Get bad channels: science lines + NaN
bad_channels = np.isnan(spectrum) | (~good_channels)
good_channels_mask = ~bad_channels

# Fit PCA components (only on good channels)
coeffs = fit_coefficients(spectrum, good_channels_mask)

# Apply correction (preserve bad channels)
corrected = spectrum.copy()
corrected[good_channels_mask] -= fitted_components[good_channels_mask]
# Science line channels: corrected[bad_channels] = spectrum[bad_channels]
```

**What this does**:
1. Marks science line channels as "bad" (from velocity window)
2. Only fits PCA on "good" channels
3. Only corrects "good" channels
4. Preserves original values in science line regions

---

## ✅ Verification: How to Check Velocity Axis

### 1. Verify VELOCITY_AXIS Exists
```python
from astropy.io import fits
hdul = fits.open('reduced_data.fits')
hdu = hdul[1]

print('Columns:', hdu.columns.names)
if 'VELOCITY_AXIS' in hdu.columns.names:
    print("✓ VELOCITY_AXIS column found")
else:
    print("✗ VELOCITY_AXIS column missing!")
```

### 2. Check Units (Should be m/s)
```python
velocity_axis = hdu.data['VELOCITY_AXIS'][0]
print(f"First value: {velocity_axis[0]}")      # Should be ~6,000,000 (m/s)
print(f"Last value: {velocity_axis[-1]}")      # Should be ~6,000,000 ± range
print(f"Range: {velocity_axis[-1] - velocity_axis[0]}")

# If values are small (50-100), they're probably in km/s already - ERROR!
```

### 3. Test Channel Mapping
```python
from oi_zeigt.pca_analysis.pca_correct_fits import get_velocity_window_channels

velocity_kms = velocity_axis / 1000.0
window = [450, 500]  # km/s

ch_min, ch_max = get_velocity_window_channels(velocity_kms, window)
print(f"Window {window} km/s → channels {ch_min}-{ch_max}")

# Verify by looking up actual velocities
print(f"Channel {ch_min} velocity: {velocity_kms[ch_min]:.1f} km/s")
print(f"Channel {ch_max} velocity: {velocity_kms[ch_max]:.1f} km/s")
```

### 4. Verify Data Type and Shape
```python
print(f"Type: {type(velocity_axis)}")           # Should be ndarray
print(f"Dtype: {velocity_axis.dtype}")          # Should be float64
print(f"Shape: {velocity_axis.shape}")          # Should be (n_channels,)
print(f"Size: {len(velocity_axis)} channels")
```

---

## ⚠️ Common Mistakes to Avoid

### ❌ Mistake 1: Forgetting the Unit Conversion
```python
# WRONG!
velocity_kms = velocity_axis  # Missing /1000

# CORRECT
velocity_kms = velocity_axis / 1000.0
```

### ❌ Mistake 2: Not Reading from FITS at All
```python
# WRONG!
velocity_axis = np.linspace(5700, 6300, n_channels)  # Assumed values!

# CORRECT
velocity_axis = data['VELOCITY_AXIS'][0]  # Read actual values
```

### ❌ Mistake 3: Using Channel Number Instead of Velocity
```python
# WRONG!
science_line_channels = [100, 200]  # What are these numbers?

# CORRECT
science_line_window = [450, 500]  # km/s - clear and unambiguous
ch_min, ch_max = get_velocity_window_channels(velocity_kms, science_line_window)
```

### ❌ Mistake 4: Applying Conversion Twice
```python
# WRONG!
velocity_kms = velocity_axis / 1000.0
ch_range = get_velocity_window_channels(velocity_kms / 1000.0, [450, 500])

# CORRECT
velocity_kms = velocity_axis / 1000.0
ch_range = get_velocity_window_channels(velocity_kms, [450, 500])  # Already converted!
```

### ❌ Mistake 5: Assuming Monotonic Increasing
```python
# WRONG!
ch_min = np.where(velocity_kms >= v_min)[0][0]  # Fails if reversed!

# CORRECT
ch_min = np.argmin(np.abs(velocity_kms - v_min))  # Works always
```

---

## 🔗 Integration Checklist

Before using any velocity-dependent feature:

- [ ] **Read VELOCITY_AXIS from FITS column** (not assumed)
- [ ] **Check units are m/s** (not km/s or channels)
- [ ] **Convert to km/s for processing** (`/ 1000.0`)
- [ ] **Use converted array for all user-facing operations**
- [ ] **Science line window in km/s** (from config.toml)
- [ ] **Channel mapping uses nearest-neighbor** (not interpolation)
- [ ] **Verify channel range is reasonable** (not 0-10 or 1000+)
- [ ] **Test with actual FITS file** (not synthetic data)
- [ ] **Check diagnostic plots** (verify protection correct)
- [ ] **Document assumptions in code** (why this velocity axis?)

---

## 📝 Code Pattern to Follow

### Template for Velocity Axis Handling
```python
# 1. Read from FITS
if 'VELOCITY_AXIS' in hdu.columns.names:
    velocity_axis_ms = data['VELOCITY_AXIS'][0]
    velocity_axis_kms = velocity_axis_ms / 1000.0
    logger.info(f"Loaded VELOCITY_AXIS: {len(velocity_axis_ms)} channels, "
                f"range {velocity_axis_kms[0]:.0f}-{velocity_axis_kms[-1]:.0f} km/s")
else:
    logger.error("VELOCITY_AXIS column not found in FITS file!")
    return None

# 2. If user provides velocity window
if user_velocity_window is not None:
    # user_velocity_window already in km/s (from config.toml or CLI)
    ch_min, ch_max = get_velocity_window_channels(velocity_axis_kms, user_velocity_window)
    
    if ch_min is None:
        logger.warning(f"Could not map window {user_velocity_window} to channels")
        return None
    
    logger.info(f"Velocity window {user_velocity_window} km/s → "
                f"channels {ch_min}-{ch_max} "
                f"({velocity_axis_kms[ch_min]:.1f}-{velocity_axis_kms[ch_max]:.1f} km/s)")

# 3. Create mask or apply protection
if ch_min is not None and ch_max is not None:
    mask = np.zeros(n_channels, dtype=bool)
    mask[ch_min:ch_max+1] = True
    logger.info(f"Protected {np.sum(mask)} channels")
```

---

## 🎯 Summary: The Golden Rules

1. **✓ X-axis is velocity** - Not wavelength, not channels
2. **✓ Read from FITS** - Don't assume or calculate
3. **✓ FITS is m/s** - Divide by 1000 for km/s
4. **✓ User parameters are km/s** - No conversion needed
5. **✓ Channel mapping is nearest-neighbor** - `argmin(abs_diff)`
6. **✓ Science lines use velocity window** - From config or detection
7. **✓ Verify with actual data** - Test on real FITS files
8. **✓ Document assumptions** - Make velocity handling clear

---

## 📚 Related Documentation

- **PCA Correction Hub**: `/docs/pca_correction/README.md`
- **Science Line Protection**: `/docs/pca_correction/PCA_COMPLETE_WORKFLOW.md`
- **Implementation Details**: `/docs/pca_correction/PCA_CORRECTION_IMPLEMENTATION.md`
- **Quick Reference**: `/docs/PCA_QUICK_REFERENCE_CARD.md`
- **Source Code**: `src/oi_zeigt/pca_analysis/pca_correct_fits.py`

---

**Version**: 1.0  
**Last Updated**: January 20, 2026  
**Status**: Production Reference ✓

**Remember**: The velocity axis is the foundation of spectral interpretation. Get this right, and everything else follows.

