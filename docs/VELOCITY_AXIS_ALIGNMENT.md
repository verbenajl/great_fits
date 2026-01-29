# Velocity Axis Alignment in PCA Decomposition and Correction

## Overview

The key insight is that **velocity alignment between PCA components and spectra happens automatically** because they use the same FITS file and therefore the same velocity axis.

## How It Works

### 1. During Decomposition (`pca_decompose_per_mission.py`)

**Input**:
- Mission-specific telluric line parameters from `mission_id_parameters.yml` (in **km/s**)
  - Example: `telluric_line_center: 592 km/s`, `telluric_line_width: 30 km/s`

**Process**:
- Extract velocity axis from FITS file: `VELOCITY_AXIS` column (in **m/s**)
- Convert velocity parameters from **km/s to channel indices**:
  ```python
  v_min = center - width/2  # = 592 - 15 = 577 km/s
  v_max = center + width/2  # = 592 + 15 = 607 km/s
  line_mask = (velocity_axis >= v_min) & (velocity_axis <= v_max)
  line_channels = np.where(line_mask)[0]  # e.g., channels 454-513
  ```

**Output**:
- PCA components (eigenvectors) stored as **channel space** arrays (700 values per component)
- The components encode spectral patterns in **channel indices**, NOT velocities
- Telluric region already masked during decomposition (masked channels set to continuum mean)

### 2. During Correction (`pca_correct_fits.py`)

**Input**:
- PCA components (eigenvectors in channel space, 700 values each)
- Spectra to correct (also from the same FITS file, 700 channels each)

**Process**:
- The spectra being corrected come from the **same FITS file** as the decomposition reference
- This means they have the **exact same velocity axis**
- Apply correction directly in channel space:
  ```python
  corrected_spectrum = spectrum - sum(coeff[i] * components[i])
  ```
- All arithmetic is in **channel indices** (0-699)

**Result**:
- Automatic velocity alignment achieved through **shared FITS velocity axis**
- No manual velocity-to-channel conversion needed
- Telluric-masked region is automatically preserved (already zeroed in components)

## Why This Works

### Example: M51 Data (2017-02-01_GR_F367)

| Aspect | Value | Unit |
|--------|-------|------|
| Total channels | 700 | channels |
| Velocity range | 350.1 - 699.6 | km/s |
| Channel spacing | ~0.5 | km/s/channel |
| Telluric center | 592 | km/s |
| Telluric width | 30 | km/s total |
| Telluric range | 577 - 607 | km/s |
| Masked channels | 454 - 513 | channel indices |
| Masked region | 60 channels | ~30 km/s |

**Key observation**: Both decomposition and correction use the **same FITS file** and therefore the **same velocity calibration**. This ensures:

1. ✅ Telluric line is masked at the correct velocity during decomposition
2. ✅ Components are uncontaminated by telluric features
3. ✅ When applied to spectra, channel indices correspond to identical velocities
4. ✅ No additional alignment needed

## Important Notes

### For Users

- **Do not** manually convert velocity parameters to channels - the code handles this
- **Do** ensure mission parameters are in **km/s** (as documented in `mission_id_parameters.yml`)
- **Do** use spectra from the **same FITS file** for both decomposition and correction

### For Developers

- Components are stored in **channel space**, not velocity space
- Velocity information is only used during:
  - Decomposition: converting km/s masking parameters to channel indices
  - Plotting: optionally adding velocity axis labels
- The pickle files do **not** include velocity axis data (not needed for correction)
- When plotting components, the x-axis is channel indices (0-699)

## Validation

The line masking is working correctly when:
```
Masking line window 577.0-607.0 km/s for decomposition
Found 60 channels in line window: indices 454-513
Corresponding velocity range: 577.1 to 606.6 km/s
✓ Masked 60 channels in line region
```

This shows:
- ✅ Correct velocity range (577-607 km/s)
- ✅ Correct number of channels (60 = 30 km/s ÷ 0.5 km/s/channel)
- ✅ Correct channel indices mapping to correct velocities
