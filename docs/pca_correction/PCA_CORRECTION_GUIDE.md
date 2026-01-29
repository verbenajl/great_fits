# PCA-Based Spectral Correction Guide

## Overview

The `pca_correct` command applies PCA-based correction to M51CENTER (or other science target) spectra using decomposition components learned from SKYCHOPDIFF (sky difference) spectra.

## Key Features

1. **Science Line Detection**: Uses OpenCV to automatically detect emission/absorption lines
2. **Fallback Window**: If line detection fails, uses configured velocity window from `config.toml`
3. **Component-Based Correction**: Subtracts PCA components weighted by noise ratios
4. **Diagnostic Plots**: Generates before/after comparison plots per scan/telescope
5. **Selective Correction**: Can filter by object and apply cutoff thresholds

## Workflow

### Step 1: PCA Decomposition (already done with `pca_decompose`)

```bash
pca_decompose --config config.toml --plot-components
```

This produces a pickle file containing:
- PCA components (from SKYCHOPDIFF spectra)
- Explained variance ratios
- Mean spectrum

Output: `output/pca_components/decomposition_<mission_id>_<date>_components.pkl`

### Step 2: Apply PCA Correction

```bash
pca_correct \
  --input /path/to/reduced_data.fits \
  --decomposition output/pca_components/decomposition_2017-02-01_GR_F367_2017021_components.pkl \
  --output /path/to/corrected_data.fits \
  --object M51CENTER \
  --plot
```

## Science Line Detection Algorithm

The correction process includes automatic detection of science lines to protect them during fitting:

### Algorithm Steps:

1. **Gaussian Blur** (kernel_size=51 by default)
   - Smooths the spectrum data to identify broad features
   - Helps distinguish science lines from noise

2. **Threshold Detection**
   - Mean + N×σ threshold (N=2.0 by default)
   - Creates binary mask of line regions

3. **Iterative Refinement**
   - Refines threshold by re-estimating statistics
   - Excludes known line regions from background estimate
   - Converges in typically 2-3 iterations

4. **Fallback to Config Window**
   - If line detection fails or produces no results
   - Uses velocity window from `config.toml` `[reduction]` section
   - Example: `window=[450, 500]` (km/s)

### Line Detection Parameters:

```bash
pca_correct \
  --input reduced_data.fits \
  --decomposition decomp.pkl \
  --output corrected.fits \
  --line-kernel-size 51          # Gaussian blur kernel (must be odd)
  --line-cutoff-std 2.0          # Threshold in sigma
  --smoothing-kernel 11          # Optional refinement kernel
  --config-window 450 500        # Fallback window (km/s)
```

### Disabling Line Detection:

```bash
pca_correct \
  --input reduced_data.fits \
  --decomposition decomp.pkl \
  --output corrected.fits \
  --no-line-detection            # Skip line detection
  --config-window 450 500        # Use config window only
```

## Component Selection & Cutoff Thresholds

The correction process uses two types of cutoffs:

### 1. Variance Cutoff

Only use components that explain more than X% of variance:

```bash
pca_correct \
  --input reduced_data.fits \
  --decomposition decomp.pkl \
  --output corrected.fits \
  --variance-cutoff 0.01         # Use only components > 1% variance
```

### 2. Noise Ratio Cutoff

Skip components where noise_ratio > threshold:

```bash
pca_correct \
  --input reduced_data.fits \
  --decomposition decomp.pkl \
  --output corrected.fits \
  --noise-ratio-cutoff 3.0       # Skip if noise_ratio > 3.0
```

**Noise Ratio Definition:**
- Noise_ratio = std(spectrum) / std(fitted_component)
- High noise_ratio = component is weak relative to spectrum noise
- Skip high noise_ratio components to avoid noise amplification

## Usage Examples

### Basic Correction with Auto Line Detection

```bash
pca_correct \
  --input reduced_data.fits \
  --decomposition decomp.pkl \
  --output corrected.fits
```

### Correct Only M51CENTER with Plots

```bash
pca_correct \
  --input reduced_data.fits \
  --decomposition decomp.pkl \
  --output corrected.fits \
  --object M51CENTER \
  --plot
```

### Conservative Correction (Higher Cutoffs)

```bash
pca_correct \
  --input reduced_data.fits \
  --decomposition decomp.pkl \
  --output corrected.fits \
  --variance-cutoff 0.05         # Only strong components
  --noise-ratio-cutoff 2.0       # More selective
  --plot
```

### Aggressive Correction (Lower Cutoffs)

```bash
pca_correct \
  --input reduced_data.fits \
  --decomposition decomp.pkl \
  --output corrected.fits \
  --variance-cutoff 0.001        # Include weak components
  --noise-ratio-cutoff 5.0       # Less selective
  --plot
```

### Manual Line Detection Tuning

```bash
pca_correct \
  --input reduced_data.fits \
  --decomposition decomp.pkl \
  --output corrected.fits \
  --line-kernel-size 31          # Smaller kernel = finer lines
  --line-cutoff-std 3.0          # Higher threshold = fewer detections
  --plot
```

## Output

### Corrected FITS File

- Same structure as input FITS
- SPECTRUM column replaced with corrected values
- All other columns preserved (OBJECT, SCAN, TELESCOP, etc.)
- VELOCITY_AXIS unchanged

### Diagnostic Plots (with `--plot`)

For each unique (scan, telescope) combination:

```
output/pca_corrected/pca_correction_<mission_id>_<scan>_<telescope>.png
```

Contains 4 panels:
1. **Mean Spectra**: Original (red) vs Corrected (green)
2. **Original Spectra**: 2D heatmap of all spectra before correction
3. **Corrected Spectra**: 2D heatmap after correction
4. **Difference**: Original - Corrected (shows removed components)

## Configuration Integration

From `config.toml`, the following parameters are used:

```toml
[reduction]
baseline = 3
window = [450, 500]    # Fallback window for line detection (km/s)

[pca]
n_components = 5       # Number of components to use for correction
pca_source = "SKYCHOPDIFF"  # Source used for decomposition
```

## Typical Workflow

```bash
# 1. Decompose SKYCHOPDIFF spectra
pca_decompose --config config.toml --plot-components

# 2. Apply correction to M51CENTER
pca_correct \
  --input /path/to/reduced_data.fits \
  --decomposition output/pca_components/decomposition_2017-02-01_GR_F367_2017021_components.pkl \
  --output /path/to/corrected_data.fits \
  --object M51CENTER \
  --plot

# 3. Verify results
# - Check output FITS file dimensions and values
# - Review diagnostic plots in output/pca_corrected/
# - Compare original vs corrected spectra
```

## Troubleshooting

### Issue: No Science Lines Detected

**Symptom**: Line detection reports 0 channels detected

**Solutions**:
1. Increase `--line-cutoff-std` (less sensitive):
   ```bash
   pca_correct ... --line-cutoff-std 1.5
   ```

2. Decrease `--line-kernel-size` (sharper detection):
   ```bash
   pca_correct ... --line-kernel-size 31
   ```

3. Use config window instead:
   ```bash
   pca_correct ... --no-line-detection --config-window 450 500
   ```

### Issue: Over-Correction (Artifacts Introduced)

**Symptom**: Corrected spectra have unrealistic features

**Solutions**:
1. Increase component cutoffs:
   ```bash
   pca_correct ... --variance-cutoff 0.05 --noise-ratio-cutoff 2.0
   ```

2. Disable weak components:
   ```bash
   pca_correct ... --variance-cutoff 0.02
   ```

### Issue: Under-Correction (Residuals Remain)

**Symptom**: Corrected spectra still contain systematic variations

**Solutions**:
1. Decrease component cutoffs:
   ```bash
   pca_correct ... --variance-cutoff 0.001 --noise-ratio-cutoff 5.0
   ```

2. Include all components:
   ```bash
   pca_correct ... --variance-cutoff 0.0 --noise-ratio-cutoff 999
   ```

### Issue: Science Features Getting Masked

**Symptom**: Real science lines are being masked during correction

**Solutions**:
1. Adjust line detection kernel:
   ```bash
   pca_correct ... --line-kernel-size 71  # Larger = broader detection
   ```

2. Use explicit config window:
   ```bash
   pca_correct ... --config-window 450 510  # Explicit science line region
   ```

## Performance Notes

- **Line Detection**: ~100 spectra analyzed (for speed)
  - Extrapolated to all spectra in file
  - Increases runtime by ~5-10% but provides robust masking

- **Component Fitting**: O(n_components × n_channels) per spectrum
  - With 5 components and 700 channels: ~3500 operations per spectrum
  - ~14000 spectra × 3500 ops = ~50M operations total (seconds to minutes)

- **Memory Usage**: ~2× input FITS file size
  - Entire FITS file loaded into memory
  - Original + corrected arrays maintained during processing

## References

- `pca_utilities.find_lines()` - Original line detection algorithm (from pca_correct.py)
- `decompose.py` - Decomposition workflow documentation
- `config.toml` - Configuration file with mission parameters
