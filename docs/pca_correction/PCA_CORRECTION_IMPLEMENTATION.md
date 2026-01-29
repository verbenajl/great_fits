# PCA Correction Implementation Summary

## What Was Built

A standalone `pca_correct` command that applies PCA-based spectral correction to M51CENTER (or other science targets) using decomposition results from SKYCHOPDIFF spectra.

## Key Features Implemented

### 1. **Science Line Detection** 
- Uses OpenCV (cv2) Gaussian blur + threshold algorithm
- Automatically identifies emission/absorption lines
- Parameters: `--line-kernel-size` (51 default), `--line-cutoff-std` (2.0 default)
- Protects detected lines during PCA component fitting

### 2. **Fallback to Config Window**
- If line detection fails, uses `reduction.window` from `config.toml`
- Supports explicit window via `--config-window V_MIN V_MAX`
- Example: `--config-window 450 500` (velocity in km/s)

### 3. **Smoothing Kernel Support**
- Optional refinement kernel via `--smoothing-kernel SIZE`
- Allows fine-tuning of line detection algorithm
- Integrated into detection pipeline

### 4. **Component Selection**
- Variance cutoff: only use components > X% of variance
- Noise ratio cutoff: skip components where noise_ratio > threshold
- Independent selective correction per spectrum

### 5. **Per-Scan/Telescope Plots**
- Diagnostic visualization with `--plot` flag
- Shows: mean spectra, 2D heatmaps before/after, difference
- Saves to `output/pca_corrected/pca_correction_<mission>_<scan>_<telescope>.png`

## Command Signature

```bash
pca_correct \
  --input FITS_FILE \
  --decomposition DECOMP_PKL \
  --output OUTPUT_FITS \
  [--object OBJECT_NAME] \
  [--variance-cutoff FLOAT] \
  [--noise-ratio-cutoff FLOAT] \
  [--line-kernel-size INT] \
  [--line-cutoff-std FLOAT] \
  [--smoothing-kernel INT] \
  [--config-window V_MIN V_MAX] \
  [--plot] \
  [--plot-dir DIR] \
  [--no-line-detection] \
  [--overwrite]
```

## File Structure

```
src/oi_zeigt/pca_analysis/
├── pca_correct_fits.py          [NEW] Core correction engine
│   ├── find_science_lines()     - OpenCV line detection
│   ├── get_velocity_window_channels() - Velocity↔channel conversion
│   ├── PCACorrector             - Main correction class
│   │   ├── fit_coefficients()   - Fit PCA to spectrum
│   │   ├── apply_correction()   - Subtract components
│   │   ├── correct_fits_file()  - Batch processing
│   │   └── _generate_plots()    - Visualization
│   └── main_correct_cli()       - CLI entry point
│
└── decompose.py                 [MODIFIED] Added pca_decompose CLI
    └── pca_decompose --plot-components generates decomp.pkl

pyproject.toml                    [MODIFIED] Added CLI entry point
├── pca_decompose = decompose:main_cli
└── pca_correct = pca_correct_fits:main_correct_cli
```

## Workflow Example

```bash
# Step 1: Decompose SKYCHOPDIFF spectra
pca_decompose --config config.toml --plot-components

# Step 2: Correct M51CENTER spectra
pca_correct \
  --input /diskB/sofiaobsdata/m51/m51_central_fits/fits_files_processed/reduced_data.fits \
  --decomposition output/pca_components/decomposition_2017-02-01_GR_F367_2017021_components.pkl \
  --output /diskB/sofiaobsdata/m51/m51_central_fits/fits_files_processed/pca_corrected.fits \
  --object M51CENTER \
  --config-window 450 500 \
  --plot
```

## Configuration Integration

The implementation respects `config.toml` settings:

```toml
[reduction]
baseline = 3
window = [450, 500]    # Fallback window for science lines (km/s)

[pca]
n_components = 5       # Components used in decomposition
pca_source = "SKYCHOPDIFF"  # Decomposition source

[output]
pcad_fits = "path/to/corrected.fits"  # Output file path
```

## Technical Details

### Science Line Detection Algorithm

1. **Input**: Sample of 100 spectra (or all if fewer)
2. **Normalize**: Scale to 0-255 uint8 range for OpenCV
3. **Gaussian Blur**: Apply kernel_size (default 51)
4. **Threshold**: Mean + cutoff_std × σ
5. **Iterate**: Refine threshold by excluding detected lines (2-3 iterations)
6. **Output**: Boolean mask (n_spectra, n_channels) of line regions

### Correction Algorithm

For each spectrum:
1. Detect/apply science line mask (protected during fitting)
2. Fit PCA components to spectrum: `coeff = components @ spectrum`
3. For each component:
   - Check variance cutoff
   - Calculate noise ratio
   - Check noise ratio cutoff
   - If accepted: subtract `coeff[i] * component[i]` from spectrum
4. Apply masks back: preserve original values in masked regions

## Usage Scenarios

### Basic Usage
```bash
pca_correct --input reduced.fits --decomposition decomp.pkl --output corrected.fits
```

### With Science Line Protection
```bash
pca_correct --input reduced.fits --decomposition decomp.pkl --output corrected.fits \
  --config-window 450 500
```

### Conservative Correction (fewer components)
```bash
pca_correct --input reduced.fits --decomposition decomp.pkl --output corrected.fits \
  --variance-cutoff 0.05 --noise-ratio-cutoff 2.0
```

### Aggressive Correction (more components)
```bash
pca_correct --input reduced.fits --decomposition decomp.pkl --output corrected.fits \
  --variance-cutoff 0.001 --noise-ratio-cutoff 5.0
```

### Fine-Tune Line Detection
```bash
pca_correct --input reduced.fits --decomposition decomp.pkl --output corrected.fits \
  --line-kernel-size 31 --line-cutoff-std 3.0 --smoothing-kernel 11
```

### Disable Line Detection
```bash
pca_correct --input reduced.fits --decomposition decomp.pkl --output corrected.fits \
  --no-line-detection --config-window 450 500
```

## Output Files

### Main Output
- **Corrected FITS**: Same structure as input, SPECTRUM column replaced

### Optional Plots (with `--plot`)
- Per scan/telescope: `pca_correction_<mission>_<scan>_<telescope>.png`
  - Panel 1: Mean spectra (original vs corrected)
  - Panel 2: Original 2D spectra array
  - Panel 3: Corrected 2D spectra array
  - Panel 4: Difference (removed components)

## Next Steps

1. **Test on M51CENTER data**:
   ```bash
   pca_correct --input reduced_data.fits \
     --decomposition output/pca_components/decomposition_2017-02-01_GR_F367_2017021_components.pkl \
     --output corrected_data.fits \
     --object M51CENTER --plot
   ```

2. **Verify results**:
   - Check output FITS file
   - Review diagnostic plots
   - Compare before/after spectra

3. **Optimize parameters**:
   - Adjust cutoffs based on results
   - Fine-tune line detection if needed
   - Save settings in `config.toml`

## Documentation

See `PCA_CORRECTION_GUIDE.md` for detailed usage guide including:
- Line detection algorithm explanation
- Cutoff threshold definitions
- Troubleshooting guide
- Performance notes
- Configuration integration details
