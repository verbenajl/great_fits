# PCA Correction Documentation

This directory contains comprehensive documentation for the `pca_correct` standalone command, which applies PCA decomposition results to correct spectra in FITS files.

## Quick Navigation

### Getting Started
- **[QUICKSTART.md](./QUICKSTART.md)** - 5-minute guide to get up and running
- **[README_PCA_CORRECTION.md](./README_PCA_CORRECTION.md)** - Documentation index and learning paths

### Detailed Guides
- **[PCA_CORRECTION_GUIDE.md](./PCA_CORRECTION_GUIDE.md)** - Complete user manual with troubleshooting
- **[PCA_CORRECTION_IMPLEMENTATION.md](./PCA_CORRECTION_IMPLEMENTATION.md)** - Technical implementation details
- **[PCA_COMPLETE_WORKFLOW.md](./PCA_COMPLETE_WORKFLOW.md)** - End-to-end workflow documentation

### Important Notes

⚠️ **Velocity Axis Handling**:
- The **x-axis is always a velocity axis** in our FITS files
- Velocity information is read from the `VELOCITY_AXIS` column in the FITS binary table
- The velocity axis is in **m/s** in the FITS file and automatically converted to **km/s** for processing
- Channel indices are determined by matching velocity values against this axis using nearest-neighbor search

### Key Implementation Details

#### Science Line Detection
- Uses OpenCV Gaussian blur + threshold detection
- Iteratively refines threshold to avoid contaminating the detection
- Kernel size configurable (default: 51 pixels)
- Cutoff threshold in standard deviations (default: 2.0σ)

#### Config Window Fallback
1. Attempts OpenCV science line detection
2. If detection fails or returns empty mask, falls back to `--config-window` parameter
3. Config window specified in km/s, automatically converted to channel indices using VELOCITY_AXIS
4. Same window applied to all spectra as fallback

#### Component Selection
- **Variance cutoff**: Only use PCA components explaining > X% of variance
- **Noise ratio cutoff**: Skip components where noise_ratio > threshold
  - noise_ratio = std(spectrum) / std(fitted_component)
  - High value indicates weak component
- Filters applied per-spectrum independently

#### Diagnostic Plots
- Generated per scan/telescope combination
- Shows: original spectra, corrected spectra, 2D heatmaps, difference
- Helps validate correction quality and identify parameter tuning needs

## CLI Usage

```bash
pca_correct \
  --input REDUCED_FITS \
  --decomposition DECOMP_PKL \
  --output OUTPUT_FITS \
  [OPTIONS]
```

## Command Reference

### Essential Parameters
```bash
--input FILE                Input FITS file with spectra
--decomposition FILE        PCA decomposition pickle file
--output FILE              Output FITS file for corrected spectra
```

### Science Line Detection
```bash
--detect-science-lines         Enable line detection (default)
--no-line-detection            Disable line detection
--line-kernel-size INT         Gaussian kernel (must be odd, default: 51)
--line-cutoff-std FLOAT        Detection threshold (default: 2.0)
--smoothing-kernel INT         Optional refinement kernel
--config-window V_MIN V_MAX    Fallback window in km/s (e.g., 450 500)
```

### Component Selection
```bash
--variance-cutoff FLOAT        Min component variance ratio
--noise-ratio-cutoff FLOAT     Max noise ratio threshold
```

### Additional Options
```bash
--object NAME                  Filter by OBJECT name
--hdu INT                     HDU index (default: 1)
--spectrum-column STR         Column name (default: SPECTRUM)
--plot                        Generate diagnostic plots
--plot-dir DIR               Plot output directory
--overwrite                   Overwrite existing output
--verbose                     Verbose logging
--debug                       Debug logging
```

## Example Workflow

### Step 1: Decompose SKYCHOPDIFF (learning set)
```bash
pca_decompose --config config.toml --plot-components
```

### Step 2: Correct M51CENTER (target observations)
```bash
pca_correct \
  --input reduced_data.fits \
  --decomposition output/pca_components/decomposition_*.pkl \
  --output pca_corrected.fits \
  --object M51CENTER \
  --config-window 450 500 \
  --variance-cutoff 0.001 \
  --noise-ratio-cutoff 2.0 \
  --plot \
  --verbose
```

### Step 3: Verify Results
- Check output FITS file structure
- Review diagnostic plots
- Validate science line protection
- Assess before/after spectra

## Architecture

### Two-Stage Process
1. **PCA Decomposition** (`pca_decompose`): Learn components from SKYCHOPDIFF
2. **PCA Correction** (`pca_correct`): Apply components to M51CENTER

### Key Functions
- `find_science_lines()` - OpenCV-based line detection
- `get_velocity_window_channels()` - Velocity→channel conversion
- `PCACorrector.correct_fits_file()` - Main batch processing engine
- `PCACorrector.apply_correction()` - Per-spectrum correction

## Data Flow

```
FITS Input
    ↓
Load Data + VELOCITY_AXIS
    ↓
Detect Science Lines (or use config window)
    ↓
Per-Spectrum Processing:
  - Create bad channel mask (NaN, science lines)
  - Fit PCA components (masked dot product)
  - Calculate noise ratios
  - Apply cutoff filters
  - Subtract components
    ↓
Generate Diagnostic Plots
    ↓
Write Corrected FITS
```

## Troubleshooting

### Line Detection Issues
- Kernel too small: increase `--line-kernel-size`
- Kernel too large: decrease `--line-kernel-size`
- Threshold issues: adjust `--line-cutoff-std`
- Use `--no-line-detection` to disable and fall back to config window

### Component Problems
- Too many components used: increase `--variance-cutoff` or `--noise-ratio-cutoff`
- No components used: decrease cutoffs
- Check diagnostic plots to visualize quality

### Output Issues
- Wrong column name: specify `--spectrum-column`
- Wrong HDU: specify `--hdu`
- File exists error: use `--overwrite`

## Configuration Files

### config.toml Integration
```toml
[reduction]
window = [450, 500]        # Used as fallback for line detection (km/s)

[pca]
n_components = 5           # From decomposition phase
pca_source = "SKYCHOPDIFF"  # Learning set
```

### FITS File Integration
```
VELOCITY_AXIS column   → Velocity values (m/s, converted to km/s)
OBJECT column          → Source name (M51CENTER, etc.)
SPECTRUM column        → Spectra to correct
SCAN, TELESCOP         → For plot organization
```

## Performance Characteristics

- **Runtime**: 2-10 minutes for ~14,000 spectra
- **Memory**: ~2× input FITS file size
- **Processing**: Per-spectrum basis with batch I/O
- **Parallelization**: Currently single-threaded, can be extended

## See Also

- [PCA Decomposition Guide](../pca_decomposition/)
- [Velocity Axis Organization](../velocity_axis/)
- [FITS File Structure](../implementation/)

