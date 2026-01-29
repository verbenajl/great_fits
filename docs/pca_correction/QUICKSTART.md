# PCA Correction Quick Start

## Installation ✓

The package has been updated with the new `pca_correct` command.

```bash
# Already installed via: pip install -e .
# Verify:
pca_correct --help
```

## Two-Command Workflow

### Step 1: Learn from Sky Data (SKYCHOPDIFF)

```bash
cd /home/verbena/software/oi_zeigt
pca_decompose --config config.toml --plot-components
```

**Output**: `output/pca_components/decomposition_2017-02-01_GR_F367_2017021_components.pkl`

### Step 2: Correct Science Data (M51CENTER)

```bash
pca_correct \
  --input /path/to/reduced_data.fits \
  --decomposition output/pca_components/decomposition_2017-02-01_GR_F367_2017021_components.pkl \
  --output /path/to/pca_corrected.fits \
  --object M51CENTER \
  --config-window 450 500 \
  --plot \
  --verbose
```

**Output**: 
- `pca_corrected.fits` - Corrected spectra
- `output/pca_corrected/pca_correction_*.png` - Diagnostic plots

## Key Parameters

| Parameter | Default | Purpose |
|-----------|---------|---------|
| `--input` | Required | Input FITS file |
| `--decomposition` | Required | Decomposition pickle |
| `--output` | Required | Output FITS file |
| `--object` | None | Filter by OBJECT (e.g., M51CENTER) |
| `--config-window` | None | Fallback window km/s (e.g., 450 500) |
| `--line-kernel-size` | 51 | OpenCV Gaussian kernel (odd number) |
| `--line-cutoff-std` | 2.0 | Line detection sensitivity |
| `--variance-cutoff` | None | Minimum component variance |
| `--noise-ratio-cutoff` | None | Maximum noise ratio |
| `--plot` | False | Generate diagnostic plots |

## Common Use Cases

### Basic Correction
```bash
pca_correct --input in.fits --decomposition decomp.pkl --output out.fits
```

### With Diagnostics
```bash
pca_correct --input in.fits --decomposition decomp.pkl --output out.fits \
  --plot --verbose
```

### Conservative (Few Components)
```bash
pca_correct --input in.fits --decomposition decomp.pkl --output out.fits \
  --variance-cutoff 0.05 --noise-ratio-cutoff 2.0 --object M51CENTER
```

### Aggressive (Many Components)
```bash
pca_correct --input in.fits --decomposition decomp.pkl --output out.fits \
  --variance-cutoff 0.001 --noise-ratio-cutoff 5.0 --plot
```

### Fine-Tune Line Detection
```bash
pca_correct --input in.fits --decomposition decomp.pkl --output out.fits \
  --line-kernel-size 31 --line-cutoff-std 3.0 --plot
```

## What Happens

1. **Loads** PCA components from pickle file
2. **Detects** science lines using OpenCV (Gaussian blur + threshold)
3. **Falls back** to config.toml window if detection fails
4. **For each spectrum:**
   - Protects detected science lines
   - Fits PCA components to remaining channels
   - Calculates noise ratio for each component
   - Subtracts components above noise threshold
   - Restores original science line values
5. **Generates** diagnostic plots per scan/telescope
6. **Saves** corrected spectra to FITS file

## Output Files

**Main Output:**
- `pca_corrected.fits` - Corrected spectra (same structure as input)

**Diagnostic Plots** (with `--plot`):
- `output/pca_corrected/pca_correction_<mission>_<scan>_<telescope>.png`
  - Panel 1: Mean spectra comparison
  - Panel 2: Original 2D array
  - Panel 3: Corrected 2D array
  - Panel 4: Difference (removed components)

## Verification

Check the results:

```bash
# Check FITS structure
fitsinfo pca_corrected.fits

# Verify spectra shape
python3 << 'PYTHON'
from astropy.io import fits
with fits.open('pca_corrected.fits') as hdul:
    print(f"Spectra shape: {hdul[1].data['SPECTRUM'].shape}")
    print(f"n_objects: {len(hdul[1].data)}")
    print(f"Objects: {set(hdul[1].data['OBJECT'])}")
PYTHON

# Review plots
ls -lh output/pca_corrected/
```

## Troubleshooting

**Issue: "No science lines detected"**
- Try: `--line-cutoff-std 1.5` (more sensitive)
- Or: `--no-line-detection --config-window 450 500` (use manual window)

**Issue: Over-correction (artifacts)**
- Try: `--variance-cutoff 0.05 --noise-ratio-cutoff 2.0` (fewer components)

**Issue: Under-correction (residuals)**
- Try: `--variance-cutoff 0.001 --noise-ratio-cutoff 5.0` (more components)

**Issue: Science features lost**
- Try: `--line-kernel-size 71` (detect broader features)
- Or: `--config-window 450 510` (explicit protection region)

## Documentation

For detailed information, see:
- `PCA_CORRECTION_GUIDE.md` - Complete user guide
- `PCA_CORRECTION_IMPLEMENTATION.md` - Technical details
- `PCA_COMPLETE_WORKFLOW.md` - End-to-end workflow
- `DELIVERY_SUMMARY.txt` - Implementation summary

## Next Steps

1. ✓ Run decomposition: `pca_decompose --config config.toml --plot-components`
2. ✓ Run correction: `pca_correct ... (see examples above)`
3. ✓ Review plots: `output/pca_corrected/`
4. ✓ Check output FITS: Use fitsinfo or Python
5. ✓ Adjust parameters if needed
6. ✓ Use corrected spectra for analysis

---
**Ready to go!** The `pca_correct` command is standalone and ready for use.
