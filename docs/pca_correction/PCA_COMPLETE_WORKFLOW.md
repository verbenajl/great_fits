# Complete PCA Analysis Workflow for M51

## Two-Stage Process

### Stage 1: PCA Decomposition (Learning from Sky Observations)

**Input**: SKYCHOPDIFF spectra from `reduced_data.fits`
- Sky difference spectra capture receiver/instrument variations
- Consistent systematic features across all observations
- Ideal for learning correlated variations

**Process**: `pca_decompose`
```bash
pca_decompose --config config.toml --plot-components
```

**What Happens**:
1. Load SKYCHOPDIFF spectra from FITS file
2. Prepare spectra (baseline removal, masking, normalization)
3. Apply telluric masking (592±15 km/s from mission_id_parameters.yml)
4. Perform PCA decomposition (5 components by default)
5. Generate visualization plots with velocity axis
6. Save decomposition to pickle file

**Output**: 
- Pickle file: `decomposition_2017-02-01_GR_F367_2017021_components.pkl`
- Plots: `pca_components_2017-02-01_GR_F367_2017021.png`
- Metrics: 13,905 valid spectra, 12.56% total variance explained

**Key Insight**: 
- Variance reduced from 96.81% → 12.56% after telluric masking
- Confirms telluric feature properly removed
- Components now capture true instrumental variations

### Stage 2: PCA Correction (Cleaning Science Observations)

**Input**: M51CENTER spectra from `reduced_data.fits`
- Science observations contain both:
  - Instrumental variations (same as in SKYCHOPDIFF)
  - Authentic science features (M51 spectrum)

**Process**: `pca_correct`
```bash
pca_correct \
  --input reduced_data.fits \
  --decomposition decomposition_2017-02-01_GR_F367_2017021_components.pkl \
  --output pca_corrected.fits \
  --object M51CENTER \
  --config-window 450 500 \
  --plot
```

**What Happens**:
1. Load M51CENTER spectra from FITS file
2. Detect science lines using OpenCV (Gaussian blur + threshold)
3. For each spectrum:
   - Mark science lines as "protected" regions
   - Fit PCA components (only to non-protected channels)
   - Calculate noise ratio for each component
   - Subtract components above noise threshold
4. Generate diagnostic plots per scan/telescope
5. Save corrected spectra to FITS file

**Output**:
- Corrected FITS: `pca_corrected.fits`
- Diagnostic plots: `output/pca_corrected/pca_correction_*.png`
- Summary statistics: n_components_used per spectrum

**Key Insight**:
- Science lines protected during fitting
- Instrumental variations removed
- Corrected spectra ready for analysis

## Data Flow Diagram

```
reduced_data.fits
│
├─ SKYCHOPDIFF spectra
│  │
│  ├─ Load & Prepare (baseline, mask, normalize)
│  │
│  ├─ Apply Telluric Mask (592±15 km/s)
│  │  └─ Variance: 96.81% → 12.56% ✓
│  │
│  ├─ PCA Decomposition (5 components)
│  │  └─ Learn instrumental variations
│  │
│  └─ Save: decomposition_*.pkl
│     ├─ components (5 × 700)
│     ├─ explained_variance_ratio
│     └─ mean_spectrum
│
├─ M51CENTER spectra
│  │
│  ├─ Load (13,986 spectra)
│  │
│  ├─ Detect Science Lines (OpenCV)
│  │  ├─ Gaussian blur (kernel=51)
│  │  ├─ Threshold detection
│  │  └─ Iterative refinement
│  │
│  ├─ For each spectrum:
│  │  ├─ Mark science lines (protected)
│  │  ├─ Fit components (to good channels)
│  │  ├─ Calculate noise ratios
│  │  ├─ Subtract accepted components
│  │  └─ Restore science lines
│  │
│  └─ Generate Plots
│     ├─ Before/after comparison
│     ├─ 2D heatmaps
│     └─ Per scan/telescope
│
└─ Output: pca_corrected.fits
   └─ Ready for further analysis
```

## Configuration Integration

The entire workflow is controlled by `config.toml`:

```toml
[input]
# M51CENTER FITS file with reduced spectra
fits_file = "M51CENTER_central_tile_..._L1.fits"

[parameters]
object = "M51"

[reduction]
baseline = 3            # Baseline order for preparation
window = [450, 500]     # Science line window (km/s) - fallback for correction

[pca]
n_components = 5        # Number of components to learn/use
pca_source = "SKYCHOPDIFF"  # Source for decomposition

[output]
reduced_fits = "reduced_data.fits"          # Input for both stages
pcad_fits = "pca_corrected.fits"            # Final output
```

## Mission-Specific Parameters

From `mission_id_parameters.yml`:

```yaml
2017-02-01_GR_F367:
  center: 592    # Telluric line center (km/s)
  width: 30      # Telluric line width (km/s)
  # Maps to channels 454-514 (60 channels)
```

Used in **Stage 1** to mask telluric absorption during decomposition learning.

## Velocity Axis Handling

FITS file contains `VELOCITY_AXIS` column:
- Values: m/s (350100-699600 = 350.1-699.6 km/s)
- 700 channels, same for all spectra
- Used for:
  - Telluric masking in decomposition
  - Science line detection in correction
  - Plotting with velocity axis (km/s)

## Typical Parameter Settings

### For M51 (conservative correction)
```bash
pca_correct \
  --input reduced_data.fits \
  --decomposition decomp.pkl \
  --output corrected.fits \
  --object M51CENTER \
  --variance-cutoff 0.02          # Only strong components
  --noise-ratio-cutoff 2.0        # Conservative
  --line-kernel-size 51           # Standard
  --line-cutoff-std 2.0           # Standard
  --config-window 450 500         # Fallback (from config)
  --plot
```

### For M51 (aggressive correction)
```bash
pca_correct \
  --input reduced_data.fits \
  --decomposition decomp.pkl \
  --output corrected.fits \
  --object M51CENTER \
  --variance-cutoff 0.001         # Include weak components
  --noise-ratio-cutoff 5.0        # Less selective
  --line-kernel-size 31           # Sharper detection
  --line-cutoff-std 3.0           # Less sensitive
  --config-window 450 500
  --plot
```

## Quality Metrics

### Stage 1 (Decomposition) Success Indicators

✓ **Telluric Masking Working**:
- Variance explained drops 8× (96.81% → 12.56%)
- Confirms telluric feature properly removed
- Components now represent instrumental variations

✓ **Component Quality**:
- Decreasing explained variance ratio across components
- No large jumps between components
- Mean spectrum visible in first component

✓ **Spectra Prepared Successfully**:
- 13,905 / 13,986 spectra (99.4%)
- Only 81 failed (likely bad data)
- No NaN values in final arrays

### Stage 2 (Correction) Success Indicators

✓ **Science Lines Protected**:
- Detected lines match expected regions
- No loss of flux in science features
- Protected regions preserved in output

✓ **Component Subtraction Working**:
- Corrected spectra show reduced instrumental features
- Noise levels reasonable (no amplification)
- Difference plot shows removed components

✓ **Per-Spectrum Consistency**:
- Components used: typically 2-5 out of 5
- Noise ratios: reasonable distribution
- No "stuck" spectra with all components used/none used

## Next Steps for Analysis

After correction:

1. **Validate Output**
   ```bash
   # Check FITS structure
   fitsinfo pca_corrected.fits
   
   # Verify spectra shape and values
   python3 -c "from astropy.io import fits; hdul=fits.open('pca_corrected.fits'); print(hdul[1].data['SPECTRUM'].shape)"
   ```

2. **Compare Before/After**
   - Review diagnostic plots
   - Check mean spectrum reduction
   - Verify science lines preserved

3. **Downstream Analysis**
   - Line fitting on corrected spectra
   - Kinematic analysis
   - Continuum level stability

4. **Parameter Tuning** (if needed)
   - Adjust cutoffs based on comparison
   - Fine-tune line detection
   - Update `config.toml` with optimal settings

## References

- **Decomposition**: `PCA_DECOMPOSITION_GUIDE.md`
- **Correction**: `PCA_CORRECTION_GUIDE.md`
- **Implementation Details**: `PCA_CORRECTION_IMPLEMENTATION.md`
- **Code**: `src/oi_zeigt/pca_analysis/`
  - `decompose.py` - Decomposition engine
  - `pca_correct_fits.py` - Correction engine
  - `line_detection.py` - Telluric masking
  - `config.py` - Configuration management
