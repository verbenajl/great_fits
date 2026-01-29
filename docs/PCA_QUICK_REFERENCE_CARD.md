# PCA Correction Quick Reference Card

**Version**: 1.0  
**Date**: January 20, 2026  
**Status**: Production Ready ✓

---

## ⚠️ CRITICAL: VELOCITY AXIS

**THE X-AXIS IS ALWAYS A VELOCITY AXIS**

| Parameter | Units | Where | Notes |
|-----------|-------|-------|-------|
| VELOCITY_AXIS (FITS) | m/s | FITS binary table | Converted internally to km/s |
| config window | km/s | config.toml [reduction.window] | [450, 500] example |
| CLI --config-window | km/s | Command line | e.g., --config-window 450 500 |
| Internal processing | km/s | Python code | velocity_kms = velocity_ms / 1000.0 |

**Channel Mapping**: 
```python
ch_min, ch_max = get_velocity_window_channels(velocity_axis_kms, [450, 500])
```
Uses nearest-neighbor: argmin(abs(velocity_axis - target_velocity))

---

## 🚀 Two-Stage Workflow

### Stage 1: Decomposition (Learning)
**Input**: Sky chopping difference spectra (SKYCHOPDIFF)  
**Process**: Learn PCA components with telluric masking  
**Output**: `output/pca_components/decomposition_*.pkl`

```bash
pca_decompose --config config.toml --plot-components
```

### Stage 2: Correction (Application)
**Input**: Science spectra (M51CENTER) + decomposition  
**Process**: Apply PCA correction with science line protection  
**Output**: Corrected FITS + diagnostic plots

```bash
pca_correct \
  --input reduced_data.fits \
  --decomposition output/pca_components/decomposition_*.pkl \
  --output pca_corrected.fits \
  --object M51CENTER \
  --config-window 450 500 \
  --plot
```

---

## 🔧 Key Parameters

### Science Line Detection
| Parameter | Default | Range | Notes |
|-----------|---------|-------|-------|
| --line-kernel-size | 51 | odd numbers | Gaussian blur kernel, must be odd |
| --line-cutoff-std | 2.0 | 0.5-5.0 | Threshold in standard deviations |
| --smoothing-kernel | - | odd numbers | Optional refinement kernel |
| --no-line-detection | - | flag | Disable, use config window only |

### Component Selection
| Parameter | Default | Range | Notes |
|-----------|---------|-------|-------|
| --variance-cutoff | - | 0-1 | Min fraction of variance explained |
| --noise-ratio-cutoff | - | 1-10 | Max noise_ratio threshold |

### Config Fallback
| Parameter | Source | Units | Notes |
|-----------|--------|-------|-------|
| --config-window | config.toml [reduction.window] | km/s | Used if line detection fails |

---

## 📊 Data Flow

```
┌─────────────────┐
│  FITS Input     │
│ + VELOCITY_AXIS │
└────────┬────────┘
         │
         ▼
┌─────────────────────────────────────┐
│  Load Data & Get Velocity Axis      │
│  velocity_ms = VELOCITY_AXIS[0]     │
│  velocity_kms = velocity_ms / 1000  │
└────────┬────────────────────────────┘
         │
         ▼
┌──────────────────────────────────────────┐
│  Detect Science Lines                    │
│  OR use config window (fallback)         │
│  Get channel indices via velocity mapping │
└────────┬─────────────────────────────────┘
         │
         ▼
┌──────────────────────────────────────────┐
│  Per-Spectrum Processing:                │
│  1. Create bad channel mask              │
│  2. Fit PCA components (masked)          │
│  3. Calculate noise ratios               │
│  4. Apply cutoff filters                 │
│  5. Subtract components                  │
└────────┬─────────────────────────────────┘
         │
         ▼
┌─────────────────────────────┐
│  Generate Diagnostic Plots  │
│  (before/after/difference)  │
└────────┬────────────────────┘
         │
         ▼
┌─────────────────────────┐
│  Write Corrected FITS   │
└─────────────────────────┘
```

---

## 🎯 Science Line Protection

### How It Works
1. **Detect lines**: OpenCV Gaussian blur + threshold OR use config window
2. **Convert velocity**: window (km/s) → channel indices via VELOCITY_AXIS
3. **Create mask**: channels marked as "bad" during fitting
4. **Protect values**: Original spectrum values preserved in science line regions
5. **Restore after**: Corrected spectrum in bad channels = original spectrum

### Code Flow
```python
# 1. Get velocity window from config (km/s)
window_kms = [450, 500]

# 2. Convert FITS velocity (m/s) to km/s
velocity_axis_kms = velocity_axis_ms / 1000.0

# 3. Map to channel indices
ch_min, ch_max = get_velocity_window_channels(velocity_axis_kms, window_kms)

# 4. Create mask
science_line_mask = np.zeros(n_channels, dtype=bool)
science_line_mask[ch_min:ch_max+1] = True

# 5. During correction
good_channels = ~science_line_mask & ~nan_mask
corrected = spectrum.copy()
corrected[good_channels] -= fitted_component[good_channels]
# Science line regions unchanged: corrected[~good_channels] = spectrum[~good_channels]
```

---

## 📚 Documentation Quick Links

| Purpose | File | Read Time |
|---------|------|-----------|
| **START HERE** | `docs/pca_correction/QUICKSTART.md` | 5 min |
| **Critical!** | `docs/VELOCITY_AXIS_REFERENCE.md` | 10 min |
| **Overview** | `docs/pca_correction/README.md` | 10 min |
| **Complete Guide** | `docs/pca_correction/PCA_CORRECTION_GUIDE.md` | 30 min |
| **Technical** | `docs/pca_correction/PCA_CORRECTION_IMPLEMENTATION.md` | 20 min |
| **Full Workflow** | `docs/pca_correction/PCA_COMPLETE_WORKFLOW.md` | 25 min |
| **Telluric (decomp)** | `docs/pca_correction/telluric_masking/TELLURIC_MASKS_IN_PCA.md` | 15 min |
| **Reference Map** | `docs/PCA_DOCUMENTATION_MAP.md` | 10 min |

---

## ✅ Pre-Execution Checklist

Before running `pca_correct`:

- [ ] Input FITS file exists with SPECTRUM column
- [ ] Input FITS has VELOCITY_AXIS column (in m/s)
- [ ] Decomposition pickle file exists
- [ ] config.toml has [reduction.window] set (km/s)
- [ ] Output directory writable
- [ ] Enough disk space (~2× input file size)

Before running `pca_decompose`:

- [ ] SKYCHOPDIFF FITS file exists
- [ ] VELOCITY_AXIS column present (m/s)
- [ ] config.toml [pca] settings correct
- [ ] Telluric masking enabled if needed
- [ ] Output directory writable

---

## 🐛 Troubleshooting

### Wrong channels protected
**Cause**: Velocity axis issue  
**Check**: 
```python
print("Velocity axis min/max:", velocity_ms.min()/1000, velocity_ms.max()/1000)
print("Config window:", [450, 500])
ch_min, ch_max = get_velocity_window_channels(velocity_ms/1000, [450, 500])
print("Mapped channels:", ch_min, ch_max)
```

### No components used
**Cause**: Cutoffs too strict  
**Solution**: Decrease --variance-cutoff or --noise-ratio-cutoff

### Too many components used
**Cause**: Cutoffs too loose  
**Solution**: Increase --variance-cutoff or --noise-ratio-cutoff

### Line detection fails
**Solution**: 
```bash
pca_correct ... --no-line-detection --config-window 450 500
```

### Missing VELOCITY_AXIS
**Cause**: FITS file doesn't have this column  
**Solution**:
```bash
pca_correct ... --no-line-detection --config-window 450 500
```

---

## 📝 Configuration Template

### config.toml
```toml
[reduction]
baseline = 3
window = [450, 500]        # km/s - used as fallback

[pca]
n_components = 5           # Set by decomposition
pca_source = "SKYCHOPDIFF" # Learning set
```

### Command Line (full example)
```bash
pca_correct \
  --input reduced_data.fits \
  --decomposition output/pca_components/decomposition_M51CENTER.pkl \
  --output pca_corrected.fits \
  --object M51CENTER \
  --hdu 1 \
  --spectrum-column SPECTRUM \
  --line-kernel-size 51 \
  --line-cutoff-std 2.0 \
  --config-window 450 500 \
  --variance-cutoff 0.001 \
  --noise-ratio-cutoff 2.0 \
  --plot \
  --verbose
```

---

## 🔍 Verification

After running `pca_correct`:

1. **Check output FITS**:
   ```python
   from astropy.io import fits
   hdul = fits.open('pca_corrected.fits')
   print("Columns:", hdul[1].columns.names)
   print("SPECTRUM shape:", hdul[1].data['SPECTRUM'].shape)
   ```

2. **Check plots generated**:
   ```bash
   ls -lh output/pca_corrected/pca_correction_*png
   ```

3. **Visual inspection**:
   - Open diagnostic plots
   - Compare original vs corrected
   - Check science lines are protected
   - Verify realistic correction

---

## 🎓 Key Concepts

**Noise Ratio**: 
- Formula: `noise_ratio = std(spectrum) / std(fitted_component)`
- High value: weak signal in component (skip it)
- Low value: strong signal in component (use it)

**Variance Cutoff**:
- Filters components explaining < X% of variance
- Low cutoff: use all components
- High cutoff: use only strong components

**Science Line Protection**:
- Regions defined by velocity window
- Automatic channel mapping via VELOCITY_AXIS
- Original spectrum values preserved in these channels

**Diagnostic Plots**:
- Original spectrum (red) vs corrected (green)
- 2D heatmaps showing before/after
- Difference heatmap highlighting changes
- Useful for parameter tuning

---

## 🚨 Remember

1. **The x-axis is ALWAYS a velocity axis**
2. **VELOCITY_AXIS in FITS is m/s** (converted internally to km/s)
3. **User parameters are in km/s** (config.toml, CLI)
4. **Channel mapping uses nearest-neighbor search** (argmin absolute difference)
5. **Science lines are protected** from fitting and correction

---

**Print This**: Easy reference for command line  
**Bookmark**: `docs/pca_correction/README.md` for quick access  
**Remember**: Read `docs/VELOCITY_AXIS_REFERENCE.md` first!

