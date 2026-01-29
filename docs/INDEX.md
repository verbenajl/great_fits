# OI_ZEIGT Documentation Index

Welcome to the comprehensive documentation for the OI_ZEIGT spectral processing pipeline. This documentation is organized by topic for easy navigation.

## 🚀 Quick Start

- **[PCA Correction Quick Start](./pca_correction/QUICKSTART.md)** - 5 minutes to your first correction
- **[Main README](./README.md)** - Project overview and status

## 📚 Documentation by Topic

### PCA Correction (Spectral Denoising)
The `pca_correct` command applies learned PCA components to remove instrumental noise from spectral observations.

- **[PCA Correction Hub](./pca_correction/README.md)** - Central documentation hub
  - Overview of science line detection
  - Velocity axis handling (m/s → km/s conversion)
  - Config window fallback mechanism
  - Component selection filters
  - Diagnostic visualization

- **[Complete User Guide](./pca_correction/PCA_CORRECTION_GUIDE.md)** - Detailed manual with examples
- **[Implementation Details](./pca_correction/PCA_CORRECTION_IMPLEMENTATION.md)** - Technical architecture
- **[Full Workflow](./pca_correction/PCA_COMPLETE_WORKFLOW.md)** - Two-stage decomposition+correction process

### Velocity Axis Organization
Understanding how velocity information is stored and processed in FITS files.

- [Velocity Axis Guide](./velocity_axis/VELOCITY_AXIS_ORGANIZATION.md) - How velocity coordinates are managed
- Velocity values in m/s from VELOCITY_AXIS column
- Automatic conversion to km/s for user-facing parameters
- Channel index mapping via nearest-neighbor lookup

### FITS File Handling
Working with FITS binary tables and spectral data.

- [Implementation Index](./implementation/IMPLEMENTATION_INDEX.md)
- [Project Status](./PROJECT_STATUS.md)
- [Repository Guide](./REPOSITORY_GUIDE.md)

### Reference Documents
Detailed technical references and investigation logs.

- [Blank Value Detection](./BLANK_VALUE_DETECTION.md)
- [PCA Decomposition Configuration](./pca_decomposition/PCA_DECOMPOSE_CONFIG.md)

## 🏗️ Project Structure

### Source Code
```
src/oi_zeigt/
├── pca_analysis/
│   ├── pca_correct_fits.py    ← Main correction engine
│   ├── pca_decompose_fits.py  ← Decomposition engine
│   └── __init__.py
├── reduction/
│   ├── core.py                ← Spectral reduction
│   └── __init__.py
├── mapping/
│   ├── gridding.py            ← Data gridding
│   └── __init__.py
├── statistics/
│   ├── quality.py             ← Quality metrics
│   └── __init__.py
└── __init__.py
```

### Configuration
- `config.toml` - Main pipeline configuration
  - Reduction parameters (window, baseline)
  - PCA settings (n_components, pca_source)
  - Velocity calibration values

### Outputs
```
output/
├── pca_components/            ← Decomposed components (pkl files)
├── reduced_spectra/           ← Corrected FITS files
└── plots/                      ← Diagnostic visualizations
```

## ⚙️ Core Commands

### 1. PCA Decomposition (Learning Phase)
Learn PCA components from sky chopping difference spectra:

```bash
pca_decompose --config config.toml --plot-components
```

Outputs: `output/pca_components/decomposition_*.pkl`

### 2. PCA Correction (Application Phase)
Apply learned components to science observations:

```bash
pca_correct \
  --input reduced_data.fits \
  --decomposition output/pca_components/decomposition_*.pkl \
  --output pca_corrected.fits \
  --object M51CENTER \
  --config-window 450 500 \
  --plot
```

Outputs: `pca_corrected.fits` + diagnostic plots

## 🔑 Key Concepts

### Velocity Axis
- **Storage**: VELOCITY_AXIS column in FITS binary table
- **Units**: Meters per second (m/s) in file
- **Processing**: Converted to km/s for user parameters
- **Mapping**: Nearest-neighbor search to find channel indices
- **Example**: config `window = [450, 500]` in km/s → channel indices via VELOCITY_AXIS

### Science Line Detection
- **Algorithm**: OpenCV Gaussian blur + threshold
- **Purpose**: Protect emission/absorption lines from correction
- **Parameters**:
  - `--line-kernel-size` (default: 51, must be odd)
  - `--line-cutoff-std` (default: 2.0σ)
- **Fallback**: Uses `reduction.window` from config.toml if detection fails

### Component Selection
- **Variance cutoff**: Filter weak components
- **Noise ratio cutoff**: Filter unreliable fits
  - Formula: noise_ratio = std(spectrum) / std(component)
  - High value = weak signal in component

### Diagnostic Plots
- **Per scan/telescope**: Organized by observational metadata
- **Content**: Original → Corrected spectra, 2D heatmaps, difference
- **Purpose**: Visual validation of correction quality

## 📖 Reading Guide

### For New Users
1. Start with [Quick Start Guide](./pca_correction/QUICKSTART.md)
2. Review [PCA Correction Hub](./pca_correction/README.md)
3. Run example command with `--verbose` to see what's happening
4. Check diagnostic plots to understand correction quality

### For Integration/Scripting
1. Read [Implementation Details](./pca_correction/PCA_CORRECTION_IMPLEMENTATION.md)
2. Check function signatures in `pca_correct_fits.py`
3. See [Complete Workflow](./pca_correction/PCA_COMPLETE_WORKFLOW.md) for data flow
4. Use `--debug` flag for detailed logging

### For Troubleshooting
1. Check [User Guide Troubleshooting](./pca_correction/PCA_CORRECTION_GUIDE.md#troubleshooting)
2. Run with `--verbose --plot` to visualize issues
3. Review diagnostic plots for pattern analysis
4. Check config.toml for parameter alignment

### For Understanding Velocity Handling
1. Read [Velocity Axis Organization](./velocity_axis/VELOCITY_AXIS_ORGANIZATION.md)
2. See how `get_velocity_window_channels()` converts km/s → channels
3. Understand VELOCITY_AXIS column role in FITS files

## 🔗 Cross-References

### Related Topics
- **[Telluric Masking](../TELLURIC_MASKS_IN_PCA.md)** - How telluric lines are handled in decomposition
- **[FITS Structure](./implementation/IMPLEMENTATION_INDEX.md)** - Binary table organization
- **[Spectral Reduction](.)** - Pipeline baseline and window parameters

### External Resources
- Astropy FITS documentation: https://docs.astropy.org/en/stable/io/fits/
- OpenCV documentation: https://docs.opencv.org/
- SciPy signal processing: https://docs.scipy.org/doc/scipy/reference/signal.html

## 📝 Document Maintenance

This documentation reflects the state of the code as of January 20, 2026.

All documents are organized as follows:
- **Quick guides**: Get started in 5-15 minutes
- **Complete guides**: Full algorithm explanation + examples
- **Implementation docs**: Technical details for developers
- **Reference docs**: Specific parameters and configurations
- **Investigation logs**: Historical decision documentation

## 🤝 Contributing

When adding new features or making changes:

1. Update relevant documentation files
2. Add examples if new parameters are introduced
3. Update this index if new documentation sections are created
4. Keep velocity axis handling documentation current

---

**Last Updated**: January 20, 2026
**Documentation Version**: 1.0
**Code Version**: PCA Correction System v1.0

