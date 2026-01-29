# PCA Correction Documentation Map

**Date**: January 20, 2026
**Status**: Production Ready ✓

## 📂 Folder Structure

```
docs/
├── pca_correction/                    # PCA Correction Documentation Hub
│   ├── README.md                      # Central navigation hub
│   ├── QUICKSTART.md                  # 5-minute quick start guide
│   ├── README_PCA_CORRECTION.md       # Comprehensive documentation index
│   ├── PCA_CORRECTION_GUIDE.md        # Complete user manual
│   ├── PCA_CORRECTION_IMPLEMENTATION.md   # Technical implementation details
│   ├── PCA_COMPLETE_WORKFLOW.md       # Two-stage workflow documentation
│   │
│   ├── reference/                     # Reference documents
│   │   ├── CHANNEL_TRIMMING_AND_PCA_PARAMETERS.md
│   │   ├── EXPLAINED_VARIANCE_INTERPRETATION.md
│   │   └── PCA_DECOMPOSE_FITS_OPTION.md
│   │
│   └── telluric_masking/              # Telluric line masking in decomposition
│       ├── TELLURIC_MASKS_IN_PCA.md
│       ├── TELLURIC_MASKS_DETAILED_ANALYSIS.md
│       ├── TELLURIC_MASKS_QUICK_REFERENCE.md
│       ├── INVESTIGATION_SUMMARY_TELLURIC_MASKS.md
│       └── ORIGINAL_VS_CURRENT_TELLURIC_MASKING.md
│
├── VELOCITY_AXIS_REFERENCE.md         # Critical: Velocity axis handling
├── INDEX.md                           # Master documentation index
│
├── velocity_axis/                     # Velocity axis organization
│   └── *.md
│
├── implementation/                    # FITS structure and implementation
│   └── *.md
│
└── ... (other documentation folders)
```

## 🎯 Reading Paths by Role

### For First-Time Users
1. Read: [`docs/pca_correction/QUICKSTART.md`](./pca_correction/QUICKSTART.md)
2. Read: [`docs/pca_correction/README.md`](./pca_correction/README.md)
3. Run example from QUICKSTART
4. Review diagnostic plots

### For Complete Understanding
1. Read: [`docs/pca_correction/README.md`](./pca_correction/README.md)
2. Read: [`docs/pca_correction/PCA_CORRECTION_GUIDE.md`](./pca_correction/PCA_CORRECTION_GUIDE.md)
3. Read: [`docs/pca_correction/PCA_COMPLETE_WORKFLOW.md`](./pca_correction/PCA_COMPLETE_WORKFLOW.md)
4. Read: [`docs/VELOCITY_AXIS_REFERENCE.md`](./VELOCITY_AXIS_REFERENCE.md) (critical!)
5. Review [`docs/pca_correction/PCA_CORRECTION_IMPLEMENTATION.md`](./pca_correction/PCA_CORRECTION_IMPLEMENTATION.md) for code structure

### For Integration/Development
1. Read: [`docs/pca_correction/PCA_CORRECTION_IMPLEMENTATION.md`](./pca_correction/PCA_CORRECTION_IMPLEMENTATION.md)
2. Read: [`docs/pca_correction/PCA_COMPLETE_WORKFLOW.md`](./pca_correction/PCA_COMPLETE_WORKFLOW.md)
3. Read: [`docs/VELOCITY_AXIS_REFERENCE.md`](./VELOCITY_AXIS_REFERENCE.md) (ESSENTIAL!)
4. Review source: `src/oi_zeigt/pca_analysis/pca_correct_fits.py`
5. Check reference docs as needed

### For Troubleshooting
1. Check: [`docs/pca_correction/PCA_CORRECTION_GUIDE.md#troubleshooting`](./pca_correction/PCA_CORRECTION_GUIDE.md)
2. Read: [`docs/VELOCITY_AXIS_REFERENCE.md`](./VELOCITY_AXIS_REFERENCE.md) for velocity issues
3. Check: [`docs/pca_correction/reference/`](./pca_correction/reference/) for specific topics

### For Understanding Telluric Masking
1. Read: [`docs/pca_correction/telluric_masking/TELLURIC_MASKS_QUICK_REFERENCE.md`](./pca_correction/telluric_masking/TELLURIC_MASKS_QUICK_REFERENCE.md)
2. Deep dive: [`docs/pca_correction/telluric_masking/TELLURIC_MASKS_IN_PCA.md`](./pca_correction/telluric_masking/TELLURIC_MASKS_IN_PCA.md)
3. Investigation details: [`docs/pca_correction/telluric_masking/INVESTIGATION_SUMMARY_TELLURIC_MASKS.md`](./pca_correction/telluric_masking/INVESTIGATION_SUMMARY_TELLURIC_MASKS.md)

## 📋 Document Summaries

### Quick Start Documents

| Document | Size | Purpose | Audience |
|----------|------|---------|----------|
| [`QUICKSTART.md`](./pca_correction/QUICKSTART.md) | 4.9 KB | 5-minute guide to run first correction | Everyone |
| [`README.md`](./pca_correction/README.md) | 6.6 KB | Navigation hub and feature overview | Everyone |

### Complete Guides

| Document | Size | Purpose | Audience |
|----------|------|---------|----------|
| [`PCA_CORRECTION_GUIDE.md`](./pca_correction/PCA_CORRECTION_GUIDE.md) | 8.5 KB | Full user manual with troubleshooting | Users |
| [`PCA_COMPLETE_WORKFLOW.md`](./pca_correction/PCA_COMPLETE_WORKFLOW.md) | 7.8 KB | Two-stage decomposition + correction | Users & Developers |
| [`README_PCA_CORRECTION.md`](./pca_correction/README_PCA_CORRECTION.md) | 8.7 KB | Comprehensive documentation index | Researchers |

### Technical References

| Document | Size | Purpose | Audience |
|----------|------|---------|----------|
| [`PCA_CORRECTION_IMPLEMENTATION.md`](./pca_correction/PCA_CORRECTION_IMPLEMENTATION.md) | 6.6 KB | Code architecture and algorithms | Developers |
| [`VELOCITY_AXIS_REFERENCE.md`](./VELOCITY_AXIS_REFERENCE.md) | TBD | **CRITICAL**: Velocity axis handling | Everyone using velocity parameters |

### Reference Documents

| Document | Size | Purpose | Audience |
|----------|------|---------|----------|
| [`reference/CHANNEL_TRIMMING_AND_PCA_PARAMETERS.md`](./pca_correction/reference/CHANNEL_TRIMMING_AND_PCA_PARAMETERS.md) | 8.3 KB | Channel handling details | Developers |
| [`reference/EXPLAINED_VARIANCE_INTERPRETATION.md`](./pca_correction/reference/EXPLAINED_VARIANCE_INTERPRETATION.md) | 4.7 KB | Understanding variance metrics | Researchers |
| [`reference/PCA_DECOMPOSE_FITS_OPTION.md`](./pca_correction/reference/PCA_DECOMPOSE_FITS_OPTION.md) | 4.4 KB | Decomposition FITS support | Developers |

### Telluric Masking (Decomposition Phase)

| Document | Size | Purpose | Audience |
|----------|------|---------|----------|
| [`telluric_masking/TELLURIC_MASKS_QUICK_REFERENCE.md`](./pca_correction/telluric_masking/TELLURIC_MASKS_QUICK_REFERENCE.md) | 4.4 KB | Quick reference for masking | Everyone |
| [`telluric_masking/TELLURIC_MASKS_IN_PCA.md`](./pca_correction/telluric_masking/TELLURIC_MASKS_IN_PCA.md) | 8.9 KB | How masking works in decomposition | Researchers |
| [`telluric_masking/TELLURIC_MASKS_DETAILED_ANALYSIS.md`](./pca_correction/telluric_masking/TELLURIC_MASKS_DETAILED_ANALYSIS.md) | 12 KB | Deep technical analysis | Developers |
| [`telluric_masking/INVESTIGATION_SUMMARY_TELLURIC_MASKS.md`](./pca_correction/telluric_masking/INVESTIGATION_SUMMARY_TELLURIC_MASKS.md) | 8.6 KB | Investigation methodology | Researchers |
| [`telluric_masking/ORIGINAL_VS_CURRENT_TELLURIC_MASKING.md`](./pca_correction/telluric_masking/ORIGINAL_VS_CURRENT_TELLURIC_MASKING.md) | 11 KB | Comparison and implementation history | Developers |

## 🔑 Critical Concepts

### ⚠️ Velocity Axis Handling (READ THIS FIRST!)

**File**: [`docs/VELOCITY_AXIS_REFERENCE.md`](./VELOCITY_AXIS_REFERENCE.md)

**Key Points**:
- X-axis is ALWAYS a velocity axis
- VELOCITY_AXIS column in FITS is in **m/s**
- User parameters (config.toml, CLI) are in **km/s**
- Conversion: `velocity_kms = velocity_ms / 1000.0`
- Channel mapping: nearest-neighbor search via argmin(abs_difference)

**Must Know**:
```python
# FITS file → Internal processing → User parameters
velocity_axis_ms = data['VELOCITY_AXIS'][0]  # m/s from FITS
velocity_axis_kms = velocity_axis_ms / 1000.0  # Convert to km/s
window_kms = [450, 500]  # User specifies in km/s
ch_min, ch_max = get_velocity_window_channels(velocity_axis_kms, window_kms)
```

### Science Line Detection
**Primary File**: [`docs/pca_correction/README.md#science-line-detection`](./pca_correction/README.md)

Algorithm: OpenCV Gaussian blur + threshold
Kernel size: 51 pixels (configurable)
Cutoff: 2.0σ (configurable)
Fallback: Uses `reduction.window` from config.toml

### Component Selection
**Primary File**: [`docs/pca_correction/PCA_CORRECTION_GUIDE.md`](./pca_correction/PCA_CORRECTION_GUIDE.md)

Variance cutoff: Filter weak components
Noise ratio cutoff: Filter unreliable fits
Per-spectrum filtering: Independent for each spectrum

## 🚀 Quick Command Reference

### Decompose SKYCHOPDIFF (learning set)
```bash
pca_decompose --config config.toml --plot-components
```
**Output**: `output/pca_components/decomposition_*.pkl`

### Correct M51CENTER (science observations)
```bash
pca_correct \
  --input reduced_data.fits \
  --decomposition output/pca_components/decomposition_*.pkl \
  --output pca_corrected.fits \
  --object M51CENTER \
  --config-window 450 500 \
  --plot
```
**Output**: `pca_corrected.fits` + diagnostic plots

## 🔗 Cross-Reference Map

### Velocity Axis Related
- Main reference: [`docs/VELOCITY_AXIS_REFERENCE.md`](./VELOCITY_AXIS_REFERENCE.md)
- Config integration: [`docs/pca_correction/README.md`](./pca_correction/README.md)
- Detailed docs: `docs/velocity_axis/`

### PCA Decomposition Related
- Telluric masking: [`docs/pca_correction/telluric_masking/`](./pca_correction/telluric_masking/)
- Config: `config.toml` [pca] section
- Source: `src/oi_zeigt/pca_analysis/pca_decompose_fits.py`

### PCA Correction Related
- Main hub: [`docs/pca_correction/README.md`](./pca_correction/README.md)
- Complete workflow: [`docs/pca_correction/PCA_COMPLETE_WORKFLOW.md`](./pca_correction/PCA_COMPLETE_WORKFLOW.md)
- Source: `src/oi_zeigt/pca_analysis/pca_correct_fits.py`

## 📊 Documentation Statistics

| Category | Files | Total Size | Status |
|----------|-------|-----------|--------|
| Quick Start | 2 | ~11.5 KB | ✓ Complete |
| Complete Guides | 3 | ~23 KB | ✓ Complete |
| Technical Refs | 2 | ~13.2 KB | ✓ Complete |
| Reference Docs | 3 | ~17.4 KB | ✓ Complete |
| Telluric Masking | 5 | ~44.9 KB | ✓ Complete |
| **Total PCA Correction** | **15** | **~109.9 KB** | **✓ READY** |

## ✅ Verification Checklist

- ✓ All documents organized in `/docs/pca_correction/`
- ✓ Quick start guides available for new users
- ✓ Complete guides for deep understanding
- ✓ Technical implementation documentation
- ✓ Velocity axis handling explicitly documented in separate reference
- ✓ Telluric masking organized in subfolder
- ✓ Clear reading paths for different user roles
- ✓ Cross-references between related documents
- ✓ Command examples and usage patterns
- ✓ Troubleshooting guides

## 🎯 Next Steps

1. Start with [`docs/pca_correction/QUICKSTART.md`](./pca_correction/QUICKSTART.md)
2. Review [`docs/VELOCITY_AXIS_REFERENCE.md`](./VELOCITY_AXIS_REFERENCE.md) (critical!)
3. Run first correction command
4. Check diagnostic plots
5. Adjust parameters as needed
6. Refer to specific guides for deep dives

---

**Navigation**: 
- Start here: [`docs/INDEX.md`](./INDEX.md)
- Quick start: [`docs/pca_correction/QUICKSTART.md`](./pca_correction/QUICKSTART.md)
- Velocity critical: [`docs/VELOCITY_AXIS_REFERENCE.md`](./VELOCITY_AXIS_REFERENCE.md)

**Last Updated**: January 20, 2026
**Maintainer**: PCA Analysis Module
**Status**: Production Ready ✓

