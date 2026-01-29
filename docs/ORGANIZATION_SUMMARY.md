# Documentation Organization Summary - January 20, 2026

## Overview

All PCA correction documentation has been organized in the `/docs` folder with a clear emphasis on velocity axis handling, which is critical for correct operation.

## What Was Done

### ✅ Folder Structure
Created a well-organized documentation hierarchy:

```
docs/
├── pca_correction/                  # Main PCA correction hub
│   ├── README.md                    # Central navigation
│   ├── QUICKSTART.md                # 5-minute quick start
│   ├── README_PCA_CORRECTION.md     # Full index
│   ├── PCA_CORRECTION_GUIDE.md      # User manual
│   ├── PCA_CORRECTION_IMPLEMENTATION.md
│   ├── PCA_COMPLETE_WORKFLOW.md
│   ├── reference/                   # Reference materials (3 docs)
│   └── telluric_masking/            # Masking docs (5 docs)
├── VELOCITY_AXIS_REFERENCE.md       # ⚠️  CRITICAL
├── PCA_DOCUMENTATION_MAP.md         # Navigation map
├── PCA_QUICK_REFERENCE_CARD.md      # Printable card
└── INDEX.md                         # Master index (updated)
```

### ✅ Files Organized

**In `/docs/pca_correction/` (6 main documents)**:
- QUICKSTART.md (5-minute guide)
- README.md (navigation hub - NEW)
- README_PCA_CORRECTION.md (full index)
- PCA_CORRECTION_GUIDE.md (user manual)
- PCA_CORRECTION_IMPLEMENTATION.md (technical details)
- PCA_COMPLETE_WORKFLOW.md (end-to-end workflow)

**In `/docs/pca_correction/reference/` (3 reference documents)**:
- CHANNEL_TRIMMING_AND_PCA_PARAMETERS.md
- EXPLAINED_VARIANCE_INTERPRETATION.md
- PCA_DECOMPOSE_FITS_OPTION.md

**In `/docs/pca_correction/telluric_masking/` (5 masking documents)**:
- TELLURIC_MASKS_IN_PCA.md
- TELLURIC_MASKS_DETAILED_ANALYSIS.md
- TELLURIC_MASKS_QUICK_REFERENCE.md
- INVESTIGATION_SUMMARY_TELLURIC_MASKS.md
- ORIGINAL_VS_CURRENT_TELLURIC_MASKING.md

**New Top-Level Files in `/docs/`**:
- VELOCITY_AXIS_REFERENCE.md ⚠️  (CRITICAL - velocity handling)
- PCA_DOCUMENTATION_MAP.md (navigation and structure map)
- PCA_QUICK_REFERENCE_CARD.md (printable quick reference)
- INDEX.md (updated with new structure)

## ⚠️ Critical: Velocity Axis Handling

A new comprehensive reference document was created: **`docs/VELOCITY_AXIS_REFERENCE.md`**

### Key Points About the X-Axis

1. **The x-axis is ALWAYS a velocity axis** in our spectral data
2. **VELOCITY_AXIS column in FITS**: meters per second (m/s)
   - Example value: 6,000,000 m/s
   - Access: `data['VELOCITY_AXIS'][0]`
3. **User parameters (config.toml, CLI)**: kilometers per second (km/s)
   - Example: `window = [450, 500]`
   - Range: typically 400-600 km/s
4. **Conversion rule**: `velocity_kms = velocity_ms / 1000.0`
5. **Channel mapping**: Uses nearest-neighbor search
   - Function: `get_velocity_window_channels(velocity_kms, window_kms)`
   - Algorithm: `argmin(abs(velocity_axis - target_velocity))`

### Code Implementation Example
```python
# From pca_correct_fits.py lines 421-423
velocity_axis = data['VELOCITY_AXIS'][0]              # m/s from FITS
velocity_axis_kms = velocity_axis / 1000.0            # Convert to km/s
logger.info(f"Using VELOCITY_AXIS: {len(velocity_axis)} channels")

# Later, when using config window
ch_range = get_velocity_window_channels(velocity_axis_kms, config_window)
```

### Science Line Protection Flow
1. Get velocity window from `config.toml [reduction.window]` (km/s)
2. Convert VELOCITY_AXIS from FITS from m/s to km/s
3. Map velocity window to channel indices via nearest-neighbor
4. Create boolean mask of protected channels
5. During fitting, treat protected channels as "bad" (don't fit)
6. After correction, restore original values in protected regions

## Documentation Organization by Audience

### 👤 New Users (15 minutes)
1. `docs/pca_correction/QUICKSTART.md` (5 min)
2. `docs/VELOCITY_AXIS_REFERENCE.md` (5 min) ⚠️  **READ THIS!**
3. Run first command (5 min)

### 👤 Scientists/Researchers (1 hour)
1. `docs/pca_correction/README.md` (10 min)
2. `docs/VELOCITY_AXIS_REFERENCE.md` (10 min) ⚠️  **ESSENTIAL**
3. `docs/pca_correction/PCA_CORRECTION_GUIDE.md` (20 min)
4. `docs/pca_correction/PCA_COMPLETE_WORKFLOW.md` (20 min)

### 👤 Developers (1.5 hours)
1. `docs/pca_correction/PCA_CORRECTION_IMPLEMENTATION.md` (20 min)
2. `docs/VELOCITY_AXIS_REFERENCE.md` (15 min) ⚠️  **ESSENTIAL**
3. `docs/pca_correction/PCA_COMPLETE_WORKFLOW.md` (25 min)
4. `src/oi_zeigt/pca_analysis/pca_correct_fits.py` (30 min)
5. `docs/pca_correction/reference/` as needed

### 🔧 Troubleshooting
1. `docs/PCA_QUICK_REFERENCE_CARD.md` (5 min)
2. `docs/pca_correction/PCA_CORRECTION_GUIDE.md#troubleshooting`
3. `docs/VELOCITY_AXIS_REFERENCE.md` for velocity-related issues

## Navigation Documents

### `docs/INDEX.md` (Master Index)
- Overview of all documentation
- Reading guides by role
- Quick command reference
- Technology inventory
- Cross-references

### `docs/PCA_DOCUMENTATION_MAP.md` (Structure Map)
- Folder structure visualization
- Document summaries table
- Reading paths by audience
- Quick start example
- Verification checklist

### `docs/PCA_QUICK_REFERENCE_CARD.md` (Printable Card)
- Critical velocity axis info
- Two-stage workflow
- Key parameters table
- Data flow diagram
- Science line protection summary
- Troubleshooting quick guide
- Configuration templates

## What Each Document Covers

| Document | Size | Purpose | Audience |
|----------|------|---------|----------|
| QUICKSTART | 4.9 KB | Get started in 5 minutes | Everyone |
| README (hub) | 6.6 KB | Navigation & overview | Everyone |
| VELOCITY_AXIS_REFERENCE ⚠️ | ~8 KB | Velocity handling (CRITICAL) | Everyone |
| CORRECTION_GUIDE | 8.5 KB | Complete user manual | Users |
| IMPLEMENTATION | 6.6 KB | Technical architecture | Developers |
| COMPLETE_WORKFLOW | 7.8 KB | Two-stage process | Users & Devs |
| Documentation Map | ~10 KB | Structure & navigation | Everyone |
| Quick Reference Card | ~12 KB | Printable quick ref | Everyone |

## Key Emphasis: Velocity Axis

The new `VELOCITY_AXIS_REFERENCE.md` file **explicitly documents**:

✓ X-axis is always velocity (not wavelength, not channel)
✓ VELOCITY_AXIS in FITS is m/s (not km/s)
✓ User parameters are km/s (not m/s)
✓ Automatic conversion: divide by 1000
✓ Channel mapping via nearest-neighbor search
✓ Science line protection using velocity window
✓ Config.toml fallback mechanism
✓ Testing and verification methods
✓ Common issues and solutions

## Implementation in Code

The `pca_correct_fits.py` module correctly handles velocity:

- **Line 421-423**: Loads VELOCITY_AXIS from FITS (m/s), converts to km/s
- **Line 123-159**: `get_velocity_window_channels()` function converts km/s window to channel indices
- **Line 455-456**: Uses velocity window to create science line mask
- **Line 419-426**: Full velocity axis loading and conversion flow

## Documentation Statistics

| Category | Count | Size | Status |
|----------|-------|------|--------|
| Main PCA docs | 6 | ~37 KB | ✓ Complete |
| Reference docs | 3 | ~17 KB | ✓ Organized |
| Telluric masking | 5 | ~45 KB | ✓ Organized |
| Navigation docs | 4 | ~30 KB | ✓ Complete |
| **TOTAL** | **18** | **~129 KB** | **✓ READY** |

## Quick Start Commands

### Step 1: Learn components (SKYCHOPDIFF)
```bash
pca_decompose --config config.toml --plot-components
```

### Step 2: Apply correction (M51CENTER)
```bash
pca_correct \
  --input reduced_data.fits \
  --decomposition output/pca_components/decomposition_M51CENTER.pkl \
  --output pca_corrected.fits \
  --object M51CENTER \
  --config-window 450 500 \
  --plot
```

## Next Steps

1. **Review**: `docs/pca_correction/README.md` (central hub)
2. **Read**: `docs/VELOCITY_AXIS_REFERENCE.md` (critical!)
3. **Test**: Run pca_correct with `--verbose --plot`
4. **Verify**: Check diagnostic plots
5. **Adjust**: Parameters based on results

## Files to Keep for Reference

Essential documents:
- `docs/pca_correction/QUICKSTART.md` - Start here
- `docs/VELOCITY_AXIS_REFERENCE.md` - Read before using velocity parameters
- `docs/pca_correction/PCA_CORRECTION_GUIDE.md` - Full reference manual
- `docs/PCA_QUICK_REFERENCE_CARD.md` - Printable quick reference

## Important Reminders

1. ⚠️  **The X-axis is ALWAYS a velocity axis**
2. ⚠️  **VELOCITY_AXIS in FITS is m/s** (not km/s)
3. ⚠️  **User parameters are km/s** (config.toml, CLI)
4. ⚠️  **Conversion rule**: velocity_kms = velocity_ms / 1000.0
5. ⚠️  **Read `VELOCITY_AXIS_REFERENCE.md` FIRST**

## Configuration Integration

```toml
# config.toml
[reduction]
window = [450, 500]        # km/s - used for science line protection

[pca]
n_components = 5           # Set by decomposition
pca_source = "SKYCHOPDIFF" # Learning source
```

## Verification Checklist

Documentation organization:
- ✓ All PCA docs in `/docs/pca_correction/`
- ✓ Reference docs organized in `reference/` subfolder
- ✓ Telluric masking docs organized in `telluric_masking/` subfolder
- ✓ Velocity axis reference at `/docs/VELOCITY_AXIS_REFERENCE.md`
- ✓ Navigation maps created (INDEX, MAP, CARD)
- ✓ Cross-references between documents
- ✓ Clear reading paths for all audiences
- ✓ Emphasis on velocity axis handling
- ✓ Examples and command reference provided

---

**Status**: ✅ Documentation organization complete and ready for use

**Last Updated**: January 20, 2026  
**Version**: 1.0  
**Maintainer**: PCA Analysis Module

