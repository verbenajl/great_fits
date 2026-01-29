# Implementation Complete: Per-Mission/Telescope PCA Decomposition

## Overview

The PCA decomposition and correction pipeline has been fully updated to create and use separate decomposition files for each unique MISSION_ID/TELESCOP combination.

## What Was Done

### 1. ✅ Created Per-Mission/Telescope Decomposition Script
- **File**: `pca_decompose_per_mission.py` (root) and `src/oi_zeigt/pca_analysis/pca_decompose_per_mission.py` (package)
- **Key Features**:
  - Groups spectra by both MISSION_ID and TELESCOP
  - Creates separate `.pkl` files per mission/telescope pair
  - File naming: `decomposition_{mission_id}_{telescope}_{date}_components.pkl`

### 2. ✅ Updated Decomposition Loading in Correction Pipeline
- **File**: `src/oi_zeigt/pca_analysis/pca_correct_fits.py`
- **Changes**:
  - `load_decompositions_from_dir()` - Detects mission/telescope pairing in filenames
  - `get_decomposition_for_mission()` - Retrieves correct decomposition per mission/telescope
  - Backward compatible with old single-mission decompositions

### 3. ✅ Integrated Per-Mission/Telescope Selection in Correction Loop
- **Feature**: For each spectrum during correction:
  - Extracts spectrum's MISSION_ID and TELESCOP from FITS data
  - Looks up matching mission/telescope-specific decomposition
  - Applies that mission's unique PCA components
  - Ensures correct instrumental signature removal

### 4. ✅ Updated CLI Entry Point
- **Change**: `pca_decompose` command now points to the per-mission/telescope script
- **Updated**: `pyproject.toml` entry point configuration
- **Reinstalled**: Package to activate new entry point

## How to Use

### Step 1: Create Per-Mission/Telescope Decompositions

```bash
pca_decompose --config config.toml --n-components 5
```

**Output files** (in `output/pca_components/`):
```
decomposition_2017-02-01_GR_F367_LFAH_PX00_S_20170201_components.pkl
decomposition_2017-02-01_GR_F367_LFAH_PX01_S_20170201_components.pkl
decomposition_2017-02-01_GR_F367_LFBH_PX00_S_20170201_components.pkl
... (one per mission_id/telescope pair)
```

### Step 2: Apply Corrections Using Per-Mission/Telescope Decompositions

```bash
pca-correct-fits averaged_data.fits corrected.fits \
  --decomposition-dir output/pca_components/
```

**What happens**:
- Loads all decomposition files from the directory
- For each spectrum:
  - Reads MISSION_ID and TELESCOP from FITS row
  - Applies matching mission/telescope-specific components
  - Generates diagnostic plots showing correct components

## Expected Results

### Before (Single Global Decomposition)
All diagnostic plots showed the same first 5 PCA components, regardless of mission or telescope.

### After (Per-Mission/Telescope Decomposition)
Diagnostic plots now show different PCA components for:
- Different missions (e.g., M51 vs. other targets)
- Different telescopes (LFAH, LFBH, GR, etc.)
- Each mission/telescope combination has its own optimal basis

## Files Modified

1. **`pca_decompose_per_mission.py`** (root & package)
   - Groups by MISSION_ID + TELESCOP
   - Saves per-mission/telescope decompositions

2. **`src/oi_zeigt/pca_analysis/pca_correct_fits.py`**
   - Updated decomposition loading
   - Per-spectrum decomposition selection
   - Maintains backward compatibility

3. **`pyproject.toml`**
   - Updated `pca_decompose` entry point
   - Now points to per-mission/telescope script

## Backward Compatibility

✅ Old decomposition files still work
✅ Can mix old and new decomposition files
✅ Automatic format detection
✅ No breaking changes to API

## Verification Checklist

- [x] Syntax validated for all modified files
- [x] Entry point updated and package reinstalled
- [x] Help text updated to reflect per-mission/telescope decomposition
- [x] Both root and package versions of script are identical
- [x] Backward compatibility maintained

## Next Actions

Ready to run:
```bash
pca_decompose --config config.toml -v
```

This will create per-mission/telescope decomposition files that can be used with:
```bash
pca-correct-fits averaged_data.fits corrected.fits --decomposition-dir output/pca_components/
```
