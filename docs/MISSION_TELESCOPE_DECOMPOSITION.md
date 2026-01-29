# Mission/Telescope Per-Decomposition Implementation

## Overview

Updated the PCA decomposition and correction pipeline to create and use separate decomposition files for each unique MISSION_ID/TELESCOP combination, rather than a single global decomposition for all data.

## Changes Made

### 1. `pca_decompose_per_mission.py` (Root Directory)

**Modified Functions:**

#### `load_spectra_by_mission(fits_file: str) -> Dict[str, Dict]`
- **Changed:** Now groups spectra by both `MISSION_ID` AND `TELESCOP` (was only by `MISSION_ID`)
- **Details:**
  - Extracts both `MISSION_ID` and `TELESCOP` columns from FITS data
  - Creates unique pairs: `(mission_id, telescope)`
  - Returns dict with keys like: `"2017-02-01_GR_F367_LFAH_PX00_S"` (mission_id_telescope combination)
  - Handles filesystem-safe telescope names: replaces `/` and spaces with `_`

#### `save_results(result: DecompositionResult, mission_id: str, telescop: str, flight_date: str) -> Path`
- **Changed:** Now accepts `telescop` parameter (new)
- **Details:**
  - Creates output filenames: `decomposition_{mission_id}_{safe_telescop}_{flight_date}_components.pkl`
  - Example: `decomposition_2017-02-01_GR_F367_LFAH_PX00_S_20170201_components.pkl`

#### Main Processing Loop (lines 342-365)
- **Changed:** Iterates over mission/telescope combinations instead of just missions
- **Details:**
  - Extracts both `mission_id` and `telescop` from data
  - Passes both to `decompose_mission_spectra()` and `save_results()`
  - Logs: `"PROCESSING X MISSION/TELESCOPE COMBINATIONS"` (was "MISSIONS")

**Output:**
- Creates one decomposition file per unique mission_id/telescope pair
- Format: `decomposition_<MISSION_ID>_<TELESCOP>_<DATE>_components.pkl`

### 2. `src/oi_zeigt/pca_analysis/pca_correct_fits.py`

**Key Changes:**

#### `load_decompositions_from_dir(decomp_dir)` - Static Method (lines 299-367)
- **Changed:** Parse telescope from filenames (new feature)
- **Details:**
  - Detects if file has mission_id/telescope pairing or just mission_id
  - Identifies telescopes by known prefixes: `LFAH`, `LFAI`, `LFBH`, `LFBI`, `GR`, `PRISM`, `HR`, `R100Q`
  - Returns dict with keys: `(mission_id, telescope)` for new format or `mission_id` for old format
  - Maintains backward compatibility with old single-mission decomposition files

#### `get_decomposition_for_mission(mission_id, telescope, mission_decompositions)` - Method (lines 395-425)
- **Changed:** Now accepts `telescope` parameter (new)
- **Details:**
  - First tries to find: `mission_decompositions[(mission_id, telescope)]` (new format)
  - Falls back to: `mission_decompositions[mission_id]` (old format, backward compatible)
  - Logs appropriate warning if neither found

#### `correct_fits_file()` - Main Correction Method
- **Added parameter:** `mission_decompositions=None` (optional)
- **Modified correction loop (lines 818-893):**
  - For each spectrum, if `mission_decompositions` provided:
    - Extracts spectrum's `MISSION_ID` and `TELESCOP` from FITS data
    - Calls `get_decomposition_for_mission(mission_id, telescope, mission_decompositions)`
    - Temporarily swaps `self.components` and `self.explained_variance_ratio` with the per-mission/telescope decomposition
    - Applies correction using mission/telescope-specific components
    - Restores original components after correction
  - Allows per-spectrum decomposition selection based on instrument metadata

#### CLI Integration (lines 1549-1575)
- **Changed:** When using per-mission decompositions, adds `get_decomposition_for_mission` method to mock PCACorrector
- **Changed:** Passes `mission_decompositions=mission_decompositions_to_use` to `correct_fits_file()`
- **Result:** Correction now uses per-mission/telescope decompositions automatically

## Workflow

### Creating Per-Mission/Telescope Decompositions

```bash
python pca_decompose_per_mission.py --config config.toml --n-components 5
```

**Output files:**
```
output/pca_components/
├── decomposition_2017-02-01_GR_F367_LFAH_PX00_S_20170201_components.pkl
├── decomposition_2017-02-01_GR_F367_LFAH_PX01_S_20170201_components.pkl
├── decomposition_2017-02-01_GR_F367_LFBH_PX00_S_20170201_components.pkl
└── ... (one per mission_id/telescope combination)
```

### Applying Per-Mission/Telescope Corrections

```bash
pca-correct-fits averaged_data.fits output_corrected.fits \
  --decomposition-dir output/pca_components/
```

**How it works:**
1. Loads all decomposition files from directory
2. For each spectrum:
   - Extracts `MISSION_ID` and `TELESCOP` from FITS row
   - Looks up matching decomposition: `(mission_id, telescope)` pair
   - Applies that mission/telescope's specific PCA components
   - Stores corrected spectrum with mission-specific corrections

## Backward Compatibility

- ✅ Old single-file decompositions still work: `--decomposition decomposition_components.pkl`
- ✅ Old per-mission-only decompositions still work (treated as `mission_id` keys)
- ✅ Can mix old and new decomposition files in same directory
- ✅ CLI automatically detects which format is being used

## Benefits

1. **Mission-Specific PCA Components:** Each mission's unique instrumental signature gets its own PCA basis
2. **Telescope-Specific Components:** Different telescopes (LFAH, LFBH, etc.) get instrument-specific corrections
3. **Improved Plot Consistency:** Diagnostic plots now show correct components for each mission/telescope combination
4. **Better Data Quality:** Avoid applying inappropriate correction components to data they weren't trained on

## Example: Impact on Plots

**Before:** All plots show same first 5 PCA components (global decomposition)
```
Plot 1 (M51, 2017-02-01, LFAH): Components 1,2,3,4,5 from global decomposition
Plot 2 (M51, 2017-02-01, LFBH): Components 1,2,3,4,5 from global decomposition  
Plot 3 (Different mission, GR): Components 1,2,3,4,5 from global decomposition
```

**After:** Each plot shows mission/telescope-specific components
```
Plot 1 (M51, 2017-02-01, LFAH): Components from decomposition_2017-02-01_LFAH_20170201_components.pkl
Plot 2 (M51, 2017-02-01, LFBH): Components from decomposition_2017-02-01_LFBH_20170201_components.pkl
Plot 3 (Different mission, GR):  Components from decomposition_MISSION2_GR_DATE_components.pkl
```

## Testing

To verify the implementation:

```bash
# 1. Create per-mission/telescope decompositions
python pca_decompose_per_mission.py --config config.toml

# 2. Check that multiple decomposition files were created
ls -la output/pca_components/

# 3. Apply corrections using the per-mission/telescope decompositions
pca-correct-fits averaged_data.fits corrected.fits --decomposition-dir output/pca_components/

# 4. Check diagnostic plots for different components per mission/telescope
# Plots should now show different first components in different plot titles
```

## Files Modified

1. `/home/verbena/software/oi_zeigt/pca_decompose_per_mission.py` - Decomposition creation
2. `/home/verbena/software/oi_zeigt/src/oi_zeigt/pca_analysis/pca_correct_fits.py` - Decomposition loading and application
