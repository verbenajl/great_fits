# Updated CLI Entry Point for Per-Mission/Telescope Decomposition

## Summary of Changes

The `pca_decompose` CLI command has been updated to use the new per-mission/telescope decomposition script instead of the older decompose.py module.

### What Changed

**File: `pyproject.toml`**
```toml
# OLD
pca_decompose = "oi_zeigt.pca_analysis.decompose:main_cli"

# NEW
pca_decompose = "oi_zeigt.pca_analysis.pca_decompose_per_mission:main_cli"
```

### New Behavior

When you run:
```bash
pca_decompose --config config.toml
```

The command now:
1. Groups spectra by **MISSION_ID and TELESCOP** (not just MISSION_ID)
2. Creates separate decomposition files per mission/telescope pair
3. Names files: `decomposition_{mission_id}_{telescope}_{date}_components.pkl`

### Files Structure

- **Root file** (kept for backward compatibility): `/home/verbena/software/oi_zeigt/pca_decompose_per_mission.py`
- **Package file** (for CLI entry point): `/home/verbena/software/oi_zeigt/src/oi_zeigt/pca_analysis/pca_decompose_per_mission.py`

Both files are identical. The package version is used for the `pca_decompose` CLI command.

### Usage Example

```bash
# Create per-mission/telescope decompositions
pca_decompose --config config.toml --n-components 5

# Verbose output
pca_decompose --config config.toml --n-components 5 -v
```

### Output Files

Decomposition files are saved to: `output/pca_components/`

Example output for M51 data:
```
output/pca_components/
├── decomposition_2017-02-01_GR_F367_LFAH_PX00_S_20170201_components.pkl
├── decomposition_2017-02-01_GR_F367_LFAH_PX01_S_20170201_components.pkl
├── decomposition_2017-02-01_GR_F367_LFBH_PX00_S_20170201_components.pkl
└── ... (one per mission_id/telescope combination)
```

### Next Step: Use With pca_correct

Once decompositions are created, apply them to your data:

```bash
pca-correct-fits averaged_data.fits corrected.fits \
  --decomposition-dir output/pca_components/
```

The correction pipeline will automatically detect and use the per-mission/telescope decompositions.
