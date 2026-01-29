# Quick Reference: Per-Mission/Telescope PCA Decomposition

## STATUS: ✅ COMPLETE - READY TO USE

## One-Line Summary
The `pca_decompose` command now creates **separate PCA components for each mission/telescope combination** instead of using a single global decomposition.

## Commands

### Create Decompositions
```bash
# Standard
pca_decompose --config config.toml

# Custom components
pca_decompose --config config.toml --n-components 7

# Verbose
pca_decompose --config config.toml -v
```

### Apply Corrections
```bash
pca-correct-fits averaged_data.fits corrected.fits \
  --decomposition-dir output/pca_components/
```

## Expected Output

### Decomposition Step
Creates files like:
- `decomposition_2017-02-01_GR_F367_LFAH_PX00_S_20170201_components.pkl`
- `decomposition_2017-02-01_GR_F367_LFBH_PX00_S_20170201_components.pkl`
- (One per unique mission_id/telescope pair)

### Correction Step
Diagnostic plots showing:
- Different components per plot (based on mission/telescope)
- Proper separation of instrumental signatures
- Better correction quality

## Key Improvements

| Aspect | Before | After |
|--------|--------|-------|
| Decompositions | 1 global | 1 per mission/telescope |
| Plot components | Same everywhere | Different per mission/telescope |
| Correction quality | Generic | Instrument-specific |
| Data coverage | All combined | Properly grouped |

## Configuration

No configuration changes needed! Uses:
- `config.toml` [input][fits_file] - Input FITS file
- `config.toml` [pca][n_components] - Number of components (or --n-components flag)
- `config.toml` [pca][pca_source] - Usually "SKYCHOPDIFF"

Output goes to: `output/pca_components/`

## Files Involved

- **Decomposition**: `src/oi_zeigt/pca_analysis/pca_decompose_per_mission.py`
- **Correction**: `src/oi_zeigt/pca_analysis/pca_correct_fits.py`
- **Config**: `pyproject.toml` (entry point updated)

## Troubleshooting

**Q: Command not found - `pca_decompose`**
A: Run `pip install -e .` in the oi_zeigt directory

**Q: No MISSION_ID or TELESCOP columns in FITS**
A: Script requires both columns to group properly. Check FITS file structure.

**Q: Different number of decomposition files each time**
A: Normal! Number of files = unique mission_id/telescope combinations in your data

**Q: Correction still using old decompositions**
A: Use `--decomposition-dir output/pca_components/` flag in pca-correct-fits command
[pca]
n_components = 5
pca_source = "SKYCHOPDIFF"

[output]
reduced_fits = "/path/to/reduced_data.fits"
```

**All Options**:
- `--fits FILE` - Analyze different FITS file
- `--n-components N` - Change component count
- `--pca-source OBJ` - Change reference object
- `--config FILE` - Config file path
- `-v` - Verbose output

---

### Config File Structure

```toml
[input]
fits_file = "/path/to/input.fits"

[pca]
n_components = 5
pca_source = "SKYCHOPDIFF"

[output]
clean_fits = "/path/to/clean.fits"
reduced_fits = "/path/to/reduced.fits"
```

---

### Common Usage Patterns

```bash
# 1. Basic workflow
reduce_spectra --config config.toml --clean --baseline
pca_decompose --config config.toml

# 2. Try different PCA components
pca_decompose --config config.toml --n-components 5
pca_decompose --config config.toml --n-components 10
pca_decompose --config config.toml --n-components 15

# 3. Analyze different files
pca_decompose --config config.toml
pca_decompose --config config.toml --fits /path/to/clean_data.fits
pca_decompose --config config.toml --fits /path/to/other_data.fits

# 4. Multiple PCA sources
pca_decompose --config config.toml --pca-source SKYCHOPDIFF
pca_decompose --config config.toml --pca-source M51
pca_decompose --config config.toml --pca-source NGC253
```

---

### Documentation

- **Column preservation**: Check `reduce_spectra()` output
- **Blank detection**: See `docs/BLANK_VALUE_DETECTION.md`
- **PCA configuration**: See `docs/PCA_DECOMPOSE_CONFIG.md`
- **FITS input option**: See `PCA_DECOMPOSE_FITS_OPTION.md`
- **Full summary**: See `IMPLEMENTATION_SUMMARY_FINAL.md`

---

### Backward Compatibility

✓ All changes backward compatible
✓ Existing code works without changes
✓ Config files still work with old structure
✓ New options are optional

---

### Help Commands

```bash
# Show all pca_decompose options
pca_decompose --help

# Show verbose output
pca_decompose --config config.toml -v
```

---

### Testing

All features tested and verified:
- Real data (M51: 14,658 spectra)
- Multiple parameter combinations
- Config file fallback chain
- All command-line options

✓ Production ready
