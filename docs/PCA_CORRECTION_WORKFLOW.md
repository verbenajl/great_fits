# PCA Correction Workflow Configuration

## Overview

The `pca_correct` command now integrates seamlessly with `config.toml` to support a fully configured correction pipeline. All three main parameters (input file, object filter, and plot directory) can be read from configuration with command-line overrides.

## Configuration-Driven Workflow

### Priority Order (Highest to Lowest)

1. **Command-line arguments** - Always override config
2. **Configuration file** (`config.toml`) - Used if CLI not provided
3. **Defaults** - Fallback if neither CLI nor config specified

## Configuration Parameters

### Input File

**Config location:** `[output][reduced_fits]`

```toml
[output]
reduced_fits = "/home/diskB/sofiaobsdata/m51/m51_central_fits/fits_files_processed/reduced_data.fits"
```

**CLI override:**
```bash
pca_correct --input /path/to/reduced_data.fits --decomposition ... --output ...
```

**Behavior:**
- If `--input` is provided on CLI, use it
- Otherwise, read from `[output][reduced_fits]` in config
- Error if neither is specified

---

### Object Filter (Substring Matching)

**Config location:** `[parameters][object]`

```toml
[parameters]
# Object name substring for filtering (PCA correction will only process spectra containing this substring)
# Example: "M51" will match "M51CENTER_UNKNOWN_...", "M51EDGE_...", etc.
object="M51"
```

**CLI override:**
```bash
pca_correct --object "M51CENTER" --input ... --decomposition ... --output ...
```

**Matching Behavior:**
- Uses **substring matching** (not exact matching)
- Case-sensitive
- Example with `object="M51"`:
  - ✓ Matches: `M51CENTER_UNKNOWN_central_tile_CB_GR_01`
  - ✓ Matches: `M51CENTER_UNKNOWN_central_tile_CB_flight1`
  - ✓ Matches: `M51EDGE_something`
  - ✗ Does NOT match: `NGC253_test`

**Behavior:**
- If `--object` is provided on CLI, use it
- Otherwise, read from `[parameters][object]` in config
- If neither specified, corrects ALL objects (with warning)

---

### Plot Output Directory

**Config location:** `[output][pca_plots_dir]`

```toml
[output]
pca_plots_dir="/home/diskB/sofiaobsdata/m51/m51_central_fits/pca_plots"
```

**CLI override:**
```bash
pca_correct --plot-dir /custom/output/path --input ... --decomposition ... --output ...
```

**Behavior:**
- If `--plot-dir` is provided on CLI, use it
- Otherwise, read from `[output][pca_plots_dir]` in config
- Fallback default: `output/pca_corrected`

---

## Example Workflows

### Minimal (Using Config Defaults)

```bash
# Requires:
# - [output][reduced_fits] in config
# - [output][pcad_fits] in config
# - [parameters][object] in config (or append --object flag)
# - [output][pca_plots_dir] in config

pca_correct \
  --decomposition /path/to/decomposition.pkl \
  --output /path/to/output.fits \
  --plot
```

Output:
```
Loading configuration from config.toml
✓ Input file = /home/diskB/.../reduced_data.fits (from config [output][reduced_fits])
✓ Object filter = M51 (from config [parameters][object])
✓ Diagnostic plots saved to: /home/diskB/.../pca_plots
✓ Corrected FITS saved to: /path/to/output.fits
```

---

### Override Input File

```bash
pca_correct \
  --input /different/reduced_data.fits \
  --decomposition /path/to/decomposition.pkl \
  --output /path/to/output.fits \
  --plot
```

Output:
```
Loading configuration from config.toml
✓ Input file = /different/reduced_data.fits (from command line)
✓ Object filter = M51 (from config [parameters][object])
✓ Diagnostic plots saved to: /home/diskB/.../pca_plots
```

---

### Override Object Filter

```bash
pca_correct \
  --object "M51CENTER" \
  --decomposition /path/to/decomposition.pkl \
  --output /path/to/output.fits \
  --plot
```

Output:
```
Loading configuration from config.toml
✓ Input file = /home/diskB/.../reduced_data.fits (from config [output][reduced_fits])
✓ Object filter = M51CENTER (from command line)
✓ Diagnostic plots saved to: /home/diskB/.../pca_plots
```

Will only correct spectra containing "M51CENTER" substring.

---

### Override Plot Directory

```bash
pca_correct \
  --plot-dir /custom/plots/location \
  --decomposition /path/to/decomposition.pkl \
  --output /path/to/output.fits \
  --plot
```

Output:
```
Loading configuration from config.toml
✓ Input file = /home/diskB/.../reduced_data.fits (from config [output][reduced_fits])
✓ Object filter = M51 (from config [parameters][object])
✓ Diagnostic plots saved to: /custom/plots/location
```

---

## Data Flow

```
config.toml ([output][reduced_fits])
    ↓
Input FITS file (all objects)
    ↓
Object Filter: "M51" substring match
    ↓
Only M51* objects pass through
    ↓
PCA Correction Applied
    ↓
Corrected Spectra
    ↓
Output FITS ([output][pcad_fits])
    ↓
Diagnostic Plots ([output][pca_plots_dir])
```

---

## Important Notes

1. **Substring Matching is Case-Sensitive**
   - `"M51"` ≠ `"m51"`
   - Configure the exact case that matches your data

2. **Config File Auto-Discovery**
   - Default config file: `config.toml` (in current directory)
   - Override: `pca_correct --config /path/to/custom.toml ...`

3. **Data Safety**
   - Original input file is unchanged
   - Output file specified via `--output` or `[output][pcad_fits]`
   - Use `--overwrite` flag to overwrite existing output files

4. **Velocity Axis Handling**
   - Velocity axis is read from FITS header: `VELOCITY_AXIS` column
   - Used for accurate science line detection and plotting
   - See `docs/VELOCITY_AXIS_DEFINITIVE.md` for details

---

## Troubleshooting

### No spectra found after object filtering
```
  Filtering to 0 spectra containing 'M51'
```
**Solution:** Check that `[parameters][object]` matches your data's OBJECT column values.

### Input file not found
```
Error: Input FITS file not found: /path/to/file.fits
```
**Solution:** 
- Verify `[output][reduced_fits]` path exists
- Use `--input` to override with correct path
- Check file permissions

### Config file not found
```
Config file not found: config.toml, using defaults
```
**Solution:**
- Ensure `config.toml` exists in current directory
- Use `--config /path/to/config.toml` to specify location
- Or provide all parameters via CLI flags

---

## See Also

- `docs/PCA_CORRECTION_IMPLEMENTATION.md` - Implementation details
- `docs/VELOCITY_AXIS_DEFINITIVE.md` - Velocity axis handling
- `docs/PCA_CORRECTION_REFERENCE.md` - Detailed PCA workflow
