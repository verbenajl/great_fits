# PCA Decompose --fits Option Summary

## Changes Made (January 20, 2026)

### New Command-Line Option: `--fits`

The `pca_decompose` command now accepts a `--fits` option to specify which FITS file to analyze.

**Purpose**: Allow flexible input file selection from command line, config file, or default locations.

### File Resolution Order

The FITS input file is determined in this priority order:

1. **Command line** (`--fits /path/to/file.fits`) - HIGHEST PRIORITY
   ```bash
   pca_decompose --config config.toml --fits /custom/path.fits
   ```

2. **Config file** `[output][reduced_fits]` - MEDIUM PRIORITY
   ```toml
   [output]
   reduced_fits = "/path/to/reduced_data.fits"
   ```

3. **Config file** `[input][fits_file]` - LOWEST PRIORITY (fallback)
   ```toml
   [input]
   fits_file = "/path/to/input_data.fits"
   ```

### Updated Help Message

```bash
$ pca_decompose --help
```

Shows:
```
options:
  --fits FITS           FITS file to analyze (overrides config file [output][reduced_fits])
  --n-components N_COMPONENTS
  --pca-source PCA_SOURCE
  -v, --verbose
```

### Usage Examples

**Example 1: Use config file settings**
```bash
pca_decompose --config config.toml
```
Uses `[output][reduced_fits]` from config file.

**Example 2: Override FITS file**
```bash
pca_decompose --config config.toml --fits /path/to/custom.fits
```
Analyzes custom.fits instead of config file setting.

**Example 3: Override multiple parameters**
```bash
pca_decompose --config config.toml \
  --fits /path/to/data.fits \
  --n-components 10 \
  --pca-source M51
```

**Example 4: Analyze different files sequentially**
```bash
# Analyze reduced data (from config)
pca_decompose --config config.toml

# Analyze clean data
pca_decompose --config config.toml --fits /path/to/clean_data.fits

# Analyze with different source
pca_decompose --config config.toml --pca-source NGC253
```

### Output Messages

The command now clearly shows which file is being analyzed:

```
Loading configuration from config.toml
✓ FITS file = /path/to/reduced_data.fits (from config [output][reduced_fits])
✓ n_components = 5 (from config file)
✓ pca_source = SKYCHOPDIFF (from config file)
Analyzing file: /path/to/reduced_data.fits
Loading SKYCHOPDIFF spectra from /path/to/reduced_data.fits
...
```

### Configuration File Structure

**Recommended setup**:
```toml
[input]
fits_file = "/path/to/input_data.fits"

[pca]
n_components = 5
pca_source = "SKYCHOPDIFF"

[output]
reduced_fits = "/path/to/reduced_data.fits"
clean_fits = "/path/to/clean_data.fits"
```

The `pca_decompose` command will use `[output][reduced_fits]` by default, falling back to `[input][fits_file]` if not specified.

### Files Modified

1. **`src/oi_zeigt/pca_analysis/decompose.py`**
   - Added `--fits` argument to argparse
   - Implemented file resolution logic (3-level precedence)
   - Added file validation (checks if file exists)
   - Enhanced logging messages

2. **`docs/PCA_DECOMPOSE_CONFIG.md`**
   - Updated with complete `--fits` documentation
   - Added file resolution order explanation
   - Included comprehensive examples
   - Updated parameter precedence section

### All Command-Line Options Now Available

```bash
pca_decompose --help
```

**Options**:
- `--config CONFIG` - Path to config file (default: config.toml)
- `--fits FITS` - **NEW** FITS file to analyze (overrides config)
- `--n-components N_COMPONENTS` - Number of PCA components (default: 5)
- `--pca-source PCA_SOURCE` - Source object for PCA (default: SKYCHOPDIFF)
- `-v, --verbose` - Verbose output

### Backward Compatibility

✓ Fully backward compatible:
- Existing scripts work without `--fits` option
- Config file `[output][reduced_fits]` is optional
- Falls back to `[input][fits_file]` if needed
- Default behavior unchanged for existing configs

### Benefits

1. **Flexibility**: Can analyze different files without editing config
2. **Clarity**: Shows which file is being used and where it came from
3. **Convenience**: No need to switch config files for different analyses
4. **Validation**: Checks file existence before processing
5. **Logging**: Clear messages about parameter sources

### Testing

```bash
# Test with help
pca_decompose --help

# Test with config file
pca_decompose --config config.toml

# Test with custom file
pca_decompose --config config.toml --fits /path/to/data.fits

# Test verbose
pca_decompose --config config.toml --fits /path/to/data.fits -v
```

All changes are production-ready and fully documented.
