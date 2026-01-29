# PCA Decomposition Configuration Guide

## Overview

The `pca_decompose` command now supports configurable parameters for PCA decomposition:

- **`--fits`**: FITS file to analyze (default: `[output][reduced_fits]` from config)
- **`n_components`**: Number of PCA components to extract (default: 5)
- **`pca_source`**: Object name to use for PCA reference spectra (default: SKYCHOPDIFF)

These parameters can be specified in three ways:
1. **Command line**: Highest priority, overrides everything
2. **Config file**: `[pca]` section and `[output][reduced_fits]` in config.toml
3. **Defaults**: Built-in defaults if neither above is specified

## Configuration File

Add the `[pca]` section to your `config.toml`:

```toml
[pca]
n_components = 5
pca_source = "SKYCHOPDIFF"

[output]
reduced_fits = "/path/to/reduced_data.fits"
```

All parameters are **optional**. If not specified, the defaults are used.

### FITS File Resolution Order

The FITS file is determined in this priority order:

1. **Command line** (`--fits` option) - highest priority
   ```bash
   pca_decompose --fits /custom/path.fits
   ```

2. **Config file** `[output][reduced_fits]`
   ```toml
   [output]
   reduced_fits = "/path/to/reduced_data.fits"
   ```

3. **Config file** `[input][fits_file]` (fallback)
   ```toml
   [input]
   fits_file = "/path/to/input_data.fits"
   ```

## Command Line Usage

### Basic Usage (uses config file values or defaults)
```bash
pca_decompose --config config.toml
```

Output:
```
Loading configuration from config.toml
✓ FITS file = /path/to/reduced_data.fits (from config [output][reduced_fits])
✓ n_components = 5 (from config file)
✓ pca_source = SKYCHOPDIFF (from config file)
Analyzing file: /path/to/reduced_data.fits
Loading SKYCHOPDIFF spectra from /path/to/reduced_data.fits
...
```

### Override FITS File
```bash
pca_decompose --config config.toml --fits /path/to/custom.fits
```

Output:
```
✓ FITS file = /path/to/custom.fits (from command line)
Analyzing file: /path/to/custom.fits
```

### Override n_components Only
```bash
pca_decompose --config config.toml --n-components 10
```

Output:
```
✓ FITS file = /path/to/reduced_data.fits (from config [output][reduced_fits])
✓ n_components = 10 (from command line)
✓ pca_source = SKYCHOPDIFF (from config file)
```

### Override pca_source Only
```bash
pca_decompose --config config.toml --pca-source M51
```

Output:
```
✓ FITS file = /path/to/reduced_data.fits (from config [output][reduced_fits])
✓ n_components = 5 (from config file)
✓ pca_source = M51 (from command line)
```

### Override Multiple Parameters
```bash
pca_decompose --config config.toml \
  --fits /custom/data.fits \
  --n-components 15 \
  --pca-source M51
```

Output:
```
✓ FITS file = /custom/data.fits (from command line)
✓ n_components = 15 (from command line)
✓ pca_source = M51 (from command line)
Analyzing file: /custom/data.fits
```

### Use All Defaults (no config sections needed)
```bash
pca_decompose --config config.toml
```
(If `[pca]` and `[output]` sections are missing from config.toml)

Output:
```
✓ FITS file = /path/from/input/section (from config [input][fits_file])
✓ n_components = 5 (default)
✓ pca_source = SKYCHOPDIFF (default)
```

### Verbose Output
Add `-v` or `--verbose` flag for detailed logging:

```bash
pca_decompose --config config.toml -v
```

## Parameter Precedence

### FITS File Precedence

```
1. Command line --fits           (highest priority)
   └─ pca_decompose --fits /path/to/file.fits

2. Config [output][reduced_fits]
   └─ reduced_fits = "/path/to/file.fits"

3. Config [input][fits_file]     (lowest priority)
   └─ fits_file = "/path/to/file.fits"
```

### n_components Precedence

```
1. Command line --n-components   (highest priority)
   └─ pca_decompose --n-components 10

2. Config [pca][n_components]
   └─ [pca]
      n_components = 5

3. Built-in default              (lowest priority)
   └─ n_components = 5
```

### pca_source Precedence

```
1. Command line --pca-source     (highest priority)
   └─ pca_decompose --pca-source M51

2. Config [pca][pca_source]
   └─ [pca]
      pca_source = "SKYCHOPDIFF"

3. Built-in default              (lowest priority)
   └─ pca_source = "SKYCHOPDIFF"
```

## Configuration Examples

### Example 1: Default Setup (Reduced Data)
```bash
# Uses reduced data from [output][reduced_fits]
pca_decompose --config config.toml
```

**config.toml**:
```toml
[pca]
n_components = 5
pca_source = "SKYCHOPDIFF"

[output]
reduced_fits = "/home/data/m51/reduced_data.fits"
```

### Example 2: Custom FITS File
```bash
# Analyze a different FITS file
pca_decompose --config config.toml --fits /custom/data.fits
```

**config.toml** (normal settings):
```toml
[pca]
n_components = 5
pca_source = "SKYCHOPDIFF"

[output]
reduced_fits = "/home/data/m51/reduced_data.fits"  # Will be overridden
```

### Example 3: Multiple Different Runs
Analyze the same config file with different settings:

```bash
# Analyze reduced_data with 5 components
pca_decompose --config config.toml

# Analyze reduced_data with 10 components
pca_decompose --config config.toml --n-components 10

# Analyze clean_data with different source
pca_decompose --config config.toml --fits /path/to/clean_data.fits --pca-source M51

# Analyze clean_data with 8 components from different source
pca_decompose --config config.toml \
  --fits /path/to/clean_data.fits \
  --n-components 8 \
  --pca-source NGC253
```

## Finding Available Source Objects

If you're not sure what source objects are available in your FITS file:

```bash
pca_decompose --config config.toml --pca-source UNKNOWN
```

The command will show you all available objects:

```
Available objects:
  SKYCHOPDIFF: 143 spectra
  M51: 8234 spectra
  NGC253: 1256 spectra
```

## Output Files

The decomposition results are saved to:
```
output/pca_components/decomposition_<MISSION_ID>_<DATE>_components.pkl
```

Example:
```
output/pca_components/decomposition_04_0116_20170201_components.pkl
```

The output filename includes:
- `MISSION_ID`: From FITS file data
- `DATE`: Flight date from DATE-OBS header

## Example Workflow

```bash
#!/bin/bash

# Step 1: Setup config with all parameters
cat > config.toml << EOF
[input]
fits_file = "/path/to/input_data.fits"

[pca]
n_components = 5
pca_source = "SKYCHOPDIFF"

[output]
clean_fits = "/path/to/clean_data.fits"
reduced_fits = "/path/to/reduced_data.fits"
datacube = "/path/to/output/datacube.fits"
EOF

# Step 2: Run PCA on reduced data (from config)
pca_decompose --config config.toml

# Step 3: Try with more components
pca_decompose --config config.toml --n-components 10

# Step 4: Try with different source
pca_decompose --config config.toml --pca-source M51

# Step 5: Analyze a different file entirely
pca_decompose --config config.toml --fits /path/to/clean_data.fits
```

## Complete Config File Example

```toml
[input]
fits_file = "/home/diskB/sofiaobsdata/m51/m51_central_region/kalibrated_all/m51_central_tile.fits"

[parameters]
object = "M51"

[filters]
blank_fraction = 0.20

[reduction]
baseline = 3
window = [450, 500]
extract = [350, 700]
group_by = "OBJECT"

[pca]
n_components = 5
pca_source = "SKYCHOPDIFF"

[output]
clean_fits = "/path/to/clean_data.fits"
reduced_fits = "/path/to/reduced_data.fits"
datacube = "/path/to/datacube.fits"

[gridding]
method = "cygrid"
beamsize_arcsec = 14.1
```

## Logging Output

The command prints which parameters are being used and where they come from:

```
Loading configuration from config.toml
✓ FITS file = /path/to/reduced_data.fits (from config [output][reduced_fits])
✓ n_components = 5 (from config file)
✓ pca_source = SKYCHOPDIFF (from config file)
Analyzing file: /path/to/reduced_data.fits
Loading SKYCHOPDIFF spectra from /path/to/reduced_data.fits
  Loaded 143 SKYCHOPDIFF spectra
  Mission ID: 04_0116
Prepared spectra: 143 spectra × 1264 channels
Performing PCA decomposition (5 components)...
✓ Results saved to output/pca_components/decomposition_04_0116_20170201_components.pkl
  Total variance explained: 87.43%
  Parameters used:
    fits_file: /path/to/reduced_data.fits
    n_components: 5
    pca_source: SKYCHOPDIFF
```

## Defaults Reference

| Parameter | Default | Config Key | Command Line |
|-----------|---------|-----------|--------------|
| FITS file | `[output][reduced_fits]` | (see [input]) | `--fits` |
| `n_components` | 5 | `[pca]n_components` | `--n-components` |
| `pca_source` | SKYCHOPDIFF | `[pca]pca_source` | `--pca-source` |

## Help

To see all available options:
```bash
pca_decompose --help
```

To see detailed analysis during execution:
```bash
pca_decompose --config config.toml -v
```

Output:
```
Loading configuration from config.toml
✓ n_components = 5 (from config file)
✓ pca_source = SKYCHOPDIFF (from config file)
Loading SKYCHOPDIFF spectra from /path/to/fits/file.fits
...
```

### Override n_components Only
```bash
pca_decompose --config config.toml --n-components 10
```

Output:
```
✓ n_components = 10 (from command line)
✓ pca_source = SKYCHOPDIFF (from config file)
```

### Override pca_source Only
```bash
pca_decompose --config config.toml --pca-source M51
```

Output:
```
✓ n_components = 5 (from config file)
✓ pca_source = M51 (from command line)
```

### Override Both Parameters
```bash
pca_decompose --config config.toml --n-components 15 --pca-source M51
```

Output:
```
✓ n_components = 15 (from command line)
✓ pca_source = M51 (from command line)
```

### Use All Defaults (no config file)
```bash
pca_decompose --config config.toml
```
(If `[pca]` section is missing from config.toml)

Output:
```
✓ n_components = 5 (default)
✓ pca_source = SKYCHOPDIFF (default)
```

### Verbose Output
Add `-v` or `--verbose` flag for detailed logging:

```bash
pca_decompose --config config.toml -v
```

## Parameter Precedence

Parameters are resolved in this order (first match wins):

```
1. Command line argument  (highest priority)
   └─ pca_decompose --n-components 10

2. Config file [pca] section
   └─ [pca]
      n_components = 5

3. Built-in default      (lowest priority)
   └─ n_components = 5
```

## Configuration Examples

### Example 1: SKYCHOPDIFF Reference (Default)
```bash
# Uses: 5 components, SKYCHOPDIFF source
pca_decompose --config config.toml
```

**config.toml**:
```toml
[pca]
n_components = 5
pca_source = "SKYCHOPDIFF"
```

### Example 2: Custom M51 Reference
```bash
# Uses: 8 components, M51 source
pca_decompose --config config.toml --n-components 8 --pca-source M51
```

**config.toml** (not needed, but shown for reference):
```toml
[pca]
# These would be overridden by command line
n_components = 5
pca_source = "SKYCHOPDIFF"
```

### Example 3: Multiple Different Runs
Run PCA with different configurations sequentially:

```bash
# First: Extract 5 components from SKYCHOPDIFF
pca_decompose --config config.toml

# Then: Extract 10 components from M51
pca_decompose --config config.toml --n-components 10 --pca-source M51

# Then: Extract 8 components from different reference
pca_decompose --config config.toml --n-components 8 --pca-source OTHER_REF
```

## Finding Available Source Objects

If you're not sure what source objects are available in your FITS file:

```bash
pca_decompose --config config.toml --pca-source UNKNOWN
```

The command will show you all available objects:

```
Available objects:
  SKYCHOPDIFF: 143 spectra
  M51: 8234 spectra
  NGC253: 1256 spectra
```

## Output Files

The decomposition results are saved to:
```
output/pca_components/decomposition_<MISSION_ID>_<DATE>_components.pkl
```

Example:
```
output/pca_components/decomposition_04_0116_20170201_components.pkl
```

The output filename includes:
- `MISSION_ID`: From FITS file data
- `DATE`: Flight date from DATE-OBS header

## Example Workflow

```bash
#!/bin/bash

# Step 1: Setup config with base parameters
cat > config.toml << EOF
[input]
fits_file = "/path/to/data.fits"

[pca]
n_components = 5
pca_source = "SKYCHOPDIFF"

[output]
datacube = "/path/to/output/datacube.fits"
EOF

# Step 2: Run PCA with config defaults
pca_decompose --config config.toml

# Step 3: Try with more components
pca_decompose --config config.toml --n-components 10

# Step 4: Try with different source
pca_decompose --config config.toml --pca-source M51
```

## Logging Output

The command prints which parameters are being used and where they come from:

```
Loading configuration from config.toml
✓ n_components = 5 (from config file)
✓ pca_source = SKYCHOPDIFF (from config file)
Loading SKYCHOPDIFF spectra from /home/data/fits/file.fits
  Loaded 143 SKYCHOPDIFF spectra
  Mission ID: 04_0116
Prepared spectra: 143 spectra × 1264 channels
Performing PCA decomposition (5 components)...
✓ Results saved to output/pca_components/decomposition_04_0116_20170201_components.pkl
  Total variance explained: 87.43%
  Parameters used:
    n_components: 5
    pca_source: SKYCHOPDIFF
```

## Defaults Reference

| Parameter | Default | Config Key | Command Line |
|-----------|---------|-----------|--------------|
| `n_components` | 5 | `[pca]n_components` | `--n-components` |
| `pca_source` | SKYCHOPDIFF | `[pca]pca_source` | `--pca-source` |

## Help

To see all available options:
```bash
pca_decompose --help
```

To see detailed analysis during execution:
```bash
pca_decompose --config config.toml -v
```
