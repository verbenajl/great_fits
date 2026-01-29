# filter_fits Command - Complete Analysis

## Overview

The `filter_fits` command filters FITS spectral data by object name, NaN content, and/or column value removal. It creates two output files:

1. **Clean FITS file**: All non-target objects + target object spectra with < nan_threshold NaNs
2. **Rejected FITS file**: Target object spectra with >= nan_threshold NaNs or matching removal criteria

---

## CLI Command Signature

```bash
filter_fits [OPTIONS]
```

### Options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--config` | Path | None | Path to config.toml file (required if no --fits) |
| `--fits` | Path | None | Path to FITS file to read directly (overrides config) |
| `--object` | String | Config | Object name to filter (substring match) |
| `--nan-threshold` | Float | Config | NaN threshold for filtering (default: 0.20 = 20%) |
| `--output-clean` | Path | Config | Output path for clean FITS file |
| `--output-rejected` | Path | Config | Output path for rejected FITS file |
| `--remove` | String | None | Column name to filter by (e.g., AOR_ID) |
| `--remove-values` | String | None | Values to remove from the column (repeatable) |

---

## Configuration File Usage

The `filter_fits` command reads parameters from `config.toml`:

### From `[filters]` section:
```toml
[filters]
blank_fraction = 0.20  # NaN threshold for filtering
```

### From `[parameters]` section:
```toml
[parameters]
object = "M51"  # Object name to filter
```

### From `[output]` section:
```toml
[output]
clean_fits = "/path/to/clean_data.fits"
rejected_fits = "/path/to/rejected_data.fits"
```

---

## Usage Examples

### 1. Basic Usage with Config File
```bash
filter_fits --config config.toml
```
**Behavior:**
- Reads object name from config: `M51`
- Reads NaN threshold from config: `0.20` (20%)
- Reads output paths from config
- Filters M51 spectra by NaN content
- Saves clean and rejected files

### 2. Custom NaN Threshold
```bash
filter_fits --config config.toml --nan-threshold 0.75
```
**Behavior:**
- Keeps spectra with < 75% NaN channels
- Rejects spectra with >= 75% NaN channels
- Other parameters from config

### 3. Different Object and Output Files
```bash
filter_fits --config config.toml --object "NGC1234" \
    --output-clean ngc1234_clean.fits \
    --output-rejected ngc1234_rejected.fits
```
**Behavior:**
- Filters for NGC1234 object instead of M51
- Saves to custom output paths
- Uses config NaN threshold

### 4. Remove Specific AOR_ID Values
```bash
filter_fits --config config.toml --remove AOR_ID \
    --remove-values 04_0116_0020609 --remove-values 04_0116_0020506
```
**Behavior:**
- Filters by NaN content (20% threshold)
- Additionally removes rows with AOR_ID matching specified values
- Removed rows go to rejected file

### 5. Direct FITS File Input
```bash
filter_fits --fits /path/to/data.fits --object M51 --nan-threshold 0.30
```
**Behavior:**
- Reads FITS file directly (bypasses config)
- Still uses M51 as target object
- Uses 30% NaN threshold
- Output files go to current directory (defaults)

---

## Processing Logic

### Step 1: Load Input Data
- Reads FITS file from either `--fits` or config
- Finds all HDUs with SPECTRUM column
- Concatenates data from all HDUs (handles multi-HDU files)

### Step 2: Apply Removal Filter (if specified)
- Checks `--remove` column for values in `--remove-values`
- Separates rows to remove from main data
- Removed rows → rejected file

### Step 3: Separate by Object
- Creates mask for target object (substring match in OBJECT column)
- Divides data into:
  - **target_data**: Spectra containing object name (e.g., "M51")
  - **other_data**: All other spectra

### Step 4: Filter Target Object by NaN Content
For each spectrum in target_data:
1. Count NaN channels
2. Calculate NaN fraction = NaN_count / total_channels
3. Keep if: NaN_fraction < nan_threshold
4. Reject if: NaN_fraction >= nan_threshold

### Step 5: Create Output Files

**Clean FITS file contains:**
- All non-target object spectra (100% kept)
- Target object spectra with NaN_fraction < threshold

**Rejected FITS file contains:**
- Target object spectra with NaN_fraction >= threshold
- Removed rows (if --remove filter applied)

---

## Output Files

### Clean FITS File
```
File: clean_data.fits
Contains:
  - All spectra from objects other than target
  - All target object spectra with acceptable NaN content
Used for: Scientific analysis
```

### Rejected FITS File
```
File: rejected_data.fits
Contains:
  - Target object spectra with too many NaNs
  - Rows matching removal criteria (if any)
Used for: Quality assessment and diagnostics
```

---

## Key Features

### 1. **Flexible Object Filtering**
- Substring match (case-insensitive)
- Can filter by object name from command line or config
- Example: "--object M51" will match "M51 NUCLEUS", "M51_EXTENDED", etc.

### 2. **Configurable NaN Threshold**
- Default: 20% (from config)
- Can override: `--nan-threshold 0.50` for 50%
- Applied only to target object spectra
- Non-target spectra kept regardless of NaN content

### 3. **Column-Based Removal**
- Remove rows where specified column has certain values
- Useful for removing problematic observation IDs (AOR_ID)
- Can specify multiple values to remove
- Removed rows go to rejected file for review

### 4. **Multi-HDU Support**
- Processes ALL HDUs with SPECTRUM columns
- Concatenates data from combined files
- Preserves header information

### 5. **Automatic Output Paths**
- Uses paths from config `[output]` section
- Can override with command-line options
- Creates parent directories if needed

---

## Function Implementation

### Location
`src/oi_zeigt/reduction/core.py: filter_and_save_fits()`

### Parameters
```python
def filter_and_save_fits(
    hdul: fits.HDUList,                    # Input FITS HDU list
    object_name: str,                      # Target object to filter
    nan_threshold: float = 0.20,           # NaN fraction threshold (0-1)
    output_clean: Optional[Path] = None,   # Clean output path
    output_rejected: Optional[Path] = None,# Rejected output path
    remove_column: Optional[str] = None,   # Column to filter by
    remove_values: Optional[list] = None   # Values to remove
) -> Tuple[Path, Path]
```

### Returns
```python
(clean_path, rejected_path)  # Paths to generated FITS files
```

---

## Current Configuration

From `config.toml`:
```toml
[parameters]
object = "M51"

[filters]
blank_fraction = 0.20

[output]
clean_fits = "/home/diskB/sofiaobsdata/m51/m51_central_fits/fits_files_processed/clean_data.fits"
rejected_fits = "/home/diskB/sofiaobsdata/m51/m51_central_fits/fits_files_processed/rejected_data.fits"
```

### Default Behavior
```bash
filter_fits --config config.toml
```
- Filters M51 object
- NaN threshold: 20%
- Outputs to paths specified in config

---

## Error Handling

### Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| "Object name not specified" | No --object and no config | Provide --object or set in config |
| "Column 'XXX' not found" | Invalid --remove column | Check column name exists in FITS |
| "No HDU with SPECTRUM column found" | Input file has no spectra | Use valid FITS spectral file |
| File not found | Invalid --fits path | Check file path exists |
| Permission denied | Can't write output | Check directory is writable |

---

## Performance Notes

- **Speed**: Linear in number of spectra
- **Memory**: Loads all spectra into memory (use for moderate file sizes)
- **Multi-HDU files**: Concatenates all data, increases memory usage

For very large files (>100k spectra):
- Consider filtering in multiple passes
- Or splitting input into smaller files

---

## Related Commands

- `clean_fits` - Alternative filtering interface
- `reduce_spectra` - Applies reduction after filtering
- `plot_sample_spectra` - Visualize spectra before/after filtering

---

## See Also

- `docs/velocity_axis/` - Velocity axis implementation
- `docs/features/` - Other feature documentation
- `config.toml` - Configuration reference

