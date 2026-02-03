# reduce_spectra Command - "Methods Applied" Messages Guide

## Overview

The `reduce_spectra` command applies spectral reduction methods to FITS data. The "Methods applied:" section shows which reduction techniques were used in the processing.

---

## Possible Messages in "Methods applied:"

### 1. **Unblank** (Fill NaN Values)
```
Methods applied:
  - Unblank (fill NaN values with linear interpolation)
```
**When shown**: Use the `--unblank` flag
**What it does**: Fills NaN (bad) values in spectra using linear interpolation between valid channels

**Example command**:
```bash
reduce_spectra --config config.toml --unblank
```

---

### 2. **Baseline Subtraction**

#### Simple baseline (default order):
```
Methods applied:
  - Baseline subtraction (order=1)
```

#### Baseline with custom polynomial order:
```
Methods applied:
  - Baseline subtraction (order=2)
```

#### Baseline with window exclusion:
```
Methods applied:
  - Baseline subtraction (order=1, window=100-200)
```

**When shown**: Use the `--baseline` flag
**What it does**: Removes polynomial baseline from spectra
  - Default order: 1 (linear) or read from config `[reduction].baseline`
  - Window: Channels to exclude from baseline fitting (e.g., strong emission lines)

**Example commands**:
```bash
reduce_spectra --config config.toml --baseline
reduce_spectra --config config.toml --baseline --baseline-order 2
reduce_spectra --config config.toml --baseline --baseline-window 100 200
reduce_spectra --config config.toml --baseline --baseline-order 3 --baseline-window 450 550
```

---

### 3. **Smoothing**
```
Methods applied:
  - Smoothing (window=5)
```

**When shown**: Use the `--smooth` flag
**What it does**: Applies boxcar smoothing to reduce noise
  - Default window: 5 channels
  - Can be customized with `--smooth-window N`

**Example commands**:
```bash
reduce_spectra --config config.toml --smooth
reduce_spectra --config config.toml --smooth --smooth-window 7
```

---

### 4. **Extraction** (Velocity Range)
This method is **NOT shown in the "Methods applied:" section**, but is applied if configured.

**When applied**: If `[reduction].extract` is defined in config.toml
**What it does**: Extracts only a specified velocity range from the spectrum

**Example in config.toml**:
```toml
[reduction]
extract = [350, 700]  # Extract 350-700 km/s
```

**Note**: This is applied internally but not listed in the output messages. You can verify it worked by checking that the output FITS has fewer channels.

---

## Combination Examples

### **Two methods** (Unblank + Baseline):
```
✓ Spectral reduction complete.
  Output: reduced_data.fits
  Methods applied:
    - Unblank (fill NaN values with linear interpolation)
    - Baseline subtraction (order=1)
```

**Command**:
```bash
reduce_spectra --config config.toml --unblank --baseline
```

---

### **Three methods** (Unblank + Baseline + Smoothing):
```
✓ Spectral reduction complete.
  Output: reduced_data.fits
  Methods applied:
    - Unblank (fill NaN values with linear interpolation)
    - Baseline subtraction (order=2, window=450-550)
    - Smoothing (window=5)
```

**Command**:
```bash
reduce_spectra --config config.toml --unblank --baseline --baseline-order 2 --baseline-window 450 550 --smooth
```

---

### **No methods** (Blank "Methods applied:"):
```
✓ Spectral reduction complete.
  Output: reduced_data.fits
  Methods applied:
```

**Meaning**: No reduction flags were specified. Data is extracted (if configured in config.toml) but no other processing was applied.

**Command**:
```bash
reduce_spectra --config config.toml
```

**Note**: This is unusual and typically means you wanted to apply some methods but forgot the flags!

---

## Complete Method Order of Application

When multiple methods are specified, they are applied in this order:

1. **Extract** (from config, if specified) - Selects velocity range
2. **Unblank** (if `--unblank`) - Fills NaN values
3. **Baseline** (if `--baseline`) - Removes baseline
4. **Smoothing** (if `--smooth`) - Applies smoothing

---

## Configuration from config.toml

The `reduce_spectra` command reads default values from the `[reduction]` section:

```toml
[reduction]
baseline = 1              # Default baseline polynomial order
extract = [350, 700]      # Extract velocity range in km/s (optional)
baseline_window = [450, 550]  # Channels to exclude from baseline fit (optional)
```

**CLI options override config values**.

---

## Input/Output File Selection

### Input file priority:
1. `--fits /path/to/file.fits` (if specified)
2. `--clean` flag → uses `[output].clean_fits` from config
3. Default → uses `[input].fits_file` from config

### Output file:
1. `--output /path/to/output.fits` (if specified)
2. Default → uses `[output].reduced_fits` from config
3. Fallback → `reduced_data.fits` (with warning if config not available)

---

## Example Complete Workflow

### Step 1: Check dataset
```bash
print_oifits_info --config config.toml
# Output shows: M51CENTER: 39,200 entries
```

### Step 2: Filter to clean data
```bash
filter_fits --config config.toml
# Output: clean_data.fits (38,000 M51 spectra)
```

### Step 3: Reduce spectra
```bash
reduce_spectra --config config.toml --clean --unblank --baseline --baseline-order 2 --smooth
```

**Expected output**:
```
Using clean FITS file: /path/to/clean_data.fits

✓ Spectral reduction complete.
  Output: /path/to/reduced_data.fits
  Methods applied:
    - Unblank (fill NaN values with linear interpolation)
    - Baseline subtraction (order=2)
    - Smoothing (window=5)
```

### Step 4: Use reduced data for analysis
```bash
map_integrated --config config.toml --reduced --plot integrated_map.png
pca_decompose --config config.toml --reduced
```

---

## Troubleshooting

### Problem: "Methods applied:" is blank

**Cause**: No reduction flags specified

**Solution**: Add appropriate flags:
```bash
# Apply baseline subtraction
reduce_spectra --config config.toml --baseline

# Or apply multiple methods
reduce_spectra --config config.toml --unblank --baseline --smooth
```

---

### Problem: Want to apply method but it doesn't appear

**Cause**: Flag might not be recognized or config issue

**Solution**: Check:
1. Correct flag name (e.g., `--baseline` not `--subtractbaseline`)
2. Config file exists and is readable
3. For baseline: Check if config has proper baseline order defined

```bash
# Verify with explicit flags
reduce_spectra --config config.toml --baseline --baseline-order 2 -v
```

---

## Return Code and Validation

After running `reduce_spectra`, check:
1. **Exit code**: `$? = 0` means success
2. **Output file**: Verify reduced FITS was created
3. **File size**: Should be similar to input (unless extraction was applied)
4. **Data integrity**: Can open with `print_oifits_info --fits reduced_data.fits`

