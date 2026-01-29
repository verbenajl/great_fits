# filter_fits Code Walkthrough

## File Structure

```
src/oi_zeigt/
├── cli.py                          # CLI command definition
└── reduction/
    └── core.py                     # Implementation function
```

---

## CLI Command Definition

### Location: `src/oi_zeigt/cli.py` (Lines 663-880)

```python
@click.command()
@click.option("--config", type=click.Path(exists=True), ...)
@click.option("--fits", type=click.Path(exists=True), ...)
@click.option("--object", default=None, ...)
@click.option("--nan-threshold", type=float, default=None, ...)
@click.option("--output-clean", type=click.Path(), ...)
@click.option("--output-rejected", type=click.Path(), ...)
@click.option("--remove", type=str, default=None, ...)
@click.option("--remove-values", type=str, multiple=True, ...)
def filter_fits(config, fits, object, nan_threshold, 
                output_clean, output_rejected, remove, remove_values):
    """Filter FITS data by object, NaN content, and/or column values."""
```

### Step 1: Load Configuration

```python
# Always try to load config
try:
    if config:
        config_data = get_config(config)
    else:
        config_data = get_config()
except FileNotFoundError:
    config_data = {}
```

**What happens:**
- Tries to load from `--config` argument, or default location
- Continues gracefully if config file not found
- Extracts parameters for later use

### Step 2: Read FITS File

```python
if fits:
    click.echo(f"Reading FITS file: {fits}")
    hdul = read_fits(fits)
elif config:
    click.echo(f"Reading config from: {config}")
    hdul = read_fits_from_config(config)
else:
    click.echo("Reading config from default location...")
    hdul = read_fits_from_config()
```

**Priority:**
1. Direct FITS file via `--fits` (highest priority)
2. FITS file from config via `--config`
3. FITS file from default config location (lowest priority)

### Step 3: Get Object Name

```python
# Get object name from argument or config
if object is None and config:
    try:
        object = config_data.get("parameters", {}).get("object")
    except (NameError, KeyError):
        pass

# Repeated logic for other case
if object is None and not fits:
    try:
        object = config_data.get("parameters", {}).get("object")
    except (NameError, KeyError):
        pass

if object is None:
    raise ValueError("Object name not specified...")
```

**Priority:**
1. `--object` command-line argument (if provided)
2. `object` from config `[parameters]` section
3. Error if neither available

### Step 4: Get NaN Threshold

```python
if nan_threshold is None:
    if config:
        try:
            nan_threshold = config_data.get("filters", {}).get("blank_fraction", 0.20)
        except (NameError, KeyError):
            nan_threshold = 0.20
    else:
        nan_threshold = 0.20
```

**Priority:**
1. `--nan-threshold` command-line argument
2. `blank_fraction` from config `[filters]` section
3. Default: 0.20 (20%)

### Step 5: Get Output Paths

```python
# Get output paths from config if not specified as arguments
if output_clean is None and config:
    try:
        output_clean = config_data.get("output", {}).get("clean_fits")
    except (NameError, KeyError):
        pass

if output_rejected is None and config:
    try:
        output_rejected = config_data.get("output", {}).get("rejected_fits")
    except (NameError, KeyError):
        pass
```

**Priority:**
1. `--output-clean` and `--output-rejected` arguments
2. `clean_fits` and `rejected_fits` from config `[output]` section
3. Defaults to None (handled by function)

### Step 6: Call Implementation Function

```python
clean_path, rejected_path = filter_and_save_fits(
    hdul,
    object_name=object,
    nan_threshold=nan_threshold,
    output_clean=output_clean,
    output_rejected=output_rejected,
    remove_column=remove,
    remove_values=list(remove_values) if remove_values else None
)
```

**Passes:**
- FITS HDU list
- Target object name
- NaN threshold
- Output paths
- Removal criteria (column + values)

### Step 7: Print Statistics

```python
hdul_clean = fits_lib.open(clean_path)
hdul_rejected = fits_lib.open(rejected_path)

n_clean = len(hdul_clean[1].data) if len(hdul_clean) > 1 else 0
n_rejected = len(hdul_rejected[1].data) if len(hdul_rejected) > 1 and hdul_rejected[1].data is not None else 0

# Print formatted results
click.echo("="*70)
click.echo(click.style("✓ Successfully created FITS files", fg="green"))
click.echo(f"Clean FITS file: {clean_path}")
click.echo(f"  Records: {n_clean}")
click.echo(f"Rejected FITS file: {rejected_path}")
click.echo(f"  Records: {n_rejected}")
```

---

## Core Implementation Function

### Location: `src/oi_zeigt/reduction/core.py` (Lines 360-550)

```python
def filter_and_save_fits(hdul, object_name, nan_threshold=0.20,
                         output_clean=None, output_rejected=None,
                         remove_column=None, remove_values=None):
```

### Step 1: Find All HDUs with SPECTRUM

```python
# Find ALL binary table HDUs with spectra
matrix_hdus = []
for idx, hdu in enumerate(hdul):
    if hasattr(hdu, 'data') and hdu.data is not None:
        if 'SPECTRUM' in hdu.data.dtype.names:
            matrix_hdus.append((idx, hdu))

if not matrix_hdus:
    raise ValueError("No HDU with SPECTRUM column found")
```

**What it does:**
- Iterates through all HDUs
- Identifies HDUs with SPECTRUM column
- Stores index and HDU reference
- Raises error if none found

**Why it matters:**
- Handles multi-HDU files (e.g., combined observations)
- Extracts all spectral data for processing

### Step 2: Combine Data from All HDUs

```python
# Combine data from all HDUs
all_data = []
for idx, hdu in matrix_hdus:
    all_data.append(hdu.data)

# Concatenate all data
data = np.concatenate(all_data)
```

**What it does:**
- Extracts data from each HDU
- Concatenates into single array
- Creates unified dataset for processing

**Example:**
```
HDU1: 1000 spectra
HDU2: 2000 spectra
HDU3: 500 spectra
Result: 3500 spectra in single array
```

### Step 3: Apply Removal Filter (if specified)

```python
removed_data = np.array([])  # Track removed rows for rejected file

if remove_column and remove_values:
    if remove_column not in data.dtype.names:
        raise ValueError(f"Column '{remove_column}' not found in FITS data")
    
    # Create a mask for values to KEEP (inverse of remove)
    keep_mask = np.ones(len(data), dtype=bool)
    for i, row in enumerate(data):
        value = row[remove_column]
        # Handle both bytes and string values
        if isinstance(value, bytes):
            value = value.decode('utf-8').strip()
        else:
            value = str(value).strip()
        
        if value in remove_values:
            keep_mask[i] = False
    
    # Save the removed data for the rejected file
    removed_data = data[~keep_mask]
    data = data[keep_mask]
```

**What it does:**
1. Validates that removal column exists
2. Creates boolean mask for rows to keep
3. Iterates through each row
4. Checks if column value matches any removal value
5. Separates removed rows for rejected file
6. Keeps only non-removed rows in main data

**Example:**
```
Input: 1000 rows
Remove AOR_ID = '04_0116_0020609'
Found: 50 rows matching
Result:
  - removed_data: 50 rows
  - data: 950 rows
```

**Important:**
- Handles both bytes and string values
- Case-sensitive value matching
- Removed rows tracked separately

### Step 4: Separate by Object

```python
# Separate data by object
target_mask = np.array([object_name.lower() in str(obj).lower() 
                       for obj in data['OBJECT']])
other_mask = ~target_mask

other_data = data[other_mask]
target_data = data[target_mask]
```

**What it does:**
1. Creates boolean mask for target object (substring match)
2. Case-insensitive comparison
3. Splits data into two groups:
   - **target_data**: Contains target object name
   - **other_data**: Does not contain target object name

**Example:**
```
object_name = "M51"
OBJECT column values: ["M51 NUCLEUS", "NGC1234", "M51_EXTENDED", ...]

target_mask: [True, False, True, ...]
target_data: All rows with "M51" in OBJECT
other_data: All rows without "M51" in OBJECT
```

### Step 5: Filter Target Object by NaN Content

```python
# Filter target object by NaN content
nan_fractions = []
for spectrum in target_data['SPECTRUM']:
    _, frac = detect_nan_channels(spectrum)
    nan_fractions.append(frac)

nan_fractions = np.array(nan_fractions)
clean_target_mask = nan_fractions < nan_threshold
rejected_target_mask = ~clean_target_mask

target_clean = target_data[clean_target_mask]
target_rejected = target_data[rejected_target_mask]
```

**What it does:**
1. Iterates through each target object spectrum
2. Calls `detect_nan_channels()` to count NaNs
3. Calculates NaN fraction for each spectrum
4. Creates two masks:
   - **clean_target_mask**: NaN_fraction < threshold
   - **rejected_target_mask**: NaN_fraction >= threshold
5. Separates into two arrays

**Function: detect_nan_channels()**
```python
def detect_nan_channels(spectrum):
    """
    Parameters:
        spectrum: 1D array of spectral values
    
    Returns:
        nan_mask: Boolean array marking NaN locations
        nan_frac: Fraction of NaN channels (0-1)
    """
    nan_mask = np.isnan(spectrum)
    nan_frac = np.sum(nan_mask) / len(spectrum)
    return nan_mask, nan_frac
```

**Example:**
```
target_data: 100 M51 spectra
nan_threshold = 0.20 (20%)

For each spectrum, calculate NaN fraction:
  Spectrum 0: 5/1264 channels = 0.0039 (keep)
  Spectrum 1: 251/1264 channels = 0.199 (keep)
  Spectrum 2: 253/1264 channels = 0.200 (REJECT)
  Spectrum 3: 400/1264 channels = 0.316 (REJECT)
  ...

Result:
  target_clean: 70 spectra (NaN < 20%)
  target_rejected: 30 spectra (NaN >= 20%)
```

### Step 6: Combine Rejected Data

```python
# Add the removed data to rejected
if len(removed_data) > 0:
    target_rejected = np.concatenate([target_rejected, removed_data])
```

**What it does:**
- Adds removed rows (from Step 3) to rejected spectra
- Creates comprehensive rejected file

**Example:**
```
target_rejected: 30 spectra (NaN >= 20%)
removed_data: 50 spectra (matched removal criteria)

Result: 80 total rejected spectra
```

### Step 7: Combine Clean Data

```python
# Combine other objects with clean target objects
clean_combined = np.concatenate([other_data, target_clean])
```

**What it does:**
- Combines non-target objects with clean target spectra
- Creates comprehensive clean file

**Logic:**
```
other_data: All non-M51 spectra (kept as-is)
+ 
target_clean: M51 spectra with < 20% NaNs

= 
clean_combined: All spectra kept for analysis
```

### Step 8: Set Default Output Paths

```python
# Set default output paths
if output_clean is None:
    output_clean = Path("clean_data.fits")
else:
    output_clean = Path(output_clean)

if output_rejected is None:
    output_rejected = Path("rejected_data.fits")
else:
    output_rejected = Path(output_rejected)
```

**What it does:**
- Uses provided paths or defaults
- Converts to Path objects for consistency

### Step 9: Create FITS Files

```python
# Copy primary HDU
primary_hdu = hdul[0].copy()

# Create binary table HDU for clean data
from astropy.table import Table
table_clean = Table(clean_combined)
hdu_clean = fits.BinTableHDU(table_clean)
hdu_clean.name = matrix_hdus[0][1].name

# Copy header information from original
for key in matrix_hdus[0][1].header:
    if key not in ['NAXIS1', 'NAXIS2', 'TFIELDS'] and key != '':
        try:
            hdu_clean.header[key] = matrix_hdus[0][1].header[key]
        except (ValueError, KeyError):
            pass

# Create FITS file for clean data
hdul_clean = fits.HDUList([primary_hdu, hdu_clean])
hdul_clean.writeto(output_clean, overwrite=True)
```

**What it does:**
1. Copies primary HDU from original
2. Converts numpy array to astropy Table
3. Creates binary table HDU
4. Preserves original HDU name
5. Copies header information (except dimension-related keys)
6. Creates HDUList with primary + table
7. Writes to file with overwrite

**Why preserve header:**
- Maintains observational metadata
- Preserves WCS information
- Keeps instrument configuration

### Step 10: Create Rejected File

```python
# Create binary table HDU for rejected data (only if there are rejected spectra)
if len(target_rejected) > 0:
    table_rejected = Table(target_rejected)
    hdu_rejected = fits.BinTableHDU(table_rejected)
    hdu_rejected.name = matrix_hdus[0][1].name
    
    # Copy header information from original
    for key in matrix_hdus[0][1].header:
        if key not in ['NAXIS1', 'NAXIS2', 'TFIELDS'] and key != '':
            try:
                hdu_rejected.header[key] = matrix_hdus[0][1].header[key]
            except (ValueError, KeyError):
                pass
    
    # Create FITS file for rejected data
    hdul_rejected = fits.HDUList([primary_hdu, hdu_rejected])
    hdul_rejected.writeto(output_rejected, overwrite=True)
else:
    # Create empty rejected file if no rejections
    primary_rejected = fits.PrimaryHDU()
    hdul_rejected = fits.HDUList([primary_rejected])
    hdul_rejected.writeto(output_rejected, overwrite=True)
```

**What it does:**
- Same as clean file creation
- Creates empty file if no rejections (for consistency)
- Ensures rejected file always exists

### Step 11: Return Paths

```python
return output_clean, output_rejected
```

---

## Data Flow Diagram

```
FITS Input File
      ↓
[Load Config]
      ↓
[Read FITS → Concatenate all HDUs]
      ↓
[Apply Removal Filter (if --remove specified)]
      ├─→ removed_data → Rejected File
      │
[Separate by Object]
      ├─→ other_data (non-target objects)
      │
[Filter Target by NaN]
      ├─→ target_clean (NaN < threshold)   ┐
      ├─→ target_rejected (NaN >= threshold)─→ Rejected File
      │
[Combine for Output]
      ├─→ other_data + target_clean → Clean File
      └─→ target_rejected + removed_data → Rejected File
```

---

## Key Variables

| Variable | Type | Purpose |
|----------|------|---------|
| `hdul` | fits.HDUList | Input FITS file |
| `matrix_hdus` | list | Found HDUs with SPECTRUM |
| `data` | numpy array | Combined spectral data |
| `target_mask` | boolean array | Identifies target object rows |
| `other_mask` | boolean array | Identifies non-target rows |
| `nan_fractions` | float array | NaN fraction for each spectrum |
| `clean_target_mask` | boolean array | Identifies spectra to keep |
| `removed_data` | numpy array | Rows matching removal criteria |
| `target_clean` | numpy array | Target spectra with < NaN threshold |
| `target_rejected` | numpy array | Target spectra with >= NaN threshold |
| `clean_combined` | numpy array | All spectra to keep |
| `output_clean` | Path | Output file path (clean) |
| `output_rejected` | Path | Output file path (rejected) |

---

## Error Handling

### Validation Points

```python
# Point 1: HDU validation
if not matrix_hdus:
    raise ValueError("No HDU with SPECTRUM column found")

# Point 2: Column validation
if remove_column not in data.dtype.names:
    raise ValueError(f"Column '{remove_column}' not found in FITS data")

# Point 3: Object validation
if object is None:
    raise ValueError("Object name not specified...")
```

### Try-Except Blocks

```python
# Config loading
try:
    config_data = get_config(config)
except FileNotFoundError:
    config_data = {}  # Continue gracefully

# Header copying
try:
    hdu_clean.header[key] = original_header[key]
except (ValueError, KeyError):
    pass  # Skip problematic keys
```

---

## Performance Considerations

### Time Complexity
- Loading HDUs: O(n) where n = number of HDUs
- Concatenating data: O(total_spectra)
- Filtering: O(total_spectra)
- Creating files: O(total_spectra)

**Overall: O(total_spectra) - linear time**

### Memory Usage
- Loads all spectra into memory
- Copies data for clean and rejected files
- **Total: ~3x input file size** (input + clean + rejected)

### Optimization Tips
- For large files (>100k spectra): Consider streaming approach
- For many small files: Process in batches
- Parallel processing: Possible with data chunking

