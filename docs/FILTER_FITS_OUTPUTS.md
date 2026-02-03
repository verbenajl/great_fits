# filter_fits Command - Outputs Guide

## Overview
The `filter_fits` command processes FITS spectral data and creates **two output FITS files**:

1. **clean_data.fits** (or custom `--output-clean` path)
2. **rejected_data.fits** (or custom `--output-rejected` path)

---

## Output File 1: Clean FITS File

### Purpose
Contains spectra that **passed** the filtering criteria.

### Content
- **If `--apply-only-to-object` is NOT used (default mode)**:
  - ALL non-target-object spectra (calibration data, other objects, etc.) - **UNFILTERED**
  - Target object spectra with NaN fraction **< nan_threshold**

- **If `--apply-only-to-object` IS used**:
  - ALL spectra (target and non-target) with NaN fraction **< nan_threshold**

### Data Removed From Clean File
- Target object spectra with NaN fraction **>= nan_threshold** (goes to rejected)
- Rows matching removal criteria if `--remove` and `--remove-values` are used (goes to rejected)

### Example Output
```
FITS File: clean_data.fits
├── Primary HDU (header info)
└── SPECTRA HDU
    ├── 39200 M51CENTER entries (passed NaN threshold)
    ├── 25704 SKYCHOPDIFF entries (kept, non-target object)
    ├── 1400 TSYS entries (kept, non-target object)
    └── ... other non-target objects
```

---

## Output File 2: Rejected FITS File

### Purpose
Contains spectra that **failed** the filtering criteria.

### Content - Default Mode (`--apply-only-to-object` NOT used)
- Target object spectra with NaN fraction **>= nan_threshold** (failed NaN filter)
- Rows matching removal criteria if `--remove` is specified

### Content - All-Spectra Mode (`--apply-only-to-object` IS used)
- ANY spectra (target or non-target) with NaN fraction **>= nan_threshold**
- Rows matching removal criteria if `--remove` is specified

### Example Output
```
FITS File: rejected_data.fits
├── Primary HDU (header info)
└── SPECTRA HDU
    ├── 1245 M51CENTER entries (failed NaN threshold >= 20%)
    └── Plus any rows removed by --remove criteria
```

---

## Filtering Logic Flowchart

### Default Mode: `filter_fits --config config.toml`
```
INPUT: Original 73,892 spectra (M51CENTER + calibration)
  │
  ├─→ Separate by object
  │   ├─→ M51CENTER: 39,200 entries
  │   └─→ Other objects: 34,692 entries (TSYS, TAU, SKYCHOPDIFF, etc.)
  │
  ├─→ Apply NaN filter ONLY to M51CENTER
  │   ├─→ NaN < 20%: 38,000 → CLEAN file
  │   └─→ NaN >= 20%: 1,200 → REJECTED file
  │
  └─→ OUTPUTS:
      ├─→ clean_data.fits: 38,000 M51 + 34,692 other = 72,692 total
      └─→ rejected_data.fits: 1,200 M51 (failed filter)
```

### All-Spectra Mode: `filter_fits --config config.toml --apply-only-to-object`
```
INPUT: Original 73,892 spectra
  │
  ├─→ Apply NaN filter to ALL spectra
  │   ├─→ NaN < 20%: 70,000 → CLEAN file
  │   └─→ NaN >= 20%: 3,892 → REJECTED file
  │
  └─→ OUTPUTS:
      ├─→ clean_data.fits: 70,000 total (all objects)
      └─→ rejected_data.fits: 3,892 total (all objects)
```

---

## Output Option: With Removal Criteria

### Command
```bash
filter_fits --config config.toml \
  --remove AOR_ID \
  --remove-values 04_0116_0020609 --remove-values 04_0116_0020506
```

### Results
- Rows matching the specified AOR_IDs are **removed first** (before object separation)
- These removed rows go to the **rejected file**
- NaN filtering is then applied to remaining data

### Example
```
Inputs: 73,892 total
├─→ Remove rows: 42 entries (matching AOR_IDs)
├─→ Process remaining: 73,850 entries
│   ├─→ Clean (NaN < 20%): 72,650
│   └─→ Rejected (NaN >= 20%): 1,200
│
Outputs:
├─→ clean_data.fits: 72,650
└─→ rejected_data.fits: 1,200 + 42 removed = 1,242
```

---

## Command Examples and Expected Output Files

### Example 1: Basic Filtering
```bash
filter_fits --config config.toml
```

**Output Files:**
- `clean_data.fits` - M51 with NaN < 20% + all other objects
- `rejected_data.fits` - M51 with NaN >= 20%

### Example 2: Custom Threshold
```bash
filter_fits --config config.toml --nan-threshold 0.75
```

**Output Files:**
- `clean_data.fits` - M51 with NaN < 75% + all other objects
- `rejected_data.fits` - M51 with NaN >= 75%

### Example 3: Filter Only Target Object
```bash
filter_fits --config config.toml --apply-only-to-object
```

**Output Files:**
- `clean_data.fits` - All spectra (any object) with NaN < 20%
- `rejected_data.fits` - All spectra (any object) with NaN >= 20%

### Example 4: Custom Output Paths
```bash
filter_fits --config config.toml \
  --output-clean /path/to/m51_clean.fits \
  --output-rejected /path/to/m51_rejected.fits
```

**Output Files:**
- `/path/to/m51_clean.fits`
- `/path/to/m51_rejected.fits`

### Example 5: Remove Specific AOR_IDs
```bash
filter_fits --config config.toml \
  --remove AOR_ID \
  --remove-values 04_0116_0020609 --remove-values 04_0116_0020506 \
  --output-clean m51_aor_filtered_clean.fits \
  --output-rejected m51_aor_filtered_rejected.fits
```

**Output Files:**
- `m51_aor_filtered_clean.fits` - M51 (NaN < 20%) minus specified AOR_IDs + other objects
- `m51_aor_filtered_rejected.fits` - M51 (NaN >= 20%) plus removed AOR_ID rows

---

## Output File Structure

### Both output files contain:

1. **Primary HDU**
   - Copy of original Primary HDU header
   - No data

2. **SPECTRA HDU** (Binary Table)
   - All original columns from input file:
     - OBJECT, LOFREQ, LINE, FOFFSET, RESTFREQ, VELDEF, DELTAV
     - EPOCH, DATE-OBS, SCAN, SUBSCAN, TELESCOP
     - Coordinates: LONGITUDE, LATITUDE, ALTITUDE, CDELT2, CDELT3
     - Temperatures: TSYS, TOUTSIDE, TCHOP, TCOLD
     - Environmental: ELEVATIO, AZIMUTH, PRESSURE, OBSTIME, TAU-ATM, MH2O
     - **SPECTRUM**: 1D spectral array
     - AOR_ID, MISSION_ID (if present in original)
     - Other instrumental parameters

### File Size Relationship
```
Original FITS file size = S
Clean file size ≈ S × (n_clean / n_total)
Rejected file size ≈ S × (n_rejected / n_total)
Clean size + Rejected size ≈ S
```

---

## Key Parameters Summary

| Parameter | Default | Behavior |
|-----------|---------|----------|
| `--nan-threshold` | 0.20 (20%) | Max acceptable NaN fraction in spectrum |
| `--output-clean` | `clean_data.fits` | Path to clean output FITS |
| `--output-rejected` | `rejected_data.fits` | Path to rejected output FITS |
| `--apply-only-to-object` | False | If True, filter all spectra; if False, only filter target object |
| `--remove` | None | Column name to filter by (e.g., AOR_ID) |
| `--remove-values` | None | Values to remove (can specify multiple) |

---

## Typical Workflow

```bash
# Step 1: Inspect dataset
print_oifits_info --config config.toml

# Step 2: Filter data
filter_fits --config config.toml

# Step 3: Check filtering results
ls -lh clean_data.fits rejected_data.fits

# Step 4: Use clean data for further processing
reduce_spectra --fits clean_data.fits
pca_decompose --fits clean_data.fits
map_integrated --fits clean_data.fits
```

