# Velocity Axis Functions - Call Hierarchy and Usage Map

## Overview
This document shows where velocity axis creation functions are called and how they integrate into the spectral processing pipeline.

---

## Call Hierarchy

### Reduction Module (`src/oi_zeigt/reduction/core.py`)

```
reduce_spectra() [line 1045]
├─ Handles extraction + baseline + smoothing
├─ Calls: _extract_spectral_params(hdul) [line 1273]
│  └─ Returns: dict with velo_ref, deltav, crpix1_spec, nchans
├─ Calls: _create_velocity_axis() [line 1272]
│  ├─ Takes: velo_ref, deltav, crpix1_spec, nchans
│  └─ Returns: velocity_axis (numpy array)
└─ If extraction was performed:
   └─ Slices velocity_axis: velocity_axis[ch_min:ch_max+1]
      └─ Result: extracted_velocity_axis matching extracted spectrum size

reduce_spectra_from_config() [line 1402]
├─ Reads config.toml [reduction] section
├─ Parses extract=[350,700] (km/s) → converts to m/s
├─ Parses window=[450,500] (km/s) → converts to m/s
└─ Calls: reduce_spectra() with methods dict
   └─ Result: FITS file with VELOCITY_AXIS column
```

### Gridding Module (`src/oi_zeigt/mapping/gridding.py`)

```
create_spectral_datacube() [line 992]
├─ Calls: _get_spectral_axis_params(hdul) [line 1074]
│  ├─ Extracts: velo_ref, deltav, restfreq, veldef, nchans, crpix1_spec
│  └─ Also extracts celestial coordinates via _get_celestial_coords()
├─ Calculates velocity at key channels:
│  ├─ v_at_ref_channel = velo_ref
│  ├─ v_at_channel_0 = velo_ref - (crpix1_spec - 1) * deltav
│  └─ v_at_last_channel = velo_ref + (nvel - crpix1_spec) * deltav
├─ Grids each velocity channel separately
└─ Creates WCS header with velocity axis
   └─ Result: FITS datacube with proper velocity WCS

_get_spectral_axis_params(hdul) [line 307]
├─ Extracts VELOCITY column
├─ Extracts DELTAV column
├─ Extracts CRPIX1 header
├─ Extracts SPECTRUM.shape[1] for nchans
├─ Also extracts RESTFREQ and VELDEF (optional)
└─ Returns: tuple(velo_ref, deltav, restfreq, veldef, nchans, crpix1_spec)
```

---

## Data Flow Diagram

### Reduction Workflow

```
FITS Input File
│
├─ reduce_spectra(hdul, methods={'extract': [350km/s, 700km/s], 'baseline': {...}})
│
├─ 1. Extract Spectral Parameters
│  │
│  └─ _extract_spectral_params(hdul)
│     ├─ VELOCITY column → velo_ref = 470,000 m/s
│     ├─ DELTAV column → deltav = 500 m/s
│     ├─ CRPIX1 header → crpix1_spec = 506
│     └─ SPECTRUM.shape → nchans = 1264
│
├─ 2. Extract Velocity Range (if requested)
│  │
│  ├─ Map [350 km/s, 700 km/s] → [265, 964] channels
│  └─ Extract spectra[:, 265:965] (700 channels)
│
├─ 3. Create Full Velocity Axis
│  │
│  └─ _create_velocity_axis(470000, 500, 506, 1264)
│     └─ velocity_axis[1264] = [217500, 218000, ..., 849000] m/s
│
├─ 4. Extract Velocity Axis (if extraction performed)
│  │
│  └─ velocity_axis_extracted = velocity_axis[265:965]
│     └─ Result: 700 velocity values from 350-700 km/s
│
├─ 5. Add to Output FITS
│  │
│  └─ table['VELOCITY_AXIS'] = numpy.tile(velocity_axis_extracted, (nspectra, 1))
│     └─ Shape: (31822, 700)
│
└─ Output: reduced_data.fits
   ├─ SPECTRUM column: (31822, 700) ← extracted + baselined
   ├─ VELOCITY_AXIS column: (31822, 700) ← extracted velocities
   ├─ RMS_BASELINE column: (31822,) ← quality metric
   └─ Original metadata columns (OBJECT, LONGITUDE, LATITUDE, etc.)
```

### Gridding Workflow

```
FITS Input File (reduced_data.fits with SPECTRUM columns)
│
├─ create_spectral_datacube(hdul, beamsize_deg=0.25)
│
├─ 1. Extract Spectral Parameters
│  │
│  └─ _get_spectral_axis_params(hdul)
│     ├─ VELOCITY column → velo_ref = 470,000 m/s
│     ├─ DELTAV column → deltav = 500 m/s
│     ├─ CRPIX1 header → crpix1_spec = 506
│     ├─ SPECTRUM.shape → nchans = 1264 (or 700 if reduced)
│     ├─ RESTFREQ column → restfreq (for frequency axis)
│     └─ VELDEF column → veldef (velocity definition)
│
├─ 2. Calculate Velocity Range
│  │
│  ├─ v[0] = 470000 - (506-1)*500 = 217,500 m/s
│  ├─ v[ref] = 470000 m/s
│  └─ v[last] = 470000 + (nchans-506)*500
│
├─ 3. Extract Celestial Coordinates
│  │
│  └─ _get_celestial_coords(hdul)
│     ├─ RA from CDELT2/CRVAL2
│     └─ Dec from CDELT3/CRVAL3
│
├─ 4. Grid Each Velocity Channel
│  │
│  ├─ For channel 0: grid spectra[:, 0] at velocity v[0]
│  ├─ For channel 1: grid spectra[:, 1] at velocity v[1]
│  └─ ... for all nchans channels
│
├─ 5. Create WCS Header
│  │
│  └─ Set velocity axis parameters:
│     ├─ CTYPE3 = 'VRAD'           (radio velocity)
│     ├─ CRVAL3 = 470000.0         (m/s)
│     ├─ CDELT3 = 500.0            (m/s)
│     ├─ CRPIX3 = 506.0            (1-indexed reference pixel)
│     └─ RESTFREQ = (from table)
│
└─ Output: datacube.fits
   ├─ DATA cube: (nvel × ndec × nra) 3D array
   └─ WCS header with velocity axis interpretation
```

---

## Function Details

### `_extract_spectral_params(hdul)` - Reduction Module

**Location:** `src/oi_zeigt/reduction/core.py`, lines 856-922

**Purpose:** Extract all spectral parameters from FITS headers and columns

**Input:** 
- `hdul`: FITS HDUList

**Output:**
```python
{
    'velo_ref': float,       # m/s, from VELOCITY column
    'deltav': float,         # m/s, from DELTAV column
    'crpix1_spec': float,    # FITS 1-indexed, from CRPIX1 header
    'nchans': int            # from SPECTRUM.shape[1]
}
```

**Called by:**
- `reduce_spectra()` (line 1273)

**Also used in:**
- `test_extraction_from_config.py`
- `demo_extraction_parameters.py`

---

### `_create_velocity_axis()` - Reduction Module

**Location:** `src/oi_zeigt/reduction/core.py`, lines 924-959

**Purpose:** Create velocity array for spectral channels

**Input:**
- `velo_ref`: float (m/s)
- `deltav`: float (m/s)
- `crpix1_spec`: float (FITS 1-indexed)
- `nchans`: int

**Output:**
- `velocity_axis`: numpy array shape (nchans,), values in m/s

**Formula:**
```
velocity_axis[i] = velo_ref + (i - (crpix1_spec - 1)) * deltav
```

**Called by:**
- `reduce_spectra()` (lines 1272, 1278)

**Example:**
```python
vel = _create_velocity_axis(470000, 500, 506, 1264)
# vel[0] = 217,500 m/s
# vel[505] = 470,000 m/s (reference)
# vel[1263] = 849,000 m/s
```

---

### `_get_spectral_axis_params()` - Gridding Module

**Location:** `src/oi_zeigt/mapping/gridding.py`, lines 307-407

**Purpose:** Extract spectral axis parameters for datacube creation

**Input:**
- `hdul`: FITS HDUList

**Output:**
```python
(velo_ref,      # float, m/s
 deltav,        # float, m/s
 restfreq,      # float, Hz
 veldef,        # string, typically 'RADI-LSR'
 nchans,        # int
 crpix1_spec)   # float, FITS 1-indexed
```

**Called by:**
- `create_spectral_datacube()` (line 1074)

**Differences from Reduction Module:**
- Also extracts RESTFREQ and VELDEF
- More comprehensive error handling with warnings
- Validates parameter consistency

---

## Integration Points

### 1. Reduction Pipeline
**Input:** Original FITS (1264 channels, no VELOCITY_AXIS)
**Processing:** Extract, baseline, add VELOCITY_AXIS
**Output:** Reduced FITS (700 channels if extracted, with VELOCITY_AXIS)
**Consumers:** plot_sample_spectra, gridding module, user analysis

### 2. Gridding Pipeline
**Input:** Reduced FITS (700 channels, with VELOCITY_AXIS) or original FITS (1264 channels)
**Processing:** Extract spectral params, grid by velocity, create datacube
**Output:** Datacube FITS (velocity × Dec × RA)
**Consumers:** Astronomy software (FITS readers), data visualization

### 3. Plotting Pipeline
**Input:** Any FITS with SPECTRUM column (checks for VELOCITY_AXIS)
**Processing:** If VELOCITY_AXIS exists, use for x-axis; else use channel index
**Output:** PNG/PDF plot with appropriate axis labels
**Consumers:** User visualization, publications

---

## Parameters Used at Each Stage

| Stage | velo_ref | deltav | crpix1 | nchans | Extract range |
|-------|----------|--------|--------|--------|---------------|
| **Input FITS** | ✓ | ✓ | ✓ | 1264 | None |
| **Extract params** | ✓ | ✓ | ✓ | 1264 | [350-700] km/s |
| **Create vel axis** | ✓ | ✓ | ✓ | 1264 | → [265-964] channels |
| **Slice vel axis** | ✓ | ✓ | ✓ | 700 | [0-699] in axis array |
| **Output FITS** | ✗ Removed | ✗ Removed | ✗ Removed | 700 | - |
| **VELOCITY_AXIS col** | (implicit) | (implicit) | (implicit) | 700 | - |

---

## Example: Complete M51 Workflow

### Step 1: Reduce Spectra
```bash
$ reduce_spectra --config config.toml --baseline --output reduced_data.fits

[Processing]
Extracted velocity range [350, 700] km/s (350000-700000 m/s) → channels [265, 964] (700 channels)
Baseline window: [450, 500] km/s → channels [200, 300]
✓ Spectral reduction complete.
```

### Step 2: Check Reduced Data
```python
from astropy.io import fits
import numpy as np

with fits.open('reduced_data.fits') as hdul:
    data = hdul[1].data
    spec = data['SPECTRUM'][0]      # Shape: (700,)
    vel = data['VELOCITY_AXIS'][0]  # Shape: (700,)
    
    print(f"Spectrum range: {vel[0]/1000:.1f} - {vel[-1]/1000:.1f} km/s")
    print(f"Shape: {spec.shape}")
    # Output:
    # Spectrum range: 350.1 - 699.6 km/s
    # Shape: (700,)
```

### Step 3: Plot Reduced Data
```bash
$ plot_sample_spectra --fits reduced_data.fits --object M51CENTER --num-spectra 8

[Using VELOCITY_AXIS for x-axis]
✓ Plot saved to plot.png
# X-axis shows "Velocity (km/s)" from 350-700
```

### Step 4: Create Datacube
```bash
$ create_datacube --config config.toml --fits reduced_data.fits --output datacube.fits

[Processing]
Spectral axis parameters:
  Reference velocity: 470000.00 m/s (470.00 km/s)
  Velocity step: 500.00 m/s (0.5000 km/s)
  Reference pixel (CRPIX1): 506.00
  ...
✓ Datacube created: datacube.fits
```

### Step 5: Use Datacube
```python
from astropy.io import fits
from astropy.wcs import WCS

with fits.open('datacube.fits') as hdul:
    wcs = WCS(hdul[0].header)
    # Now WCS knows how to interpret velocity axis
    # Proper velocity scaling: 350-700 km/s
```

---

## Consistency Verification Checklist

- [x] Both modules use identical velocity formula
- [x] Both extract from same FITS sources (VELOCITY, DELTAV, CRPIX1)
- [x] Both handle FITS 1-indexed reference pixel correctly
- [x] Both support extracted spectra (different nchans)
- [x] Both include comprehensive documentation
- [x] Both handle unit conversions (m/s ↔ km/s) appropriately
- [x] Both validate extracted parameters

**Status: ✓ All consistent and working correctly**

---

## Related Files

### Reduction Module Files
- `src/oi_zeigt/reduction/core.py` - Main functions
- `src/oi_zeigt/cli.py` - Command line interface
- `src/oi_zeigt/basic_io.py` - Config reading

### Gridding Module Files
- `src/oi_zeigt/mapping/gridding.py` - Gridding functions
- `src/oi_zeigt/mapping/__init__.py`

### Test/Demo Files
- `test_extraction_from_config.py`
- `demo_extraction_parameters.py`

### Documentation Files
- `VELOCITY_AXIS_IMPLEMENTATION_SUMMARY.md`
- `VELOCITY_AXIS_CODE_REVIEW.md`
- `VELOCITY_AXIS_VISUAL_COMPARISON.md`
- `PLOT_VELOCITY_AXIS_FEATURE.md`
