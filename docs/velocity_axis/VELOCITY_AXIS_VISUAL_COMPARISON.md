# Velocity Axis Creation - Visual Comparison

## Quick Reference

### The Formula (Both Modules)
```
velocity[channel_i] = VELOCITY_ref + (i - (CRPIX1 - 1)) * DELTAV

Where:
  VELOCITY_ref = reference velocity from FITS (m/s)
  CRPIX1       = reference pixel from FITS header (1-indexed)
  DELTAV       = velocity step per channel from FITS (m/s)
  i            = channel index (0-indexed Python convention)
```

---

## Side-by-Side Function Comparison

### Reduction Module: `_create_velocity_axis()`

```python
# File: src/oi_zeigt/reduction/core.py, lines 924-959
def _create_velocity_axis(velo_ref: float, deltav: float, crpix1_spec: float,
                         nchans: int) -> np.ndarray:
    """Create velocity axis array for spectral data."""
    channel_indices = np.arange(nchans, dtype=np.float64)
    velocity_axis = velo_ref + (channel_indices - (crpix1_spec - 1.0)) * deltav
    return velocity_axis
```

**Input:**
- `velo_ref`: Reference velocity (m/s)
- `deltav`: Velocity step (m/s)
- `crpix1_spec`: Reference pixel (1-indexed)
- `nchans`: Total number of channels

**Output:**
- Numpy array of velocities (m/s), shape: (nchans,)

**Used by:**
- `reduce_spectra()` to create VELOCITY_AXIS column

---

### Gridding Module: Velocity Calculation

```python
# File: src/oi_zeigt/mapping/gridding.py, lines 1074-1098
velo_ref, deltav, restfreq, veldef, nchans, crpix1_spec = _get_spectral_axis_params(hdul)

# Explicit velocity calculations:
v_at_ref_channel = velo_ref
v_at_channel_0 = velo_ref - (crpix1_spec - 1) * deltav
v_at_last_channel = velo_ref + (nvel - crpix1_spec) * deltav

# For full datacube creation (implicit velocity axis):
# velocity[i] = velo_ref + (i - (crpix1_spec - 1)) * deltav
```

**Input (from `_get_spectral_axis_params()`):**
- Extracted from FITS file directly
- VELOCITY, DELTAV, CRPIX1, SPECTRUM columns

**Output:**
- Individual velocity values at key channels
- Used for WCS header creation

**Used by:**
- `create_spectral_datacube()` to set up datacube WCS

---

## Data Flow Comparison

### Reduction Module Flow
```
FITS Input File
    ↓
_extract_spectral_params()
    ↓ returns: dict with velo_ref, deltav, crpix1_spec, nchans
    ↓
_create_velocity_axis()
    ↓ returns: numpy array of velocities
    ↓
VELOCITY_AXIS column (in reduced FITS output)
    ↓
Used by: plot_sample_spectra() for x-axis labeling
```

### Gridding Module Flow
```
FITS Input File
    ↓
_get_spectral_axis_params()
    ↓ returns: velo_ref, deltav, restfreq, veldef, nchans, crpix1_spec
    ↓
create_spectral_datacube()
    ↓ calculates v at key channels
    ↓
WCS Header (with CRVAL3, CDELT3, CTYPE3)
    ↓
Datacube FITS file
    ↓
Used by: Astronomy software (FITS readers) for axis interpretation
```

---

## Parameter Extraction Comparison

### Where Parameters Come From

| Module | Parameter | Source | Location |
|--------|-----------|--------|----------|
| **Reduction** | `velo_ref` | VELOCITY column | table['VELOCITY'][0] |
| **Reduction** | `deltav` | DELTAV column | table['DELTAV'][0] |
| **Reduction** | `crpix1_spec` | CRPIX1 header | hdul[1].header['CRPIX1'] |
| **Reduction** | `nchans` | SPECTRUM shape | table['SPECTRUM'].shape[1] |
| | | | |
| **Gridding** | `velo_ref` | VELOCITY column | table['VELOCITY'][0] |
| **Gridding** | `deltav` | DELTAV column | table['DELTAV'][0] |
| **Gridding** | `crpix1_spec` | CRPIX1 header | header['CRPIX1'] |
| **Gridding** | `nchans` | SPECTRUM shape | table['SPECTRUM'].shape[1] |
| **Gridding** | `restfreq` | RESTFREQ column | table['RESTFREQ'][0] |
| **Gridding** | `veldef` | VELDEF column | table['VELDEF'][0] |

✓ **All parameters come from the same FITS sources**

---

## M51 Dataset Verification

### Input Parameters (from FITS)
```
VELOCITY  = 470,000 m/s  (reference velocity)
DELTAV    = 500 m/s      (velocity step per channel)
CRPIX1    = 506.0        (reference pixel, 1-indexed FITS)
SPECTRUM  = (31822, 1264) shape (1264 channels per spectrum)
```

### Reduction Module Results
```python
velocity_axis = _create_velocity_axis(470000, 500, 506, 1264)

# Result:
velocity_axis[0]    = 470000 + (0 - 505) * 500      = 217,500 m/s
velocity_axis[505]  = 470000 + (505 - 505) * 500    = 470,000 m/s  ← Reference
velocity_axis[1263] = 470000 + (1263 - 505) * 500   = 849,000 m/s
```

### Gridding Module Results
```python
v_at_channel_0 = 470000 - (506 - 1) * 500 = 217,500 m/s      ✓ Matches
v_at_ref = 470000 m/s                                        ✓ Matches
v_at_last = 470000 + (1264 - 506) * 500 = 849,000 m/s       ✓ Matches
```

### After Extraction [350, 700] km/s

**Extraction range:**
- 350 km/s = 350,000 m/s → channel 265
- 700 km/s = 700,000 m/s → channel 964
- Extracted channels: 265 to 964 (700 channels)

**Reduction module extracts velocity axis:**
```python
velocity_axis_full = _create_velocity_axis(470000, 500, 506, 1264)
velocity_axis_extracted = velocity_axis_full[265:965]
# Result: velocity_axis_extracted[0] = 350,000 m/s (350 km/s)
#         velocity_axis_extracted[-1] = 699,500 m/s (~700 km/s)
```

---

## Consistency Checklist

### ✓ Mathematical Formula
- [x] Reduction: `velo_ref + (i - (crpix1 - 1)) * deltav`
- [x] Gridding: `velo_ref - (crpix1 - 1) * deltav` for v_at_0, etc.
- [x] Formula components match exactly

### ✓ FITS Convention Handling
- [x] Both use CRPIX1 (1-indexed) from FITS header
- [x] Both convert to 0-indexed properly with `(crpix1_spec - 1)`
- [x] Reference pixel interpretation is consistent

### ✓ Data Source Consistency
- [x] Both extract VELOCITY from column 0 (assumed constant)
- [x] Both extract DELTAV from column 0 (assumed constant)
- [x] Both use CRPIX1 from FITS header
- [x] Both count channels from SPECTRUM column

### ✓ Handling of Extracted Spectra
- [x] Reduction: Properly slices velocity axis for extracted range
- [x] Gridding: Can handle any SPECTRUM.shape[1] nchans
- [x] Both would work correctly with reduced/extracted data

### ✓ Documentation
- [x] Both document formula clearly
- [x] Both explain FITS 1-indexed convention
- [x] Both include examples
- [x] Both have comprehensive docstrings

### ✓ Error Handling
- [x] Reduction: Implicit (assumes valid FITS)
- [x] Gridding: Explicit with warnings for missing data
- [x] Both have fallback defaults

---

## Usage Examples

### Example 1: Plotting Reduced Data with Velocity Axis

```python
# Reduction module creates VELOCITY_AXIS column
from astropy.io import fits

with fits.open('reduced_data.fits') as hdul:
    data = hdul[1].data
    spectrum = data['SPECTRUM'][0]         # Shape: (700,)
    vel_axis = data['VELOCITY_AXIS'][0]    # Shape: (700,)
    
    # Plot with velocity axis
    import matplotlib.pyplot as plt
    plt.plot(vel_axis / 1000, spectrum)    # Convert m/s to km/s
    plt.xlabel('Velocity (km/s)')
    plt.ylabel('Intensity')
```

**Result:** X-axis shows 350-700 km/s for extracted spectrum

---

### Example 2: Creating Datacube WCS

```python
# Gridding module uses velocity parameters for WCS
from oi_zeigt.mapping.gridding import create_spectral_datacube

datacube, header, fig = create_spectral_datacube(
    hdul,
    beamsize_deg=0.25,
    output_file='datacube.fits'
)

# Header now contains proper velocity WCS:
# CTYPE3   = 'VRAD    '           / Velocity (radio convention)
# CRVAL3   = 470000.0             / Reference velocity (m/s)
# CDELT3   = 500.0                / Velocity step (m/s)
# CRPIX3   = 506.0                / Reference pixel (1-indexed)
```

**Result:** Datacube has proper velocity axis for WCS-aware astronomy software

---

## Summary Table

| Aspect | Reduction | Gridding | Status |
|--------|-----------|----------|--------|
| **Formula** | `v[i] = velo_ref + (i - (crpix1-1)) * deltav` | Same | ✓ Identical |
| **Parameters** | From FITS columns/header | From FITS columns/header | ✓ Consistent |
| **Output** | Velocity array (m/s) | WCS header values | ✓ Compatible |
| **Extracted spectra** | Explicitly supported | Implicit support | ✓ Both work |
| **Documentation** | ✓ Complete | ✓ Complete | ✓ Good |
| **Error handling** | Implicit | Explicit | ✓ Adequate |

---

## Recommendation

**No changes needed.** Both modules correctly implement velocity axis creation with:
1. ✓ Identical mathematical formula
2. ✓ Consistent parameter extraction from FITS
3. ✓ Proper FITS convention handling
4. ✓ Support for extracted spectra
5. ✓ Clear documentation

The implementation is correct, consistent, and production-ready.
