# Velocity Axis Quick Reference Guide

## TL;DR

The output FITS file from `reduce_spectra()` now includes a **VELOCITY_AXIS** column that matches the SPECTRUM column, making it easy to plot and analyze spectra.

```python
# Load and plot
from astropy.io import fits
import matplotlib.pyplot as plt

with fits.open('reduced.fits') as f:
    spec = f[1].data['SPECTRUM'][0]
    vel = f[1].data['VELOCITY_AXIS'][0]
    
plt.plot(vel, spec)
plt.xlabel('Velocity (m/s)')
plt.ylabel('Flux')
plt.show()
```

---

## Velocity Axis Formula

The velocity at each spectral channel is calculated as:

$$v(i) = v_{\text{ref}} + (i - (p_{\text{ref}} - 1)) \cdot \Delta v$$

Where:
- $v(i)$ = velocity at channel $i$ [m/s]
- $v_{\text{ref}}$ = reference velocity from VELOCITY column [m/s]
- $\Delta v$ = velocity spacing from DELTAV column [m/s/channel]
- $p_{\text{ref}}$ = reference pixel from CRPIX1 [1-indexed, FITS convention]
- $i$ = channel index [0-indexed]

---

## Visual Example

### Input Parameters
```
VELOCITY = 470000 m/s
DELTAV = 500 m/s/channel
CRPIX1 = 506 (1-indexed)
Channels = 1264 total
```

### Resulting Velocity Axis

```
Channel index (0-indexed)      Velocity (m/s)
            0          →        217500
          100          →        320000
          200          →        420000
          300          →        520000
          400          →        620000
          505          →        470000  ← Reference (CRPIX1 - 1)
          506          →        470500
          600          →        520000
          700          →        620000
          800          →        720000
         1000          →        920000
         1263          →       849000
```

### Graph
```
Velocity (m/s)
    ↑
950000 |                              •
850000 |                      •
750000 |              •
650000 |      •
550000 |•
470000 |•————————— Reference (channel 505)
350000 |•
250000 |
150000 |
    0  +————————————————————————————→ Channel
       0       256      512      768    1264
```

---

## Data Structure in FITS File

### Before reduce_spectra()
```
Binary Table HDU:
  OBJECT         [n_spectra]     → Object name
  VELOCITY       [n_spectra]     → Reference velocity
  DELTAV         [n_spectra]     → Velocity spacing
  SPECTRUM       [n_spectra, n_channels] → Spectral data
```

### After reduce_spectra()
```
Binary Table HDU:
  OBJECT         [n_spectra]     → Object name (preserved)
  VELOCITY       [n_spectra]     → Reference velocity (preserved)
  DELTAV         [n_spectra]     → Velocity spacing (preserved)
  SPECTRUM       [n_spectra, n_channels] → Reduced spectra
  VELOCITY_AXIS  [n_spectra, n_channels] ← NEW! Velocity for each channel
```

---

## Common Operations

### Plot Single Spectrum with Velocity Axis
```python
from astropy.io import fits
import matplotlib.pyplot as plt

with fits.open('reduced.fits') as hdul:
    data = hdul[1].data
    idx = 0  # First spectrum
    
    spectrum = data['SPECTRUM'][idx]
    velocity = data['VELOCITY_AXIS'][idx]
    
    plt.figure(figsize=(12, 5))
    plt.plot(velocity, spectrum, 'b-', linewidth=1.5)
    plt.xlabel('Velocity (m/s)', fontsize=12)
    plt.ylabel('Flux', fontsize=12)
    plt.title(f"Spectrum {idx}: {data['OBJECT'][idx]}")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.show()
```

### Extract Velocity Range
```python
# Get data in velocity range 450000-490000 m/s
v_min, v_max = 450000, 490000

with fits.open('reduced.fits') as hdul:
    data = hdul[1].data
    spectrum = data['SPECTRUM'][0]
    velocity = data['VELOCITY_AXIS'][0]
    
    mask = (velocity >= v_min) & (velocity <= v_max)
    spectrum_slice = spectrum[mask]
    velocity_slice = velocity[mask]
    
    print(f"Extracted {np.sum(mask)} channels")
```

### Calculate Line Properties
```python
# Find velocity of maximum intensity
max_idx = np.argmax(spectrum)
velocity_of_max = velocity[max_idx]

print(f"Peak intensity at velocity: {velocity_of_max:.0f} m/s")
print(f"Peak intensity value: {spectrum[max_idx]:.3f}")
```

### Stack Spectra by Velocity
```python
# Stack multiple spectra for source averaging
with fits.open('reduced.fits') as hdul:
    data = hdul[1].data
    
    # Get all M51 spectra
    m51_mask = (data['OBJECT'] == 'M51')
    spectra = data['SPECTRUM'][m51_mask]
    velocity = data['VELOCITY_AXIS'][0]  # Same for all
    
    # Stack (mean)
    stacked = np.mean(spectra, axis=0)
    std = np.std(spectra, axis=0)
    
    plt.errorbar(velocity, stacked, yerr=std, fmt='o-')
    plt.xlabel('Velocity (m/s)')
    plt.ylabel('Mean Flux')
    plt.show()
```

---

## Properties of VELOCITY_AXIS

| Property | Value |
|----------|-------|
| **Data Type** | float64 |
| **Shape** | (n_spectra, n_channels) |
| **Units** | meters/second (m/s) |
| **Range** | ~217500 to ~849000 m/s (typical) |
| **Spacing** | Uniform: deltav per channel |
| **Reference Point** | Channel (CRPIX1 - 1) = velocity VELOCITY |
| **Convention** | LSR (Local Standard of Rest) |

---

## Checking Your Data

### Verify VELOCITY_AXIS exists
```python
with fits.open('reduced.fits') as hdul:
    cols = hdul[1].columns
    if 'VELOCITY_AXIS' in cols.names:
        print("✓ VELOCITY_AXIS found")
    else:
        print("✗ VELOCITY_AXIS not found")
```

### Print axis statistics
```python
with fits.open('reduced.fits') as hdul:
    vel = hdul[1].data['VELOCITY_AXIS'][0]
    
    print(f"Velocity range: {vel.min():.0f} to {vel.max():.0f} m/s")
    print(f"Number of channels: {len(vel)}")
    print(f"Channel spacing: {vel[1] - vel[0]:.1f} m/s")
```

---

## Troubleshooting

### "VELOCITY_AXIS not found"
- Check that VELOCITY, DELTAV, and SPECTRUM columns exist in input
- Check that reduce_spectra() ran without errors
- Try running without reduction methods to isolate issue

### "Velocity values don't look right"
- Verify CRPIX1 is correct (check header: `hdul[1].header['CRPIX1']`)
- Check VELOCITY and DELTAV values in table
- Calculate manually: v[0] = VELOCITY + (0 - (CRPIX1-1)) * DELTAV

### File size too large
- VELOCITY_AXIS adds ~10MB for 70k spectra × 1264 channels
- Consider alternative: store only header keywords, reconstruct on-the-fly
- Or: reduce precision (float32 instead of float64)

---

## Technical Details

### Implementation Location
- File: `src/oi_zeigt/reduction/core.py`
- Functions:
  - `_extract_spectral_params()` - extracts FITS parameters
  - `_create_velocity_axis()` - computes velocity array
  - `reduce_spectra()` - adds column to output

### Testing
```bash
python3 test_reduce_spectra_velocity.py
```

Expected output:
```
✓ VELOCITY_AXIS column found in output
✓ Velocity axis values are correct!
```

---

## See Also

- `VELOCITY_AXIS_IMPLEMENTATION.md` - Design details
- `VELOCITY_AXIS_IMPLEMENTATION_SUMMARY.md` - Full summary
- `src/oi_zeigt/mapping/gridding.py` - How datacube uses velocity axis
