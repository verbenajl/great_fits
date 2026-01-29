# Adding Velocity Axis to FITS Output from reduce_spectra

## Overview

When you run `reduce_spectra()`, the output FITS file should include a VELOCITY column that provides the velocity axis for each spectrum. This allows you to correlate spectrum values with their corresponding velocities.

## How Velocity Axis is Constructed

Based on the gridding.py implementation, the velocity axis is calculated from FITS header parameters:

```python
# Parameters from FITS header/table:
velo_ref     = VELOCITY[0]      # Reference velocity (m/s)
deltav       = DELTAV[0]        # Velocity spacing per channel (m/s)
crpix1_spec  = CRPIX1           # Reference pixel (1-indexed)
nchans       = len(SPECTRUM)    # Number of spectral channels

# Velocity axis for channel i (0-indexed):
v[i] = velo_ref + (i - (crpix1_spec - 1)) * deltav
```

## Example

With these FITS parameters:
- VELOCITY = 470000 m/s (reference)
- DELTAV = 500 m/s (per channel)
- CRPIX1 = 506 (reference pixel, 1-indexed)
- n_channels = 1264

The velocity axis would be:
```
Channel 0:    v = 470000 + (0 - 505) * 500 = 217500 m/s
Channel 505:  v = 470000 + (505 - 505) * 500 = 470000 m/s (reference)
Channel 1263: v = 470000 + (1263 - 505) * 500 = 849500 m/s
```

## Implementation Plan

### 1. Create a helper function to extract spectral parameters:

```python
def _extract_spectral_params(hdul: fits.HDUList) -> dict:
    """Extract spectral axis parameters from FITS file."""
    table = hdul[1].data
    
    # Get velocity parameters
    velo_ref = float(table['VELOCITY'][0]) if 'VELOCITY' in table.names else 0.0
    deltav = float(table['DELTAV'][0]) if 'DELTAV' in table.names else 1.0
    
    # Get reference pixel
    table_header = hdul[1].header
    crpix1_spec = float(table_header.get('CRPIX1', 1.0))
    
    # Get number of channels
    nchans = table['SPECTRUM'].shape[1] if 'SPECTRUM' in table.names else 0
    
    return {
        'velo_ref': velo_ref,
        'deltav': deltav,
        'crpix1_spec': crpix1_spec,
        'nchans': nchans,
    }
```

### 2. Create velocity axis for each spectrum:

```python
def _create_velocity_axis(velo_ref, deltav, crpix1_spec, nchans):
    """Create velocity axis array."""
    channel_indices = np.arange(nchans)
    velocity_axis = velo_ref + (channel_indices - (crpix1_spec - 1)) * deltav
    return velocity_axis
```

### 3. Add VELOCITY column to output FITS:

In the `reduce_spectra()` function, after creating the reduced spectra, add a VELOCITY column:

```python
# Create velocity axes for all spectra
spectral_params = _extract_spectral_params(hdul)
velocity_axis = _create_velocity_axis(
    spectral_params['velo_ref'],
    spectral_params['deltav'],
    spectral_params['crpix1_spec'],
    spectral_params['nchans']
)

# Add VELOCITY column (same for all spectra since velocity axis is universal)
velocity_column = np.tile(velocity_axis, (len(spectra), 1))
table['VELOCITY_AXIS'] = velocity_column
```

## Data Structure

The output FITS will have:
- SPECTRUM: [n_spectra, n_channels] - Spectral values
- VELOCITY_AXIS: [n_spectra, n_channels] - Corresponding velocities

This allows for easy plotting and calculations:
```python
import matplotlib.pyplot as plt

# Load FITS data
with fits.open('output.fits') as hdul:
    data = hdul[1].data
    spectrum = data['SPECTRUM'][0]
    velocity = data['VELOCITY_AXIS'][0]
    
# Plot spectrum vs velocity
plt.plot(velocity, spectrum)
plt.xlabel('Velocity (m/s)')
plt.ylabel('Intensity')
plt.show()
```

## Implementation Notes

1. **Column Name**: Use 'VELOCITY_AXIS' to avoid confusion with scalar VELOCITY column
2. **Data Type**: Store as float32 or float64 depending on precision needs
3. **Header Keywords**: Copy VELOCITY, DELTAV, VELDEF, RESTFREQ, CRPIX1 to output header
4. **Consistency**: Same velocity axis for all spectra (unless per-spectrum velocities vary)
5. **File Size**: Each VELOCITY_AXIS column adds ~10 MB for 1264 channels × 70,000 spectra

## Alternative: Store as Header Keywords Only

If file size is a concern, store parameters in header and reconstruct axis on-the-fly:
```python
# This would allow users to reconstruct:
velocity = velo_ref + np.arange(nchans) * deltav
```

But the explicit column is cleaner for users and enables direct correlation.
