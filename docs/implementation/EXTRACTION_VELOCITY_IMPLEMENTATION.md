# Spectral Extraction with Absolute Velocity Parameters

## Overview

Implemented velocity-based spectral extraction in the reduction pipeline. Both the extraction range and baseline window are specified as **absolute velocities in km/s** in the config file, which are automatically converted to the appropriate channel indices.

## Configuration Parameters

Add to `config.toml` in the `[reduction]` section:

```toml
[reduction]
baseline = 3                # Polynomial order for baseline subtraction
window = [450, 500]        # Baseline window: absolute velocities in km/s
extract = [350, 700]       # Extraction range: absolute velocities in km/s
```

### Parameter Meanings

- **`baseline`**: Order of polynomial for baseline subtraction
- **`window`**: Velocity range [min, max] in km/s where baseline is fit. Values are converted to channel indices using the velocity axis.
  - Example: `[450, 500]` → fits baseline using channels corresponding to 450-500 km/s
- **`extract`**: Velocity range [min, max] in km/s to extract from spectra
  - Example: `[350, 700]` → extracts channels corresponding to 350-700 km/s
  - Extraction happens BEFORE baseline subtraction

## Technical Implementation

### Velocity-to-Channel Conversion

The velocity axis is defined by three parameters (constant across all spectra):
- **VELOCITY** (velo_ref): Reference velocity (m/s) 
- **DELTAV**: Channel spacing (m/s/channel)
- **CRPIX1**: Reference pixel (1-indexed)

Formula for velocity-to-channel conversion:
```python
i = (velocity - velo_ref) / deltav + (crpix1 - 1)
```

For the M51 dataset:
- velo_ref = 470,000 m/s
- deltav = 500 m/s/channel
- crpix1 = 506 (1-indexed)
- Full velocity range: 217,600 - 849,100 m/s

### Processing Order

1. **Extraction** (if `extract` in config)
   - Velocity range [350, 700] km/s → [350,000, 700,000] m/s
   - Convert to channel indices [265, 964]
   - Extract 700 channels from each spectrum
   
2. **Baseline Subtraction** (if `baseline` in config)
   - Window [450, 500] km/s → [450,000, 500,000] m/s
   - Convert to channel indices in EXTRACTED coordinate system
   - Fit polynomial baseline using window channels
   - Subtract baseline from extracted spectra
   
3. **Create Velocity Axis**
   - Generate velocity axis for extracted channel range
   - Add VELOCITY_AXIS column to output

### Code Changes

**File: `src/oi_zeigt/reduction/core.py`**

#### `reduce_spectra_from_config()` (Lines ~1384-1418)
- Parses `window` and `extract` from config as km/s values
- Converts to m/s for internal processing
- Stores conversion information in methods dict

#### `reduce_spectra()` (Lines ~1125-1210)
- **Extraction step** (Lines ~1125-1152):
  - Gets full velocity axis from FITS
  - Extracts spectra in velocity range
  - Stores channel indices for baseline window adjustment
  
- **Baseline step** (Lines ~1154-1203):
  - Accepts window in m/s (from config)
  - Converts window to channel indices using velocity axis
  - Adjusts for extraction if applied
  - Fits and subtracts baseline
  
- **Velocity axis creation** (Lines ~1209-1210):
  - Generates velocity axis for extracted range only

#### Helper Functions
- `_extract_spectral_params()`: Extract VELOCITY, DELTAV, CRPIX1 from FITS
- `_create_velocity_axis()`: Generate velocity array
- `_velocity_to_channel_index()`: Convert velocity to channel index
- `_extract_velocity_range()`: Extract spectra in velocity range

## Usage Example

```python
from oi_zeigt.reduction.core import reduce_spectra_from_config

# Reads config.toml with extract=[350, 700] and window=[450, 500]
output_path = reduce_spectra_from_config(
    config_path='config.toml',
    output_path='reduced_data.fits',
    overwrite=True
)
```

## Output

The output FITS file contains:
- **SPECTRUM**: Extracted spectra (31,822 × 700 channels)
- **VELOCITY_AXIS**: Velocity values for each channel (m/s)
- All original columns preserved (52 total)

## Test Results

Running `test_extraction_from_config.py`:

```
Config [reduction] section:
  baseline: 3
  window: [450, 500]
  extract: [350, 700]

Extraction range from config: [350, 700] km/s
  Converted to m/s: [350000, 700000]
  Mapped to channel indices: [265, 964]
  Number of channels after extraction: 700

Baseline window: [450, 500] km/s → channels [200, 300]

Output spectra shape: (31822, 700)
VELOCITY_AXIS column present:
  Length: 700
  Velocity range: 350100.0 - 699600.0 m/s
  All spectra have identical velocity axis: True
```

## Key Points

- **Both `extract` and `window` parameters are absolute velocities in km/s**
  - `extract=[350, 700]` means extract channels corresponding to 350-700 km/s
  - `window=[450, 500]` means fit baseline using channels corresponding to 450-500 km/s
  - Window velocities remain constant even after extraction
  
- **Processing Order:**
  1. Extract spectra to velocity range [350, 700] km/s → 700 channels
  2. Fit baseline polynomial using window [450, 500] km/s 
     - Maps to channels [465, 565] in full spectrum
     - Becomes channels [200, 300] in extracted spectrum (relative to extraction start)
  3. Subtract baseline from entire extracted spectrum
  4. Generate VELOCITY_AXIS for extracted range only

- **Output Properties:**
  - 700 channels per spectrum (from 350-700 km/s)
  - VELOCITY_AXIS column reflects extracted range only (350.1-699.6 km/s)
  - All 31,822 spectra share identical velocity axis
  - Original 52 FITS columns preserved

## Demonstration

Run `demo_extraction_parameters.py` to see a detailed walkthrough of how parameters are mapped to channels.
