# Velocity Axis Reconstruction for Spectral Datacube

## Summary of Changes

I have successfully modified the `create_spectral_datacube()` function to extract and use the **proper velocity axis parameters** from the FITS file. The datacube now has a fully WCS-compliant velocity axis instead of simple channel indices.

### What Was Changed

**File Modified:** `src/oi_zeigt/mapping/gridding.py`

#### 1. New Helper Function: `_get_spectral_axis_params()`

A new function was added to extract spectral axis parameters from FITS table columns:

```python
def _get_spectral_axis_params(hdul: fits.HDUList) -> Tuple[float, float, float, str, float]:
    """
    Extract spectral axis parameters from FITS file for velocity axis reconstruction.
    
    Returns:
    - velo_ref: Reference velocity in m/s (from VELOCITY column)
    - deltav: Velocity spacing per channel in m/s (from DELTAV column)
    - restfreq: Rest frequency in Hz (from RESTFREQ column)
    - veldef: Velocity definition string (from VELDEF column)
    - nchans: Number of spectral channels
    """
```

**Parameters extracted:**
- **VELOCITY**: Reference velocity (470,000 m/s for M51)
- **DELTAV**: Velocity step per channel (500 m/s for M51)
- **VELDEF**: Velocity definition ("RADI-LSR" for radial velocity in Local Standard of Rest)
- **RESTFREQ**: Rest frequency of the line (1.9005369×10¹² Hz for M51)
- **SPECTRUM shape**: Number of channels (1,264 for M51)

#### 2. Modified `create_spectral_datacube()` Function

The function now:
1. **Extracts spectral parameters** from the FITS table using `_get_spectral_axis_params()`
2. **Prints spectral axis information** for user feedback
3. **Uses proper velocity spacing** when creating the WCS header
4. **Sets CDELT3 to the actual velocity step** (not 1.0 like before!)

#### 3. Updated WCS Header Creation

The key fix in WCS header creation:

**Before:**
```python
cdelt3=1.0,  # ← WRONG: Just channel indices!
```

**After:**
```python
cdelt3=deltav,  # ← CORRECT: Actual velocity spacing from FITS (500 m/s for M51)
```

### Velocity Axis Reconstruction

The velocity of each channel is now properly calculated as:

$$v(i) = \text{CRVAL3} + i \times \text{CDELT3}$$

Where:
- **CRVAL3** = 470,000 m/s (reference velocity)
- **CDELT3** = 500 m/s (velocity per channel)
- **i** = channel index (0 to 1,263)

**Example velocities for M51:**
- Channel 0: 470,000 m/s (470.0 km/s)
- Channel 1: 470,500 m/s (470.5 km/s)
- Channel 1,263: 1,101,500 m/s (1,101.5 km/s)
- **Total span: 632 km/s**

### Test Results

A test was run on the M51 data:

```
Filtered to 39,200 observations matching 'M51'
Spectral axis parameters:
  Reference velocity: 470000.00 m/s (470.00 km/s)
  Velocity step: 500.00 m/s (0.5000 km/s)
  Velocity definition: RADI-LSR
  Rest frequency: 1.9005369000e+12 Hz
  Total velocity range: 632000.00 m/s (632.00 km/s)

Datacube created: 1264 channels × 12 × 11 pixels

WCS Header Output:
  CTYPE3: VRAD (Radial velocity)
  CUNIT3: m/s (meters per second)
  CRVAL3: 470000.00 m/s (470.00 km/s)
  CRPIX3: 1.0 (reference pixel)
  CDELT3: 500.00 m/s (0.5000 km/s per channel)
```

✅ **Test completed successfully!**

### Advantages

1. **WCS Compliance**: Output FITS files are now fully WCS-compliant
2. **Scientific Accuracy**: Velocities properly calibrated to LSR frame
3. **Tool Compatibility**: Can be read by CASA, spectral-cube, and other radio astronomy tools
4. **Self-Documenting**: Velocity information stored in FITS headers
5. **Backward Compatible**: Function signature and behavior unchanged

### Files Modified/Created

1. **Modified:** `src/oi_zeigt/mapping/gridding.py`
   - Added `_get_spectral_axis_params()` function
   - Updated `create_spectral_datacube()` to extract and use velocity parameters
   - Updated docstring to document velocity axis reconstruction

2. **Created:** `test_velocity_axis.py`
   - Test script demonstrating the new functionality
   - Verifies velocity axis parameters are correctly extracted

3. **Created:** `VELOCITY_AXIS_RECONSTRUCTION.md`
   - Detailed documentation of the changes
   - Velocity axis reconstruction formula
   - Usage examples and future enhancements

### How to Use

The usage is **identical** to before. Just use the CLI command as usual:

```bash
create_datacube --config config.toml --output my_datacube.fits
```

Or programmatically:

```python
from oi_zeigt.mapping.gridding import create_spectral_datacube
from astropy.io import fits

hdul = fits.open('your_data.fits')
datacube, wcs_header, fig = create_spectral_datacube(
    hdul,
    beamsize_deg=0.05,
    pixsize=0.01,
    output_file='datacube.fits'
)

# Now the velocity axis is properly calibrated!
print(f"Velocity per channel: {wcs_header['CDELT3']} {wcs_header['CUNIT3']}")
```

### What's Next?

The datacube now has:
- ✅ Proper spatial axes (RA, Dec) - already implemented
- ✅ Proper velocity axis (with correct scaling) - **just implemented**
- Ready for advanced analysis (moment maps, spectral line fitting, etc.)

Would you like to:
1. **Save a test datacube** with the proper velocity axis?
2. **Analyze the datacube** (create moment maps, integrated intensity maps)?
3. **Perform spectral analysis** (line fitting, velocity channel selection)?
