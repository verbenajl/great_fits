# Adding Spectral Axis Information to Datacube - Options & Best Practices

## Current State

**What we have:**
- ✅ Datacube shape: (1264 channels, 12 dec pixels, 11 ra pixels)
- ✅ WCS header with proper CRPIX3, CDELT3, CRVAL3 keywords
- ✅ Users can reconstruct velocity using WCS library: `wcs.pixel_to_world_values()`

**What's missing:**
- ❌ Explicit velocity array (array of velocities for each channel)
- ❌ Easy access for users who don't use WCS
- ❌ Direct velocity information in the FITS structure

---

## Option 1: Add Spectral Axis as Binary Table Extension (⭐ RECOMMENDED)

### What it does:
Creates a new FITS Binary Table extension containing:
- Channel number (0-1263)
- Velocity (in m/s)
- Frequency (in Hz, if needed)

### Advantages:
- ✅ Non-intrusive (doesn't modify datacube data)
- ✅ Standard FITS format
- ✅ Easy to read with any FITS tool
- ✅ Metadata (units, descriptions) included
- ✅ Compatible with spectral-cube library
- ✅ Can be extended with more information later

### Disadvantages:
- ⚠️ Adds extra HDU to FITS file
- ⚠️ File size increases slightly
- ⚠️ Users need to know to look for it

### Example FITS Structure:
```
HDU 0 (PRIMARY):
  - Datacube array (1264, 12, 11)
  - WCS header

HDU 1 (BINTABLE - SPECTRUM):
  - CHANNEL: array of channel indices (0-1263)
  - VELOCITY: array of velocities in m/s
  - FREQUENCY: array of frequencies in Hz
  - Comments: CRPIX1, RESTFREQ, VELDEF, etc.
```

---

## Option 2: Add Spectral Axis as FITS Table Keywords

### What it does:
Adds extra keywords to PRIMARY header:
- VELAXIS: Array of velocities (if possible)
- Or structured keywords: VEL0001, VEL0002, ...

### Advantages:
- ✅ Everything in one HDU
- ✅ No extra complexity

### Disadvantages:
- ❌ FITS keywords have 80-char limit
- ❌ Can't store 1264 velocity values as keywords
- ❌ Very inefficient
- ❌ Not recommended by FITS standard

### Verdict: **NOT RECOMMENDED** - Too limiting

---

## Option 3: Save Spectral Axis to Separate File

### What it does:
Create a companion file with velocity information:
```
datacube.fits          → Main datacube
datacube_spec.txt      → Text file with channel | velocity | frequency
```

### Advantages:
- ✅ Simple, human-readable
- ✅ Doesn't modify main FITS file
- ✅ Easy to parse in scripts

### Disadvantages:
- ❌ Separate file to manage
- ❌ Easy to lose/forget
- ❌ Not automatically available with datacube

### Verdict: **OK as backup**, but not primary solution

---

## Option 4: Use spectral-cube Library Format

### What it does:
Creates a FITS file compatible with astropy's `spectral-cube` library:
- Proper WCS in PRIMARY header ✅ (we already have this)
- Optionally add BEAMS HDU for beam information
- Optionally add per-channel metadata

### Advantages:
- ✅ Compatible with popular radio astronomy tools
- ✅ Automatic velocity axis support
- ✅ Advanced features (moment maps, slicing)
- ✅ Industry standard

### Disadvantages:
- ⚠️ Requires WCS to be set up correctly (we have this!)
- ⚠️ More complex FITS structure

### Verdict: **GOOD** - but focuses on WCS interpretation, not explicit axis storage

---

## Recommendation: Hybrid Approach ⭐⭐⭐

### Implement **Option 1 (Binary Table)** + leverage **Option 4 (spectral-cube)**

**Why this works:**
1. **Primary method**: Use proper WCS header (already done!)
   - Users with spectral-cube can do: `cube.spectral_axis` 
   - Users with astropy.wcs can do: `wcs.pixel_to_world_values()`
   
2. **Convenience method**: Add Binary Table with explicit velocities
   - For users who want direct access
   - For tools that don't use WCS
   - For verification/debugging

**Implementation Steps:**
1. Keep datacube as PRIMARY HDU with WCS header (no changes needed!)
2. Add SPECTRUM Binary Table as HDU 1 with:
   - CHANNEL column (channel index)
   - VELOCITY column (m/s, from FITS)
   - FREQUENCY column (Hz, calculated from velocity)
3. Add header comments in SPECTRUM HDU explaining the axis

---

## How Each Tool Would Use It

### Using `spectral-cube` library:
```python
from spectral_cube import SpectralCube
cube = SpectralCube.read('datacube.fits')
spec_axis = cube.spectral_axis  # Automatically extracted from WCS!
print(spec_axis[505])  # Get velocity at channel 505
```

### Using `astropy.wcs`:
```python
from astropy.wcs import WCS
from astropy.io import fits

hdul = fits.open('datacube.fits')
wcs = WCS(hdul[0].header)
velocities = wcs.pixel_to_world_values(np.arange(1264), 0, 0)[2]
print(velocities[505])  # Channel 505
```

### Using direct table access:
```python
from astropy.io import fits

hdul = fits.open('datacube.fits')
spec_table = hdul['SPECTRUM'].data  # New table HDU
velocities = spec_table['VELOCITY']  # Direct velocity array
print(velocities[505])  # Channel 505
```

---

## Code Implementation Preview

### What we'd add to `create_spectral_datacube()`:

```python
# After creating WCS header and before saving FITS

# Create spectral axis table
from astropy.table import Table

# Calculate velocities for all channels
channels = np.arange(nvel)
velocities = velo_ref + (channels - (crpix1_spec - 1)) * deltav
frequencies = restfreq_from_table + (velocities / 3e8) * restfreq_from_table

# Create table
spec_table = Table()
spec_table['CHANNEL'] = channels
spec_table['VELOCITY'] = velocities * u.m / u.s
spec_table['FREQUENCY'] = frequencies * u.Hz

# Create Binary Table HDU
spec_bintable = fits.BinTableHDU(spec_table, name='SPECTRUM')
spec_bintable.header['COMMENT'] = 'Spectral axis information'
spec_bintable.header['CRPIX1'] = crpix1_spec
spec_bintable.header['RESTFREQ'] = restfreq_from_table
spec_bintable.header['VELDEF'] = veldef

# Save with both PRIMARY and SPECTRUM HDUs
hdul = fits.HDUList([
    fits.PrimaryHDU(data=datacube, header=wcs_header),
    spec_bintable
])
hdul.writeto(output_file, overwrite=True)
```

---

## Decision Points Before Implementation

**Question 1: Should we add the Binary Table?**
- [ ] Yes - Add explicit spectral axis table for easy access
- [ ] No - Rely on WCS only (current state)
- [ ] Maybe - Add as optional parameter to `create_spectral_datacube()`

**Question 2: Should we add FREQUENCY column?**
- [ ] Yes - Include both velocity and frequency
- [ ] No - Just velocity (simpler, smaller file)

**Question 3: Should we add metadata to the table?**
- [ ] Yes - Include CRPIX1, RESTFREQ, VELDEF, etc. as keywords
- [ ] No - Keep table simple

**Question 4: Should we make this compatible with spectral-cube?**
- [ ] Yes - Follow spectral-cube conventions
- [ ] No - Just add the table without special structure

---

## My Recommendation

✅ **Add Binary Table extension with:**
- CHANNEL column (for reference)
- VELOCITY column (required, m/s)
- FREQUENCY column (optional, Hz)
- Metadata in header keywords

✅ **Make it optional:**
- Add parameter `include_spectral_axis=True` to function
- Users can disable if they want smaller files

✅ **Document clearly:**
- Add comments about how to access the data
- Show examples for spectral-cube and astropy.wcs

---

## Would you like to proceed with this approach?

If yes, I can:
1. ✅ Show the exact code changes needed
2. ✅ Implement the Binary Table addition
3. ✅ Test it with real M51 data
4. ✅ Update documentation
5. ✅ Create examples for different use cases

Let me know your thoughts!
