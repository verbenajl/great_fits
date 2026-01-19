#!/usr/bin/env python3
"""
Documentation: Velocity Axis Reconstruction in Spectral Datacube

This document explains the modifications made to create_spectral_datacube()
to properly reconstruct the velocity axis from FITS table parameters.
"""

# ==============================================================================
# SUMMARY OF CHANGES
# ==============================================================================

"""
Modified Files:
  - src/oi_zeigt/mapping/gridding.py

Changes Made:
  1. Added new helper function: _get_spectral_axis_params()
  2. Modified create_spectral_datacube() to extract and use velocity parameters
  3. Updated WCS header creation to use proper velocity spacing

Key Improvements:
  ✓ Velocity axis is now a proper WCS axis (not just channel indices)
  ✓ Velocity information extracted from FITS table columns
  ✓ Proper CRVAL3/CDELT3 values in output FITS headers
"""

# ==============================================================================
# NEW HELPER FUNCTION: _get_spectral_axis_params()
# ==============================================================================

"""
Location: src/oi_zeigt/mapping/gridding.py (lines ~305-378)

Purpose:
  Extract spectral axis parameters from FITS table columns to reconstruct
  the velocity axis in the output datacube.

Parameters extracted from FITS table:
  1. VELOCITY column
     - Reference velocity in m/s
     - LSR (Local Standard of Rest) convention
     - Example: 470,000 m/s = 470 km/s
  
  2. DELTAV column
     - Velocity spacing per spectral channel in m/s
     - Example: 500 m/s = 0.5 km/s per channel
  
  3. VELDEF column
     - Velocity definition string
     - Example: 'RADI-LSR' (Radial velocity in LSR)
  
  4. RESTFREQ column
     - Rest frequency of the observed line in Hz
     - Example: 1.9005369e12 Hz (CO 1→0 transition)
  
  5. SPECTRUM column shape
     - Number of spectral channels
     - Example: 1264 channels

Returns:
  tuple: (velo_ref, deltav, restfreq, veldef, nchans)
  - velo_ref: Reference velocity (m/s)
  - deltav: Velocity spacing per channel (m/s)
  - restfreq: Rest frequency (Hz)
  - veldef: Velocity definition string
  - nchans: Number of spectral channels
"""

# ==============================================================================
# VELOCITY AXIS RECONSTRUCTION FORMULA
# ==============================================================================

"""
The velocity of spectral channel i is calculated as:

    v(i) = CRVAL3 + (i - CRPIX3 + 1) * CDELT3

Where:
  - v(i): Velocity of channel i in m/s
  - CRVAL3: Reference velocity from VELOCITY column (m/s)
  - i: Channel index (0-based)
  - CRPIX3: Reference pixel in FITS convention (1-based), usually 1.0
  - CDELT3: Velocity spacing from DELTAV column (m/s)

For CRPIX3 = 1.0 (standard convention):
    v(i) = CRVAL3 + i * CDELT3

Example (M51 data):
  - CRVAL3 = 470,000 m/s (470 km/s)
  - CDELT3 = 500 m/s (0.5 km/s)
  - Channel 0: v = 470,000 m/s
  - Channel 1: v = 470,500 m/s
  - Channel 1263: v = 1,101,500 m/s
  - Total span: 632 km/s
"""

# ==============================================================================
# WCS HEADER KEYWORDS (BEFORE AND AFTER)
# ==============================================================================

"""
BEFORE (Old Implementation):
  CTYPE3 = 'VRAD'
  CUNIT3 = 'm/s'
  CRVAL3 = 0.0 (from header.get('VELO-LSR', 0.0), often missing)
  CRPIX3 = 1.0
  CDELT3 = 1.0  ← THIS WAS WRONG! (channel indices, not velocity)

AFTER (New Implementation):
  CTYPE3 = 'VRAD'
  CUNIT3 = 'm/s'
  CRVAL3 = 470000.0  ← From VELOCITY column
  CRPIX3 = 1.0
  CDELT3 = 500.0  ← From DELTAV column (m/s per channel)

The WCS library can now correctly convert between:
  - Pixel coordinates (channel indices)
  - Velocity coordinates (m/s in LSR frame)
"""

# ==============================================================================
# TEST RESULTS
# ==============================================================================

"""
Test file: test_velocity_axis.py

Sample output from M51 observations:
  Input: 39,200 M51 observations, 1,264 spectral channels each
  Output datacube: 1264 channels × 12 pixels × 11 pixels

Extracted Parameters:
  Reference velocity: 470,000 m/s (470.00 km/s)
  Velocity step: 500 m/s (0.5000 km/s)
  Velocity definition: RADI-LSR
  Rest frequency: 1.9005369e12 Hz
  Total velocity range: 632,000 m/s (632.00 km/s)

WCS Header Output:
  CTYPE3: VRAD
  CUNIT3: m/s
  CRVAL3: 470000.00 m/s
  CRPIX3: 1.0
  CDELT3: 500.00 m/s
  
  First channel: 470,000 m/s (470.00 km/s)
  Last channel: 1,101,500 m/s (1101.50 km/s)
  Total span: 631,500 m/s (631.50 km/s)
"""

# ==============================================================================
# ADVANTAGES OF PROPER VELOCITY AXIS
# ==============================================================================

"""
1. WCS Compliance
   - Output FITS files are now fully WCS-compliant
   - Can be read by CASA, spectral-cube, and other tools
   - Automatic coordinate transformations work correctly

2. Scientific Accuracy
   - Velocity coordinates are properly calibrated to LSR frame
   - Consistent with original SOFIA observations
   - Supports velocity range queries and slicing

3. Interoperability
   - Can stack datacubes from different observations
   - Velocity axes are aligned across multiple datasets
   - Compatible with standard radio astronomy tools

4. Documentation
   - Velocity parameters stored in FITS header
   - Self-documenting format for data analysis
   - Traceable back to original observations
"""

# ==============================================================================
# HOW TO USE THE MODIFIED FUNCTION
# ==============================================================================

"""
Usage is identical to before:

    from oi_zeigt.mapping.gridding import create_spectral_datacube
    from astropy.io import fits
    
    # Open FITS file
    hdul = fits.open('your_data.fits')
    
    # Create datacube with proper velocity axis
    datacube, wcs_header, fig = create_spectral_datacube(
        hdul,
        beamsize_deg=0.05,
        pixsize=0.01,
        object_filter='M51',
        output_file='m51_datacube.fits'
    )
    
    # The datacube now has a proper velocity axis in the WCS header!
    print(wcs_header['CRVAL3'], wcs_header['CDELT3'], wcs_header['CUNIT3'])
    # Output: 470000.0 500.0 m/s
"""

# ==============================================================================
# BACKWARD COMPATIBILITY
# ==============================================================================

"""
✓ The changes are fully backward compatible:
  - Function signature unchanged
  - Default behavior unchanged
  - Existing code continues to work
  - Output format compatible with old datacubes
  - Only the velocity axis is now properly calibrated
"""

# ==============================================================================
# FUTURE ENHANCEMENTS
# ==============================================================================

"""
Potential improvements:
  1. Add velocity axis conversion utilities
  2. Support for frequency axis (CTYPE3 = 'FREQ')
  3. Doppler shift calculation tools
  4. Velocity channel selection/slicing utilities
  5. Multi-line observations with different rest frequencies
"""

print(__doc__)
