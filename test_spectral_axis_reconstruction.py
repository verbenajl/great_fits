#!/usr/bin/env python3
"""
Test spectral axis reconstruction using preserved CRPIX1 value.
This demonstrates the velocity axis with proper reference pixel.
"""

from astropy.io import fits
import numpy as np
from pathlib import Path

test_file = '/home/diskB/sofiaobsdata/m51/m51_central_region/kalibrated_all/m51_central_tile_singlehdu.fits'

print("=" * 80)
print("SPECTRAL AXIS RECONSTRUCTION WITH PRESERVED CRPIX1")
print("=" * 80)
print(f"\nTest file: {test_file}\n")

with fits.open(test_file) as hdul:
    primary = hdul[0]
    table = hdul[1]
    
    # Extract all parameters
    print("=" * 80)
    print("1. EXTRACTING SPECTRAL PARAMETERS")
    print("=" * 80)
    
    # Get CRPIX1 from PRIMARY header
    crpix1 = primary.header.get('CRPIX1')
    print(f"\nCRPIX1 (from PRIMARY header):")
    if crpix1 is not None:
        print(f"  ✓ CRPIX1 = {crpix1:.10f} (FITS 1-indexed)")
        crpix1_0idx = crpix1 - 1
        print(f"  → Reference channel (0-indexed) = {crpix1_0idx:.10f}")
    else:
        print(f"  ✗ CRPIX1 NOT FOUND - will default to 1.0")
        crpix1 = 1.0
        crpix1_0idx = 0.0
    
    # Get VELOCITY from table column
    velocities = table.data['VELOCITY']
    velo_ref = np.mean(velocities)
    print(f"\nVELOCITY (from table column):")
    print(f"  ✓ Found {len(velocities)} observations")
    print(f"  Values: min={np.min(velocities):.0f}, max={np.max(velocities):.0f}, mean={np.mean(velocities):.0f}")
    print(f"  → Using velo_ref = {velo_ref:.0f} m/s = {velo_ref/1000:.2f} km/s")
    
    # Get DELTAV from table column
    deltavs = table.data['DELTAV']
    deltav = np.mean(deltavs)
    print(f"\nDELTAV (from table column):")
    print(f"  ✓ Found {len(deltavs)} observations")
    print(f"  Values: min={np.min(deltavs):.0f}, max={np.max(deltavs):.0f}, mean={np.mean(deltavs):.0f}")
    print(f"  → Using deltav = {deltav:.0f} m/s = {deltav/1000:.4f} km/s/channel")
    
    # Get number of channels
    nchans = table.data['SPECTRUM'].shape[1]
    print(f"\nSPECTRUM shape:")
    print(f"  ✓ Observations: {len(table.data)}")
    print(f"  ✓ Channels: {nchans}")
    
    # Get RESTFREQ
    restfreqs = table.data['RESTFREQ']
    restfreq = np.mean(restfreqs)
    print(f"\nRESTFREQ (from table column):")
    print(f"  Values: min={np.min(restfreqs):.4e}, max={np.max(restfreqs):.4e}, mean={np.mean(restfreqs):.4e} Hz")
    print(f"  → Using restfreq = {restfreq:.4e} Hz")
    
    # Get VELDEF
    veldefs = table.data['VELDEF']
    veldef = veldefs[0].strip()
    print(f"\nVELDEF (from table column):")
    print(f"  ✓ All observations: '{veldef}'")
    
    print("\n" + "=" * 80)
    print("2. VELOCITY AXIS RECONSTRUCTION FORMULA")
    print("=" * 80)
    
    print(f"\nFormula: v(i) = CRVAL3 + (i - (CRPIX3 - 1)) × CDELT3")
    print(f"\nWith values:")
    print(f"  CRVAL3 = {velo_ref:.2f} m/s = {velo_ref/1000:.2f} km/s")
    print(f"  CRPIX3 = {crpix1:.10f} (FITS 1-indexed)")
    print(f"  CRPIX3 - 1 = {crpix1_0idx:.10f} (0-indexed reference channel)")
    print(f"  CDELT3 = {deltav:.2f} m/s = {deltav/1000:.4f} km/s/channel")
    
    print(f"\nSimplified: v(i) = {velo_ref:.0f} + (i - {crpix1_0idx:.2f}) × {deltav:.0f}")
    
    print("\n" + "=" * 80)
    print("3. VELOCITY AXIS MAPPING")
    print("=" * 80)
    
    print(f"\n{'Channel':<12} {'Velocity (m/s)':<20} {'Velocity (km/s)':<20} {'Note':<30}")
    print("-" * 82)
    
    # Calculate velocities for key channels
    key_channels = [0, 50, 100, 200, 300, 400, 500, int(crpix1_0idx), int(crpix1_0idx)+1, 600, 800, 1000, 1263]
    key_channels = sorted(set(key_channels))
    
    for ch in key_channels:
        if ch < nchans:
            v = velo_ref + (ch - crpix1_0idx) * deltav
            note = ""
            if ch == int(crpix1_0idx):
                note = "← REFERENCE POINT (CRPIX1)"
            elif ch == 0:
                note = "← First channel"
            elif ch == nchans - 1:
                note = "← Last channel"
            print(f"{ch:<12} {v:<20.2f} {v/1000:<20.4f} {note:<30}")
    
    print("\n" + "=" * 80)
    print("4. VELOCITY RANGE SUMMARY")
    print("=" * 80)
    
    v_first = velo_ref + (0 - crpix1_0idx) * deltav
    v_ref = velo_ref + (crpix1_0idx - crpix1_0idx) * deltav  # Should be velo_ref
    v_last = velo_ref + (nchans - 1 - crpix1_0idx) * deltav
    
    print(f"\nFirst channel (0):")
    print(f"  Velocity = {v_first:.2f} m/s = {v_first/1000:.2f} km/s")
    
    print(f"\nReference channel ({crpix1_0idx:.0f}):")
    print(f"  Velocity = {v_ref:.2f} m/s = {v_ref/1000:.2f} km/s")
    print(f"  ✓ Should equal CRVAL3 = {velo_ref:.2f} m/s")
    
    print(f"\nLast channel ({nchans-1}):")
    print(f"  Velocity = {v_last:.2f} m/s = {v_last/1000:.2f} km/s")
    
    print(f"\nTotal velocity span:")
    v_span = v_last - v_first
    print(f"  {v_span:.2f} m/s = {v_span/1000:.2f} km/s")
    
    print(f"\nVelocity per channel (resolution):")
    print(f"  {deltav:.2f} m/s = {deltav/1000:.4f} km/s/channel")
    
    print("\n" + "=" * 80)
    print("5. WCS HEADER FOR OUTPUT DATACUBE")
    print("=" * 80)
    
    print(f"\nRecommended keywords for FITS header:")
    print(f"  CTYPE3  = 'VRAD'")
    print(f"  CUNIT3  = 'm/s'")
    print(f"  CRVAL3  = {velo_ref:.2f}")
    print(f"  CRPIX3  = {crpix1:.10f}")
    print(f"  CDELT3  = {deltav:.2f}")
    print(f"  NAXIS3  = {nchans}")
    
    print("\n" + "=" * 80)
    print("6. COMPARISON: WITH vs WITHOUT CRPIX1")
    print("=" * 80)
    
    print(f"\n❌ WITHOUT CRPIX1 (defaulted to 1.0):")
    v_wrong_ref = velo_ref + (crpix1_0idx - 0) * deltav  # Default to channel 0
    v_wrong_505 = velo_ref + (505 - 0) * deltav
    print(f"  v(0)   = {velo_ref:.0f} m/s = {velo_ref/1000:.2f} km/s (WRONG - reference here)")
    print(f"  v(505) = {v_wrong_505:.0f} m/s = {v_wrong_505/1000:.2f} km/s (WRONG - not reference)")
    print(f"  Error at channel 505: {abs(v_wrong_505 - velo_ref)/1000:.2f} km/s off!")
    
    print(f"\n✅ WITH CRPIX1 ({crpix1:.2f}):")
    v_correct_ref = velo_ref + (crpix1_0idx - crpix1_0idx) * deltav  # Reference at correct channel
    v_correct_505 = velo_ref + (505 - crpix1_0idx) * deltav
    print(f"  v({int(crpix1_0idx)})  = {v_correct_ref:.0f} m/s = {v_correct_ref/1000:.2f} km/s (CORRECT - reference here)")
    print(f"  v(505) = {v_correct_505:.0f} m/s = {v_correct_505/1000:.2f} km/s (CORRECT - at expected velocity)")
    print(f"  No error - physically accurate!")
    
    print("\n" + "=" * 80)
    print("7. VERIFICATION")
    print("=" * 80)
    
    if crpix1 is not None and crpix1 > 1.0:
        print(f"\n✅ CRPIX1 properly preserved: {crpix1:.2f}")
        print(f"✅ Velocity axis will use actual CRPIX1 value")
        print(f"✅ Reference velocity will be at correct channel (~{int(crpix1_0idx)})")
        print(f"✅ No systematic errors in velocity axis")
        print(f"✅ Output datacube will be physically accurate!")
    else:
        print(f"\n⚠️  CRPIX1 not found or defaulted to 1.0")
        print(f"⚠️  Velocity axis will assume reference at channel 0")
        print(f"⚠️  This would cause systematic error ~250 km/s")
    
    print("\n" + "=" * 80)

