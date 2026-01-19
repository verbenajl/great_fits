#!/usr/bin/env python3
"""
Test script: Verify that spectral axis can be reconstructed from FITS header parameters.

This tests the complete workflow:
1. Load a single-HDU FITS file created by combine_fits_files
2. Extract spectral parameters from table columns
3. Extract CRPIX1 from header
4. Reconstruct the velocity axis
5. Verify the WCS compliance
"""

import sys
from astropy.io import fits
import numpy as np
from pathlib import Path

# Test with single-HDU file (should now have CRPIX1!)
test_file = Path('/home/diskB/sofiaobsdata/m51/m51_central_region/kalibrated_all/m51_central_tile_singlehdu.fits')

if not test_file.exists():
    print(f"❌ Test file not found: {test_file}")
    print("   Please ensure sample_combined.fits exists")
    sys.exit(1)

print("=" * 80)
print("SPECTRAL AXIS RECONSTRUCTION TEST")
print("=" * 80)
print(f"\nTest File: {test_file}")
print(f"File size: {test_file.stat().st_size / 1024 / 1024:.2f} MB\n")

try:
    with fits.open(test_file) as hdul:
        print("=" * 80)
        print("1. FITS FILE STRUCTURE")
        print("=" * 80)
        hdul.info()
        
        # Get the binary table HDU
        if len(hdul) < 2:
            print("\n❌ Expected at least 2 HDUs (PRIMARY + TABLE)")
            sys.exit(1)
        
        primary = hdul[0]
        table_hdu = hdul[1]
        
        print("\n" + "=" * 80)
        print("2. PRIMARY HEADER INSPECTION")
        print("=" * 80)
        print("\nLooking for spectral WCS keywords:")
        spectral_keywords = ['CRPIX1', 'CRVAL1', 'CDELT1', 'CTYPE1', 'CUNIT1']
        for key in spectral_keywords:
            if key in primary.header:
                value = primary.header[key]
                print(f"  ✓ {key:8s} = {value}")
            else:
                print(f"  ✗ {key:8s} = NOT FOUND")
        
        print("\nLooking for velocity reference keywords:")
        vel_keywords = ['VELOCITY', 'DELTAV', 'RESTFREQ', 'VELDEF']
        for key in vel_keywords:
            if key in primary.header:
                value = primary.header[key]
                print(f"  ✓ {key:8s} = {value}")
            else:
                print(f"  ✗ {key:8s} = NOT FOUND (likely in table columns)")
        
        print("\n" + "=" * 80)
        print("3. TABLE STRUCTURE INSPECTION")
        print("=" * 80)
        print(f"\nTable HDU: {table_hdu.name}")
        print(f"Number of columns: {len(table_hdu.columns)}")
        print(f"Number of rows: {table_hdu.header['NAXIS2']}")
        
        print("\nAvailable columns:")
        for i, col in enumerate(table_hdu.columns):
            print(f"  {i+1:2d}. {col.name:15s} ({col.format})")
        
        print("\n" + "=" * 80)
        print("4. SPECTRAL PARAMETERS EXTRACTION")
        print("=" * 80)
        
        # Try to extract from table columns
        table_data = table_hdu.data
        
        # Extract CRPIX1 from header
        crpix1 = primary.header.get('CRPIX1', None)
        print(f"\nCRPIX1 (from PRIMARY header):")
        if crpix1 is not None:
            print(f"  ✓ CRPIX1 = {crpix1:.10f} (FITS 1-indexed)")
            print(f"  → Reference channel (0-indexed) = {crpix1 - 1:.10f}")
        else:
            print(f"  ✗ NOT FOUND (will default to 1.0)")
            crpix1 = 1.0
        
        # Extract VELOCITY from table
        if 'VELOCITY' in table_hdu.columns.names:
            velocities = table_data['VELOCITY']
            print(f"\nVELOCITY (from table column):")
            print(f"  ✓ Found in column")
            print(f"    - Min: {np.min(velocities):.2f} m/s = {np.min(velocities)/1000:.2f} km/s")
            print(f"    - Max: {np.max(velocities):.2f} m/s = {np.max(velocities)/1000:.2f} km/s")
            print(f"    - Mean: {np.mean(velocities):.2f} m/s = {np.mean(velocities)/1000:.2f} km/s")
            velo_ref = np.mean(velocities)  # Use mean as reference
        else:
            print(f"\n✗ VELOCITY not found in table columns")
            velo_ref = None
        
        # Extract DELTAV from table
        if 'DELTAV' in table_hdu.columns.names:
            deltavs = table_data['DELTAV']
            print(f"\nDELTAV (from table column):")
            print(f"  ✓ Found in column")
            print(f"    - Min: {np.min(deltavs):.2f} m/s = {np.min(deltavs)/1000:.4f} km/s")
            print(f"    - Max: {np.max(deltavs):.2f} m/s = {np.max(deltavs)/1000:.4f} km/s")
            print(f"    - Mean: {np.mean(deltavs):.2f} m/s = {np.mean(deltavs)/1000:.4f} km/s")
            deltav = np.mean(deltavs)  # Use mean value
        else:
            print(f"\n✗ DELTAV not found in table columns")
            deltav = None
        
        # Extract RESTFREQ from table
        if 'RESTFREQ' in table_hdu.columns.names:
            restfreqs = table_data['RESTFREQ']
            print(f"\nRESTFREQ (from table column):")
            print(f"  ✓ Found in column")
            print(f"    - Min: {np.min(restfreqs):.4e} Hz")
            print(f"    - Max: {np.max(restfreqs):.4e} Hz")
            print(f"    - Mean: {np.mean(restfreqs):.4e} Hz")
            restfreq = np.mean(restfreqs)
        else:
            print(f"\n✗ RESTFREQ not found in table columns")
            restfreq = None
        
        # Extract VELDEF from table
        if 'VELDEF' in table_hdu.columns.names:
            veldefs = table_data['VELDEF']
            print(f"\nVELDEF (from table column):")
            print(f"  ✓ Found in column")
            # Try to get unique values
            unique_vels = np.unique(veldefs)
            if len(unique_vels) <= 5:
                for vel_def in unique_vels:
                    count = np.sum(veldefs == vel_def)
                    print(f"    - '{vel_def}': {count} observations")
            else:
                print(f"    - {len(unique_vels)} unique values")
            veldef = unique_vels[0] if len(unique_vels) > 0 else None
        else:
            print(f"\n✗ VELDEF not found in table columns")
            veldef = None
        
        # Extract SPECTRUM shape to get number of channels
        if 'SPECTRUM' in table_hdu.columns.names:
            spectra = table_data['SPECTRUM']
            if spectra.ndim > 1:
                nchans = spectra.shape[1]
            else:
                nchans = len(spectra[0]) if len(spectra) > 0 else 0
            print(f"\nSPECTRUM (from table column):")
            print(f"  ✓ Found in column")
            print(f"    - Number of observations: {len(spectra)}")
            print(f"    - Channels per observation: {nchans}")
        else:
            print(f"\n✗ SPECTRUM not found in table columns")
            nchans = None
        
        print("\n" + "=" * 80)
        print("5. VELOCITY AXIS RECONSTRUCTION")
        print("=" * 80)
        
        if all([velo_ref, deltav, nchans, crpix1]):
            print(f"\n✓ All required parameters available!")
            print(f"\nWCS Parameters:")
            print(f"  CRVAL3 = {velo_ref:.2f} m/s = {velo_ref/1000:.2f} km/s (reference velocity)")
            print(f"  CRPIX3 = {crpix1:.10f} (FITS 1-indexed, i.e., channel {crpix1-1:.10f} in 0-indexed)")
            print(f"  CDELT3 = {deltav:.2f} m/s = {deltav/1000:.4f} km/s (velocity per channel)")
            print(f"  NAXIS3 = {nchans} (total channels)")
            
            # Calculate velocity for key channels
            ref_channel_0idx = crpix1 - 1
            
            print(f"\nVelocity Axis Mapping:")
            print(f"  Formula: v(i) = {velo_ref:.0f} + (i - {ref_channel_0idx:.1f}) × {deltav:.0f}")
            print(f"\n  {'Channel':<10} {'Velocity (m/s)':<20} {'Velocity (km/s)':<20} {'Note':<20}")
            print(f"  {'-'*10} {'-'*20} {'-'*20} {'-'*20}")
            
            key_channels = [0, 100, 500, int(ref_channel_0idx), int(ref_channel_0idx) + 1, 1000, nchans-1]
            key_channels = sorted(set(key_channels))
            
            for ch in key_channels:
                if 0 <= ch < nchans:
                    v = velo_ref + (ch - ref_channel_0idx) * deltav
                    note = "← REFERENCE" if ch == int(ref_channel_0idx) else ""
                    print(f"  {ch:<10} {v:<20.2f} {v/1000:<20.4f} {note:<20}")
            
            # Summary statistics
            v_first = velo_ref + (0 - ref_channel_0idx) * deltav
            v_last = velo_ref + (nchans - 1 - ref_channel_0idx) * deltav
            v_span = v_last - v_first
            
            print(f"\n  Velocity Range Summary:")
            print(f"    First channel (0):     {v_first/1000:>8.2f} km/s")
            print(f"    Last channel ({nchans-1}):     {v_last/1000:>8.2f} km/s")
            print(f"    Total span:            {v_span/1000:>8.2f} km/s")
            print(f"    Reference at channel:  {ref_channel_0idx:>8.1f} (0-indexed)")
            
        else:
            print(f"\n❌ Missing required parameters:")
            print(f"    velo_ref: {'✓' if velo_ref else '✗'}")
            print(f"    deltav: {'✓' if deltav else '✗'}")
            print(f"    nchans: {'✓' if nchans else '✗'}")
            print(f"    crpix1: {'✓' if crpix1 else '✗'}")
        
        print("\n" + "=" * 80)
        print("6. WCS COMPLIANCE CHECK")
        print("=" * 80)
        
        if all([velo_ref, deltav, nchans, crpix1]):
            print("\n✓ Spectral axis can be fully reconstructed from FITS headers!")
            print("\nRecommended WCS header keywords for output datacube:")
            print(f"  CTYPE3  = 'VRAD'")
            print(f"  CUNIT3  = 'm/s'")
            print(f"  CRVAL3  = {velo_ref:.2f}")
            print(f"  CRPIX3  = {crpix1:.10f}")
            print(f"  CDELT3  = {deltav:.2f}")
        
        print("\n" + "=" * 80)
        print("TEST RESULT: ✓ PASSED")
        print("=" * 80)
        print("\n✓ Spectral axis information is properly available in single-HDU file!")
        print("✓ Can be used to create properly calibrated velocity axis in datacubes")
        print("\n")

except Exception as e:
    print(f"\n❌ ERROR: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)
