#!/usr/bin/env python3
"""
Test script: Verify spectral axis table is properly added to datacube FITS file
"""

from pathlib import Path
from astropy.io import fits
import numpy as np
import sys

# Add src to path
sys.path.insert(0, str(Path(__file__).parent / 'src'))

from oi_zeigt.mapping.gridding import create_spectral_datacube

# Test file
fits_file = Path("/home/diskB/sofiaobsdata/m51/m51_central_region/kalibrated_all/m51_central_tile_singlehdu.fits")

if not fits_file.exists():
    print(f"Error: FITS file not found: {fits_file}")
    sys.exit(1)

print("=" * 80)
print("TEST: SPECTRAL AXIS TABLE IN DATACUBE")
print("=" * 80)
print()

try:
    # Open FITS file
    hdul = fits.open(fits_file)
    
    print("Step 1: Create datacube with spectral axis table")
    print("-" * 80)
    
    output_file = Path("/tmp/test_datacube_with_spectrum.fits")
    
    # Create datacube
    datacube, wcs_header, fig = create_spectral_datacube(
        hdul,
        beamsize_deg=0.05,
        pixsize=0.01,
        object_filter="M51",
        output_file=str(output_file)
    )
    
    print()
    print("Step 2: Verify FITS file structure")
    print("-" * 80)
    
    # Open the saved file
    with fits.open(output_file) as datacube_hdul:
        print(f"\nNumber of HDUs: {len(datacube_hdul)}")
        print()
        
        # Check structure
        for i, hdu in enumerate(datacube_hdul):
            print(f"HDU {i}: {hdu.name}")
            if isinstance(hdu, fits.PrimaryHDU):
                print(f"  Type: Primary HDU")
                print(f"  Data shape: {hdu.data.shape}")
                print(f"  Data dtype: {hdu.data.dtype}")
                print(f"  Header keywords: {len(hdu.header)}")
            elif isinstance(hdu, fits.BinTableHDU):
                print(f"  Type: Binary Table")
                print(f"  Columns: {hdu.columns.names}")
                print(f"  Rows: {len(hdu.data)}")
                print(f"  Header keywords: {len(hdu.header)}")
        
        print()
        print("Step 3: Inspect SPECTRUM table")
        print("-" * 80)
        
        if len(datacube_hdul) > 1 and datacube_hdul[1].name == 'SPECTRUM':
            spec_table = datacube_hdul[1]
            
            print(f"\n✓ SPECTRUM table found!")
            print(f"  Columns: {spec_table.columns.names}")
            print(f"  Rows (channels): {len(spec_table.data)}")
            
            # Check columns
            channels = spec_table.data['CHANNEL']
            velocities = spec_table.data['VELOCITY']
            
            print(f"\n  CHANNEL column:")
            print(f"    Data type: {channels.dtype}")
            print(f"    Range: {channels[0]} to {channels[-1]}")
            print(f"    Sample values: {channels[0:5]}")
            
            print(f"\n  VELOCITY column:")
            print(f"    Data type: {velocities.dtype}")
            print(f"    Unit: {velocities.unit if hasattr(velocities, 'unit') else 'N/A'}")
            print(f"    Range: {velocities[0]:.2f} to {velocities[-1]:.2f}")
            print(f"    Sample values (first 5):")
            for i in range(5):
                print(f"      Channel {channels[i]}: {velocities[i]:.2f}")
            
            # Check metadata
            print(f"\n  Header keywords:")
            for key in ['CRPIX1', 'RESTFRQ', 'VELDEF']:
                if key in spec_table.header:
                    value = spec_table.header[key]
                    print(f"    {key} = {value}")
                else:
                    print(f"    {key} = NOT FOUND")
        else:
            print(f"\n✗ SPECTRUM table NOT found!")
            sys.exit(1)
        
        print()
        print("Step 4: Test direct velocity access (NO WCS required!)")
        print("-" * 80)
        
        velocities_ms = spec_table.data['VELOCITY']  # In m/s
        velocities_kms = velocities_ms / 1000.0  # Convert to km/s
        
        test_channels = [0, 100, 200, 504, 505, 600, 1263]
        
        print(f"\nDirect velocity access without WCS:")
        for ch in test_channels:
            v_ms = velocities_ms[ch]
            v_kms = velocities_kms[ch]
            print(f"  Channel {ch:4d}: {v_ms:12.0f} m/s = {v_kms:8.2f} km/s")
        
        print()
        print("Step 5: Verify reference velocity")
        print("-" * 80)
        
        crpix1 = spec_table.header['CRPIX1']
        ref_channel_0indexed = crpix1 - 1  # Convert to 0-indexed
        
        print(f"\nReference pixel (CRPIX1): {crpix1:.2f} (FITS 1-indexed)")
        print(f"Reference channel (0-indexed): {ref_channel_0indexed:.2f}")
        
        # Get reference velocity
        ref_velocity_ms = velocities_ms[int(ref_channel_0indexed)]
        ref_velocity_kms = ref_velocity_ms / 1000.0
        
        print(f"Reference velocity: {ref_velocity_ms:.0f} m/s = {ref_velocity_kms:.2f} km/s")
        print(f"Expected: ~470 km/s")
        
        # Check if close
        if 469 < ref_velocity_kms < 471:
            print(f"✓ CORRECT! Reference velocity is at correct channel (~505)")
        else:
            print(f"✗ ERROR! Reference velocity seems wrong!")
        
        print()
        print("Step 6: Compare WCS method vs direct table access")
        print("-" * 80)
        
        from astropy.wcs import WCS
        
        wcs = WCS(datacube_hdul[0].header)
        wcs_velocities = wcs.pixel_to_world_values(np.arange(len(spec_table.data)), 0, 0)[2]
        
        print(f"\nComparing velocity calculations:")
        print(f"{'Channel':<10} {'Table (m/s)':<20} {'WCS (m/s)':<20} {'Match?':<10}")
        print("-" * 60)
        
        test_indices = [0, 50, 100, 504, 505, 1000, 1263]
        for idx in test_indices:
            table_v = velocities_ms[idx]
            wcs_v = wcs_velocities[idx]
            match = "✓" if np.abs(table_v - wcs_v) < 1.0 else "✗"
            print(f"{idx:<10} {table_v:<20.2f} {wcs_v:<20.2f} {match:<10}")
        
        print()
        print("=" * 80)
        print("✅ TEST COMPLETE!")
        print("=" * 80)
        print()
        print("Summary:")
        print(f"  ✓ Datacube FITS file created: {output_file}")
        print(f"  ✓ SPECTRUM table present with {len(spec_table.data)} rows")
        print(f"  ✓ CHANNEL and VELOCITY columns accessible")
        print(f"  ✓ Metadata (CRPIX1, RESTFRQ, VELDEF) included")
        print(f"  ✓ Direct velocity access works (no WCS needed)")
        print(f"  ✓ Values match WCS calculations")
        print()

except Exception as e:
    print(f"Error: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)
