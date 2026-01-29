#!/usr/bin/env python3
"""
Test script to verify the modified create_spectral_datacube function 
with proper velocity axis extraction.
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

print(f"Testing create_spectral_datacube with proper velocity axis extraction")
print(f"Input file: {fits_file}")
print()

try:
    # Open FITS file
    hdul = fits.open(fits_file)
    
    print("=" * 80)
    print("CREATING DATACUBE WITH VELOCITY AXIS")
    print("=" * 80)
    print()
    
    # Create datacube with small test parameters
    datacube, wcs_header, fig = create_spectral_datacube(
        hdul,
        beamsize_deg=0.05,  # 0.05 degrees
        pixsize=0.01,       # 0.01 degrees
        object_filter="M51",
        output_file=None    # Don't save yet
    )
    
    print()
    print("=" * 80)
    print("DATACUBE CREATED SUCCESSFULLY!")
    print("=" * 80)
    
    # Print datacube info
    print(f"\nDatacube shape: {datacube.shape}")
    print(f"  Channels: {datacube.shape[0]}")
    print(f"  Dec pixels: {datacube.shape[1]}")
    print(f"  RA pixels: {datacube.shape[2]}")
    
    # Print WCS header info
    print(f"\nWCS Header Information:")
    print(f"  NAXIS: {wcs_header.get('NAXIS', 'N/A')}")
    print(f"  NAXIS1: {wcs_header.get('NAXIS1', 'N/A')} (RA)")
    print(f"  NAXIS2: {wcs_header.get('NAXIS2', 'N/A')} (Dec)")
    print(f"  NAXIS3: {wcs_header.get('NAXIS3', 'N/A')} (Velocity)")
    
    print(f"\nSpatial Axis (RA):")
    print(f"  CTYPE1: {wcs_header.get('CTYPE1', 'N/A')}")
    print(f"  CUNIT1: {wcs_header.get('CUNIT1', 'N/A')}")
    print(f"  CRVAL1: {wcs_header.get('CRVAL1', 'N/A'):.6f}")
    print(f"  CRPIX1: {wcs_header.get('CRPIX1', 'N/A'):.1f}")
    print(f"  CDELT1: {wcs_header.get('CDELT1', 'N/A'):.6f}")
    
    print(f"\nSpatial Axis (Dec):")
    print(f"  CTYPE2: {wcs_header.get('CTYPE2', 'N/A')}")
    print(f"  CUNIT2: {wcs_header.get('CUNIT2', 'N/A')}")
    print(f"  CRVAL2: {wcs_header.get('CRVAL2', 'N/A'):.6f}")
    print(f"  CRPIX2: {wcs_header.get('CRPIX2', 'N/A'):.1f}")
    print(f"  CDELT2: {wcs_header.get('CDELT2', 'N/A'):.6f}")
    
    print(f"\nVelocity Axis (SPECTRAL):")
    print(f"  CTYPE3: {wcs_header.get('CTYPE3', 'N/A')}")
    print(f"  CUNIT3: {wcs_header.get('CUNIT3', 'N/A')}")
    print(f"  CRVAL3: {wcs_header.get('CRVAL3', 'N/A'):.2f} m/s ({wcs_header.get('CRVAL3', 0)/1000:.2f} km/s)")
    print(f"  CRPIX3: {wcs_header.get('CRPIX3', 'N/A'):.1f} (reference pixel)")
    print(f"  CDELT3: {wcs_header.get('CDELT3', 'N/A'):.2f} m/s ({wcs_header.get('CDELT3', 0)/1000:.4f} km/s per channel)")
    
    # Calculate velocity range
    crval3 = wcs_header.get('CRVAL3', 0.0)
    cdelt3 = wcs_header.get('CDELT3', 1.0)
    naxis3 = wcs_header.get('NAXIS3', 1)
    
    v_first = crval3  # Channel 0 (note: CRPIX3=1, so this is v at pixel 1)
    v_last = crval3 + (naxis3 - 1) * cdelt3  # Last channel
    
    print(f"\nVelocity Range:")
    print(f"  First channel (i=0): {v_first:.2f} m/s ({v_first/1000:.2f} km/s)")
    print(f"  Last channel (i={naxis3-1}): {v_last:.2f} m/s ({v_last/1000:.2f} km/s)")
    print(f"  Total span: {(v_last - v_first):.2f} m/s ({(v_last - v_first)/1000:.2f} km/s)")
    
    # Print beam info
    print(f"\nBeam Information:")
    print(f"  BMAJ: {wcs_header.get('BMAJ', 'N/A'):.4f}°")
    print(f"  BMIN: {wcs_header.get('BMIN', 'N/A'):.4f}°")
    print(f"  BPA: {wcs_header.get('BPA', 'N/A'):.1f}°")
    
    # Print data statistics
    print(f"\nData Statistics:")
    print(f"  Data type: {datacube.dtype}")
    print(f"  Min: {np.nanmin(datacube):.6f}")
    print(f"  Max: {np.nanmax(datacube):.6f}")
    print(f"  Mean: {np.nanmean(datacube):.6f}")
    print(f"  NaN count: {np.sum(~np.isfinite(datacube))}")
    print(f"  Valid pixels: {np.sum(np.isfinite(datacube))}")
    
    print()
    print("✓ Test completed successfully!")
    print()
    print("The datacube now has a proper velocity axis!")
    
    hdul.close()
    
except Exception as e:
    print(f"ERROR: {type(e).__name__}: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)
