#!/usr/bin/env python3
"""
Check CRPIX1 and other WCS parameters in the original FITS file.
"""

from pathlib import Path
from astropy.io import fits

fits_file = Path("/home/diskB/sofiaobsdata/m51/m51_central_region/kalibrated_all/m51_central_tile_singlehdu.fits")

if not fits_file.exists():
    print(f"Error: FITS file not found: {fits_file}")
    exit(1)

with fits.open(fits_file) as hdul:
    # Check both headers
    print("PRIMARY HDU HEADER:")
    print("=" * 80)
    primary_header = hdul[0].header
    for key in ['CRPIX1', 'CRVAL1', 'CDELT1', 'CTYPE1', 'CUNIT1', 
                'CRPIX2', 'CRVAL2', 'CDELT2', 'CTYPE2', 'CUNIT2',
                'CRPIX3', 'CRVAL3', 'CDELT3', 'CTYPE3', 'CUNIT3',
                'NAXIS1', 'NAXIS2', 'NAXIS3', 'NAXIS',
                'RESTFREQ', 'RESTFRQ']:
        if key in primary_header:
            value = primary_header[key]
            print(f"  {key:15s} = {value}")
        else:
            print(f"  {key:15s} = NOT FOUND")
    
    print("\n" + "=" * 80)
    print("BINARY TABLE HDU HEADER:")
    print("=" * 80)
    table_header = hdul[1].header
    for key in ['CRPIX1', 'CRVAL1', 'CDELT1', 'CTYPE1', 'CUNIT1',
                'CRPIX2', 'CRVAL2', 'CDELT2', 'CTYPE2', 'CUNIT2',
                'CRPIX3', 'CRVAL3', 'CDELT3', 'CTYPE3', 'CUNIT3',
                'NAXIS1', 'NAXIS2', 'NAXIS3', 'NAXIS',
                'RESTFREQ', 'RESTFRQ']:
        if key in table_header:
            value = table_header[key]
            print(f"  {key:15s} = {value}")
        else:
            print(f"  {key:15s} = NOT FOUND")
    
    # Also check table columns
    print("\n" + "=" * 80)
    print("TABLE STRUCTURE:")
    print("=" * 80)
    table = hdul[1].data
    print(f"Number of rows: {len(table)}")
    print(f"SPECTRUM column shape: {table['SPECTRUM'].shape}")
    print(f"Number of spectral channels: {table['SPECTRUM'].shape[1]}")
    
    # Check for CRPIX-like columns
    print("\nColumns that might contain reference pixel info:")
    for col in table.names:
        if 'CRPIX' in col.upper() or 'REFPIX' in col.upper() or 'REFCHAN' in col.upper():
            print(f"  {col}: {table[col][:5]}")
