#!/usr/bin/env python3
"""
Detailed inspection of FITS header parameters, especially spectral axis information.
"""

from pathlib import Path
from astropy.io import fits
import sys

# Path to the input FITS file
fits_file = Path("/home/diskB/sofiaobsdata/m51/m51_central_region/kalibrated_all/m51_central_tile_singlehdu.fits")

if not fits_file.exists():
    print(f"Error: FITS file not found: {fits_file}")
    sys.exit(1)

print(f"Inspecting FITS file: {fits_file}")
print(f"File size: {fits_file.stat().st_size / (1024**2):.1f} MB\n")

# Open the FITS file
hdul = fits.open(fits_file)

print("=" * 80)
print("FITS FILE STRUCTURE")
print("=" * 80)
hdul.info()
print()

# Inspect Primary HDU header
print("=" * 80)
print("PRIMARY HDU (HDU 0) HEADER")
print("=" * 80)
if len(hdul) > 0 and hdul[0].header:
    header = hdul[0].header
    print(f"Total header cards: {len(header)}\n")
    
    # Print all header cards
    for card in header:
        comment = header.comments[card] if card in header.comments else ''
        print(f"{card:8s} = {str(header[card]):40s} / {comment}")
    print()

# Inspect Binary Table HDU header (if exists)
if len(hdul) > 1:
    print("=" * 80)
    print("BINARY TABLE HDU (HDU 1) HEADER")
    print("=" * 80)
    header = hdul[1].header
    print(f"Total header cards: {len(header)}\n")
    
    # Print all header cards
    for card in header:
        comment = header.comments[card] if card in header.comments else ''
        print(f"{card:8s} = {str(header[card]):40s} / {comment}")
    print()

# Now look for specific spectral axis parameters
print("=" * 80)
print("SPECTRAL AXIS PARAMETERS")
print("=" * 80)

# Try to find spectral information in either header
header_primary = hdul[0].header if len(hdul) > 0 else None
header_table = hdul[1].header if len(hdul) > 1 else None

print("\nSearching for spectral axis keywords:")
spectral_keywords = [
    # Standard WCS keywords
    'CRVAL1', 'CRPIX1', 'CDELT1', 'CTYPE1', 'CUNIT1',
    # Alternative names
    'CRVAL', 'CRPIX', 'CDELT', 'CTYPE', 'CUNIT',
    # Velocity/frequency related
    'VELO-LSR', 'RESTFREQ', 'RESTFRQ', 
    'FREQ-OBS', 'FELO-LSR', 'SPECSYS',
    # SOFIA specific
    'NAXIS1', 'NAXIS2', 'NAXIS3',
    'EQUINOX', 'RADESYS',
    'NCHANS', 'NLINES',
]

print("\nFrom PRIMARY HDU:")
if header_primary:
    for key in spectral_keywords:
        if key in header_primary:
            comment = header_primary.comments[key] if key in header_primary.comments else ''
            print(f"  {key:15s} = {str(header_primary[key]):40s} / {comment}")

print("\nFrom BINARY TABLE HDU:")
if header_table:
    for key in spectral_keywords:
        if key in header_table:
            comment = header_table.comments[key] if key in header_table.comments else ''
            print(f"  {key:15s} = {str(header_table[key]):40s} / {comment}")

# Check table structure
print("\n" + "=" * 80)
print("BINARY TABLE STRUCTURE")
print("=" * 80)
if len(hdul) > 1 and hdul[1].data is not None:
    table = hdul[1].data
    print(f"\nTable shape: {table.shape}")
    print(f"Columns: {table.names}")
    print()
    
    # Check SPECTRUM column
    if 'SPECTRUM' in table.names:
        spectrum_shape = table['SPECTRUM'].shape
        print(f"SPECTRUM column shape: {spectrum_shape}")
        print(f"  - Number of observations: {spectrum_shape[0]}")
        print(f"  - Number of spectral channels: {spectrum_shape[1]}")
        print()
    
    # Check for any column with 'VELO' or 'FREQ' in name
    print("Columns with spectral information:")
    for col in table.names:
        if any(x in col.upper() for x in ['VELO', 'FREQ', 'CHANNEL', 'WAVE']):
            col_data = table[col]
            if hasattr(col_data, 'shape'):
                print(f"  {col:20s}: shape={col_data.shape}, dtype={col_data.dtype}")
            else:
                print(f"  {col:20s}: scalar")

hdul.close()
