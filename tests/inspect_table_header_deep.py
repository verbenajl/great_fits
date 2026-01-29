#!/usr/bin/env python3
"""
Deep inspection of binary table header for reference channel information.
Check all keywords systematically.
"""

from pathlib import Path
from astropy.io import fits

fits_file = Path("/home/diskB/sofiaobsdata/m51/m51_central_region/kalibrated_all/m51_central_tile_singlehdu.fits")

if not fits_file.exists():
    print(f"Error: FITS file not found: {fits_file}")
    exit(1)

with fits.open(fits_file) as hdul:
    table_header = hdul[1].header
    
    print("=" * 80)
    print("BINARY TABLE HDU (HDU 1) - FULL HEADER")
    print("=" * 80)
    print(f"Total cards: {len(table_header)}\n")
    
    # Print all header cards
    for i, card in enumerate(table_header):
        value = table_header[card]
        comment = table_header.comments[card] if card in table_header.comments else ''
        print(f"{i:3d}: {card:12s} = {str(value):40s} / {comment}")
    
    print("\n" + "=" * 80)
    print("KEYWORDS RELATED TO SPECTRAL AXIS")
    print("=" * 80)
    
    # Look specifically for anything that might encode channel info
    for card in table_header:
        card_upper = card.upper()
        value = table_header[card]
        comment = table_header.comments[card] if card in table_header.comments else ''
        
        # Print cards that look relevant
        if any(x in card_upper for x in ['CRPIX', 'CRVAL', 'CDELT', 'CTYPE', 'CUNIT',
                                           'NAXIS', 'NCHA', 'NCHAN', 'FREQ', 'REST',
                                           'VELO', 'VRAD', 'REFP', 'REFC', 'CHAN', 'SPEC']):
            print(f"{card:15s} = {str(value):50s}")
    
    print("\n" + "=" * 80)
    print("TABLE DATA - FIRST OBSERVATION")
    print("=" * 80)
    
    table = hdul[1].data
    
    # Print first row of all columns
    first_row = table[0]
    for col_name in table.names:
        col_value = first_row[col_name]
        if hasattr(col_value, '__len__') and len(col_value) > 10:
            print(f"{col_name:20s} = [array of {len(col_value)} elements]")
        else:
            print(f"{col_name:20s} = {col_value}")
