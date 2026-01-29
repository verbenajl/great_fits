#!/usr/bin/env python3
"""
Detailed check of PRIMARY HDU header for spectral WCS parameters.
Even if NAXIS=0, the header might contain WCS keywords.
"""

from pathlib import Path
from astropy.io import fits

fits_file = Path("/home/diskB/sofiaobsdata/m51/m51_central_region/kalibrated_all/m51_central_tile_singlehdu.fits")

if not fits_file.exists():
    print(f"Error: FITS file not found: {fits_file}")
    exit(1)

with fits.open(fits_file) as hdul:
    primary_header = hdul[0].header
    
    print("=" * 80)
    print("PRIMARY HDU (HDU 0) - FULL HEADER")
    print("=" * 80)
    print(f"Total cards: {len(primary_header)}\n")
    
    # Print all header cards
    for card in primary_header:
        value = primary_header[card]
        comment = primary_header.comments[card] if card in primary_header.comments else ''
        print(f"{card:8s} = {str(value):40s} / {comment}")
    
    print("\n" + "=" * 80)
    print("SEARCH FOR SPECTRAL WCS KEYWORDS IN PRIMARY HDU")
    print("=" * 80)
    
    spectral_keywords = [
        'CRPIX1', 'CRVAL1', 'CDELT1', 'CTYPE1', 'CUNIT1',
        'CRPIX', 'CRVAL', 'CDELT', 'CTYPE', 'CUNIT',
        'RESTFREQ', 'RESTFRQ',
        'REFPIX', 'REFCHAN',
        'NCHANS', 'NCHAN'
    ]
    
    found = False
    for key in spectral_keywords:
        if key in primary_header:
            value = primary_header[key]
            comment = primary_header.comments[key] if key in primary_header.comments else ''
            print(f"\n{key:15s} = {value:40s} / {comment}")
            found = True
    
    if not found:
        print("No spectral WCS keywords found in primary header")
    
    # Also search by substring
    print("\n" + "=" * 80)
    print("ALL KEYWORDS CONTAINING 'CR', 'CD', 'FREQ', 'REST'")
    print("=" * 80)
    
    for card in primary_header:
        card_upper = card.upper()
        if any(x in card_upper for x in ['CR', 'CD', 'FREQ', 'REST', 'NCHA', 'REFP']):
            value = primary_header[card]
            comment = primary_header.comments[card] if card in primary_header.comments else ''
            print(f"{card:15s} = {str(value):40s} / {comment}")
