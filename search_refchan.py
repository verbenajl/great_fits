#!/usr/bin/env python3
"""
Search for reference channel parameters in FITS headers.
Look for REFCHAN and other similar keywords.
"""

from pathlib import Path
from astropy.io import fits
import numpy as np

fits_file = Path("/home/diskB/sofiaobsdata/m51/m51_central_region/kalibrated_all/m51_central_tile_singlehdu.fits")

if not fits_file.exists():
    print(f"Error: FITS file not found: {fits_file}")
    exit(1)

print("Searching for reference channel parameters in FITS headers")
print("=" * 80)

with fits.open(fits_file) as hdul:
    # Check primary header
    print("\nPRIMARY HDU HEADER (HDU 0):")
    print("-" * 80)
    primary_header = hdul[0].header
    
    # Print all keywords
    for card in primary_header:
        print(f"  {card:15s} = {str(primary_header[card]):40s}")
    
    # Check binary table header
    print("\nBINARY TABLE HDU HEADER (HDU 1):")
    print("-" * 80)
    table_header = hdul[1].header
    
    # Search for specific keywords
    keywords_to_search = [
        'REFCHAN', 'REFPIX', 'CRPIX', 'CRVAL', 'CDELT', 'CTYPE', 'CUNIT',
        'CRPIX1', 'CRVAL1', 'CDELT1', 'CTYPE1', 'CUNIT1',
        'NAXIS1', 'NAXIS2', 'NAXIS3',
        'RESTFREQ', 'RESTFRQ', 'FREQ', 'FELO', 'VELO',
        'VSYS', 'VLSR', 'VREF', 'VRADIO',
        'RFREQ', 'OBSFREQ',
    ]
    
    print("\nSearching for reference/spectral keywords:")
    print("-" * 80)
    found_keys = []
    for key in keywords_to_search:
        if key in table_header:
            value = table_header[key]
            comment = table_header.comments[key] if key in table_header.comments else ''
            print(f"  {key:15s} = {str(value):40s} / {comment}")
            found_keys.append(key)
    
    if not found_keys:
        print("  No matching keywords found")
    
    # Also search for any keyword containing 'REF' or 'CHAN'
    print("\nAll keywords containing 'REF' or 'CHAN':")
    print("-" * 80)
    for card in table_header:
        if 'REF' in card.upper() or 'CHAN' in card.upper():
            value = table_header[card]
            comment = table_header.comments[card] if card in table_header.comments else ''
            print(f"  {card:15s} = {str(value):40s} / {comment}")
    
    # Check table columns for reference information
    print("\nTable Columns (first 5):")
    print("-" * 80)
    table = hdul[1].data
    for i, col in enumerate(table.names[:5]):
        print(f"  {col:20s}")
    
    # Check if any columns have 'REF' or 'CHAN' in name
    print("\nTable Columns with 'REF', 'CHAN', or 'FREQ' in name:")
    print("-" * 80)
    for col in table.names:
        if any(x in col.upper() for x in ['REF', 'CHAN', 'FREQ', 'VELO']):
            col_data = table[col]
            if hasattr(col_data, 'shape'):
                shape_str = str(col_data.shape)
                dtype_str = str(col_data.dtype)
                print(f"  {col:20s} shape={shape_str:20s} dtype={dtype_str}")
                # Print first few values
                try:
                    if len(col_data.shape) == 1 and len(col_data) > 0:
                        print(f"    First 3 values: {col_data[:3]}")
                except:
                    pass
            else:
                print(f"  {col:20s} (scalar)")
    
    # Try to find frequency/velocity information more systematically
    print("\n" + "=" * 80)
    print("FREQUENCY/VELOCITY INFORMATION FROM TABLE COLUMNS")
    print("=" * 80)
    
    critical_cols = ['LOFREQ', 'RESTFREQ', 'FOFFSET', 'DELTAV', 'VELOCITY', 'SPECTRUM']
    for col_name in critical_cols:
        if col_name in table.names:
            col_data = table[col_name]
            print(f"\n{col_name}:")
            print(f"  Shape: {col_data.shape}")
            print(f"  Data type: {col_data.dtype}")
            if len(col_data.shape) == 1:
                print(f"  Min: {np.min(col_data)}")
                print(f"  Max: {np.max(col_data)}")
                print(f"  First 5 values: {col_data[:5]}")
            else:
                print(f"  First value shape: {col_data[0].shape if len(col_data) > 0 else 'empty'}")
                if len(col_data) > 0:
                    print(f"  First value: {col_data[0]}")
