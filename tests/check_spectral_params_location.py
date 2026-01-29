#!/usr/bin/env python3
"""
Check where spectral parameters (VELOCITY, DELTAV, etc.) are stored
in the multi-HDU FITS file.
"""

from pathlib import Path
from astropy.io import fits
from astropy.table import Table

fits_file = Path("/home/diskB/sofiaobsdata/m51/m51_central_region/kalibrated_all/m51_central_tile.fits")

if fits_file.exists():
    print(f"Examining: {fits_file}\n")
    
    with fits.open(fits_file) as hdul:
        print(f"{'='*70}")
        print("HDU STRUCTURE")
        print(f"{'='*70}\n")
        hdul.info()
        
        # Check first spectrum HDU
        spectrum_hdu = hdul[1]
        print(f"\n{'='*70}")
        print("FIRST SPECTRUM HDU (HDU 1) HEADER KEYWORDS")
        print(f"{'='*70}\n")
        
        spectral_keywords = ['VELOCITY', 'DELTAV', 'RESTFREQ', 'VELDEF', 'CRPIX1', 'CDELT3']
        for key in spectral_keywords:
            if key in spectrum_hdu.header:
                print(f"✓ {key:<12} = {spectrum_hdu.header[key]}")
            else:
                print(f"✗ {key:<12} NOT FOUND in header")
        
        # Check table columns
        if hasattr(spectrum_hdu, 'columns'):
            print(f"\n{'='*70}")
            print("TABLE COLUMNS IN FIRST SPECTRUM HDU")
            print(f"{'='*70}\n")
            
            print("Available columns:")
            for i, col_name in enumerate(spectrum_hdu.columns.names):
                col = spectrum_hdu.columns[i]
                print(f"  {i:2d}. {col_name:<20} - Format: {col.format}, Shape: {col.shape if hasattr(col, 'shape') else 'scalar'}")
            
            # Look for spectral parameters in columns
            target_cols = ['VELOCITY', 'DELTAV', 'RESTFREQ', 'VELDEF']
            print(f"\nLooking for spectral parameter columns:")
            for col_name in target_cols:
                if col_name in spectrum_hdu.columns.names:
                    col_idx = spectrum_hdu.columns.names.index(col_name)
                    col = spectrum_hdu.columns[col_idx]
                    # Get first few values
                    col_data = spectrum_hdu.data[col_name]
                    if isinstance(col_data, str):
                        print(f"  ✓ {col_name:<12} = {col_data}")
                    else:
                        print(f"  ✓ {col_name:<12} = {col_data[0]} (first row)")
                else:
                    print(f"  ✗ {col_name:<12} NOT in columns")
        
        # Also check if it's in a different HDU
        print(f"\n{'='*70}")
        print("CHECKING ALL HDUs FOR SPECTRAL PARAMETERS")
        print(f"{'='*70}\n")
        
        for hdu_idx, hdu in enumerate(hdul):
            print(f"HDU {hdu_idx}: {hdu.name}")
            
            # Check header
            found_in_header = []
            for key in spectral_keywords:
                if key in hdu.header:
                    found_in_header.append(key)
            
            if found_in_header:
                print(f"  Header keywords: {', '.join(found_in_header)}")
            
            # Check columns if it's a binary table
            if hasattr(hdu, 'columns'):
                found_in_cols = []
                for col_name in ['VELOCITY', 'DELTAV', 'RESTFREQ', 'VELDEF']:
                    if col_name in hdu.columns.names:
                        found_in_cols.append(col_name)
                
                if found_in_cols:
                    print(f"  Table columns: {', '.join(found_in_cols)}")
            print()
else:
    print(f"File not found: {fits_file}")
