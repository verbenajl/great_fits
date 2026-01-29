#!/usr/bin/env python3
"""
Check what columns are in the reduced FITS file and understand the impact of removing VELOCITY/DELTAV.
"""
import sys
from pathlib import Path
from astropy.io import fits

fits_file = Path('/home/diskB/sofiaobsdata/m51/m51_central_fits/fits_files_processed/reduced_data.fits')

if not fits_file.exists():
    print(f"Error: File not found: {fits_file}")
    sys.exit(1)

print("="*80)
print("FITS FILE STRUCTURE ANALYSIS")
print("="*80)
print(f"\nFile: {fits_file.name}\n")

with fits.open(fits_file) as hdul:
    hdul.info()
    
    print("\n" + "="*80)
    print("Column Details")
    print("="*80)
    
    for hdu_idx, hdu in enumerate(hdul):
        if hasattr(hdu, 'data') and hdu.data is not None:
            if hasattr(hdu.data, 'dtype'):
                print(f"\nHDU {hdu_idx} ({hdu.name}): {len(hdu.data)} rows")
                print("-" * 80)
                
                for col_name in hdu.data.dtype.names:
                    col_dtype = hdu.data.dtype.fields[col_name][0]
                    col_data = hdu.data[col_name]
                    
                    if hasattr(col_data[0], '__len__'):
                        # Array column
                        shape = col_data[0].shape
                        size_mb = col_dtype.itemsize * len(hdu.data) * np.prod(shape) / (1024**2)
                        print(f"  {col_name:20s}: {str(col_dtype):15s} shape={shape}  [{size_mb:8.3f} MB]")
                    else:
                        # Scalar column
                        size_mb = col_dtype.itemsize * len(hdu.data) / (1024**2)
                        print(f"  {col_name:20s}: {str(col_dtype):15s}              [{size_mb:8.3f} MB]")

print("\n" + "="*80)
print("COLUMN PRESENCE CHECK")
print("="*80)

with fits.open(fits_file) as hdul:
    for hdu_idx, hdu in enumerate(hdul):
        if hasattr(hdu, 'data') and hdu.data is not None:
            if hasattr(hdu.data, 'dtype'):
                colnames = hdu.data.dtype.names
                
                print(f"\nHDU {hdu_idx}:")
                print(f"  VELOCITY column present:       {'✓ YES' if 'VELOCITY' in colnames else '✗ NO'}")
                print(f"  DELTAV column present:         {'✓ YES' if 'DELTAV' in colnames else '✗ NO'}")
                print(f"  CRPIX1 column present:         {'✓ YES' if 'CRPIX1' in colnames else '✗ NO'}")
                print(f"  VELOCITY_AXIS column present:  {'✓ YES' if 'VELOCITY_AXIS' in colnames else '✗ NO'}")
                
                if 'VELOCITY_AXIS' in colnames:
                    vel_axis = hdu.data['VELOCITY_AXIS']
                    print(f"\n  VELOCITY_AXIS shape: {vel_axis[0].shape}")
                    print(f"  First few values (Row 0): {vel_axis[0][:5]}")

import numpy as np
