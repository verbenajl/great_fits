#!/usr/bin/env python
"""Check all header keys."""

from astropy.io import fits
from pathlib import Path

fits_file = Path('/home/diskB/sofiaobsdata/m51/m51_central_fits/fits_files_processed/reduced_data.fits')

with fits.open(fits_file) as hdul:
    print("Primary HDU header (first 50 keys):")
    for i, key in enumerate(hdul[0].header):
        if i < 50:
            print(f"  {key}: {hdul[0].header[key]}")
    
    print("\n\nData HDU header columns:")
    print(hdul[1].data.dtype.names[:20])
