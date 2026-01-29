#!/usr/bin/env python
"""Check mission date from FITS header."""

from astropy.io import fits
from pathlib import Path

fits_file = Path('/home/diskB/sofiaobsdata/m51/m51_central_fits/fits_files_processed/reduced_data.fits')

with fits.open(fits_file) as hdul:
    header = hdul[0].header
    
    print("Available header keywords related to date/mission:")
    for key in header:
        if any(x in key.lower() for x in ['date', 'miss', 'obs', 'utc', 'time']):
            print(f"  {key}: {header[key]}")
    
    # Check data header too
    print("\nData HDU header:")
    data_header = hdul[1].header
    for key in data_header:
        if any(x in key.lower() for x in ['date', 'miss', 'obs', 'utc', 'time']):
            print(f"  {key}: {data_header[key]}")
