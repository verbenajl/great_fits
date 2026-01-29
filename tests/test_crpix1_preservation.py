#!/usr/bin/env python3
"""
Test script to verify CRPIX1 and other spectral parameters are preserved
when combining FITS files into a single HDU.
"""

from pathlib import Path
from astropy.io import fits
import tempfile
import shutil

# Use the multi-HDU file as input
fits_file = Path("/home/diskB/sofiaobsdata/m51/m51_central_region/kalibrated_all/m51_central_tile.fits")

if not fits_file.exists():
    print(f"Error: FITS file not found: {fits_file}")
    print("Trying local test file instead...")
    fits_file = None

if fits_file:
    print(f"\n{'='*70}")
    print("CHECKING ORIGINAL MULTI-HDU FILE")
    print(f"{'='*70}")
    
    with fits.open(fits_file) as hdul:
        print(f"\nFile: {fits_file.name}")
        print(f"Number of HDUs: {len(hdul)}")
        
        # Check first spectrum HDU (should be HDU 1)
        spectrum_hdu = hdul[1]
        print(f"\nSpectrum HDU (HDU 1): {spectrum_hdu.name}")
        print(f"Data shape: {spectrum_hdu.data.shape if hasattr(spectrum_hdu, 'data') else 'N/A'}")
        
        print(f"\nSpectral parameters in first spectrum HDU:")
        spectral_keys = ['CRVAL1', 'CRPIX1', 'CTYPE1', 'CUNIT1', 'CDELT1',
                         'CRVAL2', 'CRPIX2', 'CTYPE2', 'CUNIT2', 'CDELT2',
                         'CRVAL3', 'CRPIX3', 'CTYPE3', 'CUNIT3', 'CDELT3',
                         'VELOCITY', 'DELTAV', 'RESTFREQ', 'VELDEF']
        
        for key in spectral_keys:
            if key in spectrum_hdu.header:
                value = spectrum_hdu.header[key]
                print(f"  {key:<12} = {value}")
        
        # Extract CRPIX1 for comparison
        original_crpix1 = spectrum_hdu.header.get('CRPIX1', None)
        original_velocity = spectrum_hdu.header.get('VELOCITY', None)
        original_deltav = spectrum_hdu.header.get('DELTAV', None)
        original_restfreq = spectrum_hdu.header.get('RESTFREQ', None)

# Now test the combine_fits function
print(f"\n{'='*70}")
print("TESTING SINGLE-HDU COMBINATION")
print(f"{'='*70}\n")

from src.oi_zeigt.basic_io import combine_fits_files

# For testing, we'll create a temporary test list with just the first file
# and combine it using the single_hdu option
if fits_file:
    # Import numpy for creating test data
    import numpy as np
    
    # Create a test with the actual multi-HDU file
    test_fits_list = [str(fits_file)]
    
    print(f"Combining {len(test_fits_list)} FITS file(s) with single_hdu=True...")
    
    try:
        combined_hdul = combine_fits_files(test_fits_list, single_hdu=True)
        
        print(f"✓ Successfully combined FITS files")
        print(f"\nCombined HDUList structure:")
        combined_hdul.info()
        
        # Get the spectrum HDU from combined file (should be HDU 1)
        combined_spectrum_hdu = combined_hdul[1]
        
        print(f"\n{'='*70}")
        print("CHECKING SINGLE-HDU OUTPUT FILE")
        print(f"{'='*70}")
        print(f"\nSpectrum HDU (HDU 1): {combined_spectrum_hdu.name}")
        
        print(f"\nSpectral parameters in combined single-HDU file:")
        print(f"  {'Key':<12} {'Original':<20} {'Combined':<20} {'Status':<10}")
        print(f"  {'-'*12} {'-'*20} {'-'*20} {'-'*10}")
        
        for key in spectral_keys:
            original_val = None
            combined_val = None
            
            if fits_file:
                with fits.open(fits_file) as hdul:
                    if key in hdul[1].header:
                        original_val = hdul[1].header[key]
            
            if key in combined_spectrum_hdu.header:
                combined_val = combined_spectrum_hdu.header[key]
            
            # Check if preserved
            if original_val is not None and combined_val is not None:
                status = "✓ PRESERVED" if original_val == combined_val else "⚠ CHANGED"
            elif original_val is None and combined_val is not None:
                status = "⚠ ADDED"
            elif original_val is not None and combined_val is None:
                status = "✗ MISSING"
            else:
                status = "- N/A"
            
            original_str = str(original_val) if original_val is not None else "N/A"
            combined_str = str(combined_val) if combined_val is not None else "MISSING!"
            
            print(f"  {key:<12} {original_str:<20} {combined_str:<20} {status:<10}")
        
        # Special check for CRPIX1
        print(f"\n{'='*70}")
        print("CRITICAL CHECK: CRPIX1 PRESERVATION")
        print(f"{'='*70}")
        
        if fits_file:
            with fits.open(fits_file) as hdul:
                orig_crpix1 = hdul[1].header.get('CRPIX1', None)
        
        comb_crpix1 = combined_spectrum_hdu.header.get('CRPIX1', None)
        
        if orig_crpix1 is not None:
            print(f"✓ Original CRPIX1 found: {orig_crpix1}")
        else:
            print(f"✗ Original CRPIX1 NOT found")
        
        if comb_crpix1 is not None:
            print(f"✓ Combined CRPIX1 found: {comb_crpix1}")
            if orig_crpix1 is not None and orig_crpix1 == comb_crpix1:
                print(f"✓ CRPIX1 CORRECTLY PRESERVED!")
            elif orig_crpix1 is not None:
                print(f"⚠ WARNING: CRPIX1 value changed!")
        else:
            print(f"✗ Combined CRPIX1 MISSING - FIX FAILED!")
        
        # Test writing to file
        print(f"\n{'='*70}")
        print("TESTING FILE WRITING")
        print(f"{'='*70}\n")
        
        with tempfile.NamedTemporaryFile(suffix='.fits', delete=False) as tmp:
            tmp_path = Path(tmp.name)
        
        print(f"Writing combined FITS to temporary file: {tmp_path}")
        combined_hdul.writeto(tmp_path, overwrite=True)
        print(f"✓ Successfully wrote file")
        
        # Read back and verify
        print(f"\nReading back combined FITS file...")
        with fits.open(tmp_path) as verify_hdul:
            verify_spectrum_hdu = verify_hdul[1]
            verify_crpix1 = verify_spectrum_hdu.header.get('CRPIX1', None)
            verify_velocity = verify_spectrum_hdu.header.get('VELOCITY', None)
            verify_deltav = verify_spectrum_hdu.header.get('DELTAV', None)
            
            print(f"✓ Successfully read file back")
            print(f"\nVerification of written file:")
            print(f"  CRPIX1: {verify_crpix1}")
            print(f"  VELOCITY: {verify_velocity}")
            print(f"  DELTAV: {verify_deltav}")
            
            if verify_crpix1 == comb_crpix1:
                print(f"\n✓ CRPIX1 successfully persisted to file!")
            else:
                print(f"\n⚠ WARNING: CRPIX1 changed when writing to file!")
        
        # Clean up
        tmp_path.unlink()
        print(f"✓ Cleaned up temporary file")
        
    except Exception as e:
        print(f"✗ Error during combination: {e}")
        import traceback
        traceback.print_exc()
else:
    print("✗ Cannot run test: FITS file not found")

print(f"\n{'='*70}")
print("TEST COMPLETE")
print(f"{'='*70}\n")
