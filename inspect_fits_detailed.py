#!/usr/bin/env python3
"""
Compare FITS file structure between working and combined files.
"""
import sys
from astropy.io import fits
from pathlib import Path

def inspect_fits_detailed(filepath):
    """Detailed inspection of FITS file."""
    filepath = Path(filepath)
    if not filepath.exists():
        print(f"File not found: {filepath}")
        return
    
    print(f"\n{'='*80}")
    print(f"File: {filepath.name}")
    print(f"{'='*80}\n")
    
    with fits.open(filepath) as hdul:
        # Show primary HDU info
        print(f"Primary HDU:")
        print(f"  Type: {type(hdul[0]).__name__}")
        if hdul[0].header:
            print(f"  Header keywords (first 10):")
            for i, key in enumerate(list(hdul[0].header.keys())[:10]):
                print(f"    {key}: {hdul[0].header[key]}")
        print()
        
        # Show each extension HDU
        for idx, hdu in enumerate(hdul[1:], 1):
            print(f"Extension HDU {idx} ({hdu.name}):")
            print(f"  Type: {type(hdu).__name__}")
            
            if hasattr(hdu, 'data') and hdu.data is not None:
                print(f"  Rows: {len(hdu.data)}")
                print(f"  Header CRVAL2: {hdu.header.get('CRVAL2', 'MISSING')}")
                print(f"  Header CRVAL3: {hdu.header.get('CRVAL3', 'MISSING')}")
                
                # Check CDELT columns
                if 'CDELT2' in hdu.data.dtype.names:
                    cdelt2_vals = hdu.data['CDELT2']
                    print(f"  CDELT2 in data: min={cdelt2_vals.min()}, max={cdelt2_vals.max()}, unique values={len(set(cdelt2_vals))}")
                
                if 'CDELT3' in hdu.data.dtype.names:
                    cdelt3_vals = hdu.data['CDELT3']
                    print(f"  CDELT3 in data: min={cdelt3_vals.min()}, max={cdelt3_vals.max()}, unique values={len(set(cdelt3_vals))}")
            print()

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python inspect_fits_detailed.py <file1> [file2] ...")
        sys.exit(1)
    
    for filepath in sys.argv[1:]:
        inspect_fits_detailed(filepath)
