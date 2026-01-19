#!/usr/bin/env python3
"""
Inspect FITS file columns and data types.
"""
import sys
from astropy.io import fits
from pathlib import Path

def inspect_fits(filepath):
    """Inspect a FITS file and show detailed information."""
    filepath = Path(filepath)
    if not filepath.exists():
        print(f"File not found: {filepath}")
        return
    
    print(f"\n{'='*70}")
    print(f"File: {filepath.name}")
    print(f"Size: {filepath.stat().st_size / 1024 / 1024:.2f} MB")
    print(f"{'='*70}\n")
    
    with fits.open(filepath) as hdul:
        # Show all HDUs
        for i, hdu in enumerate(hdul):
            print(f"HDU {i}: {hdu.name if hasattr(hdu, 'name') else 'Primary'} ({type(hdu).__name__})")
            
            if hasattr(hdu, 'data') and hdu.data is not None:
                data = hdu.data
                print(f"  Rows: {len(data)}")
                
                if hasattr(data, 'dtype'):
                    print(f"  Columns ({len(data.dtype.names)} total):")
                    for col_name in data.dtype.names:
                        col_dtype = data.dtype.fields[col_name][0]
                        col_shape = data[col_name].shape
                        
                        # Calculate size
                        if col_shape[0] > 1:
                            size_bytes = col_dtype.itemsize * len(data) * col_shape[1] if len(col_shape) > 1 else col_dtype.itemsize * len(data) * col_shape[0]
                        else:
                            size_bytes = col_dtype.itemsize * len(data)
                        
                        size_mb = size_bytes / 1024 / 1024
                        
                        if col_name == 'SPECTRUM':
                            # Show details for SPECTRUM column
                            spec_shape = data[col_name][0].shape if hasattr(data[col_name][0], 'shape') else ()
                            print(f"    - {col_name:20s}: {str(col_dtype):15s} shape={spec_shape}  [{size_mb:6.2f} MB]")
                        else:
                            print(f"    - {col_name:20s}: {str(col_dtype):15s}              [{size_mb:6.2f} MB]")
            print()

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python inspect_fits.py <file1> [file2] ...")
        print("\nExample:")
        print("  python inspect_fits.py clean_data.fits reduced_data.fits")
        sys.exit(1)
    
    for filepath in sys.argv[1:]:
        inspect_fits(filepath)
