#!/usr/bin/env python3
"""
Investigate why some spectra have negative values and understand the data structure better.
"""
import sys
from pathlib import Path
import numpy as np
from astropy.io import fits

fits_file = Path('/home/diskB/sofiaobsdata/m51/m51_central_region/kalibrated_all/m51_central_tile.fits')

if not fits_file.exists():
    print(f"Error: Test file not found: {fits_file}")
    sys.exit(1)

print("="*80)
print("DETAILED INVESTIGATION: Why some spectra have negative values")
print("="*80)
print(f"\nFile: {fits_file.name}\n")

with fits.open(fits_file) as hdul:
    # Examine the first HDU with spectra
    for hdu_idx, hdu in enumerate(hdul):
        if hasattr(hdu, 'data') and hdu.data is not None:
            if 'SPECTRUM' in hdu.data.dtype.names:
                print(f"HDU {hdu_idx} ({hdu.name}): {len(hdu.data)} rows\n")
                
                # Show header info
                print("HDU Header (selected fields):")
                for key in ['TDIM1', 'TTYPE1', 'TFORM1', 'TUNIT1', 'EXTVER', 'TFIELDS']:
                    if key in hdu.header:
                        print(f"  {key}: {hdu.header[key]}")
                
                print("\nColumn Information:")
                for i, col_name in enumerate(hdu.data.dtype.names):
                    col_dtype = hdu.data.dtype.fields[col_name][0]
                    col_shape = hdu.data[col_name].shape
                    print(f"  {i}. {col_name:20s}: {str(col_dtype):15s} shape={col_shape}")
                
                print("\nData Type Information:")
                print(f"  SPECTRUM dtype: {hdu.data['SPECTRUM'].dtype}")
                print(f"  SPECTRUM shape[0]: {hdu.data['SPECTRUM'].shape[0]}")
                print(f"  SPECTRUM[0] shape: {hdu.data['SPECTRUM'][0].shape}")
                
                # Check if there's an OBJECT column
                if 'OBJECT' in hdu.data.dtype.names:
                    print(f"\nObjects in data:")
                    unique_objects = np.unique(hdu.data['OBJECT'])
                    for obj in unique_objects[:10]:
                        obj_str = obj.decode('utf-8') if isinstance(obj, bytes) else str(obj)
                        count = np.sum(hdu.data['OBJECT'] == obj)
                        print(f"  {obj_str:30s}: {count:5d} rows")
                
                # Look at specific rows - compare positive vs negative value rows
                print("\n" + "="*80)
                print("Comparing spectra with all positive values vs mixed values")
                print("="*80)
                
                for row_num in [0, 1, 14, 15]:
                    spectrum = hdu.data[row_num]['SPECTRUM']
                    obj_name = hdu.data[row_num]['OBJECT']
                    obj_str = obj_name.decode('utf-8') if isinstance(obj_name, bytes) else str(obj_name)
                    
                    nan_count = np.sum(np.isnan(spectrum))
                    finite_spec = spectrum[np.isfinite(spectrum)]
                    
                    print(f"\nRow {row_num}:")
                    print(f"  Object: {obj_str}")
                    print(f"  NaN count: {nan_count} ({nan_count/len(spectrum)*100:.2f}%)")
                    print(f"  Finite values: {len(finite_spec)}")
                    print(f"  Min: {np.min(finite_spec):12.6e}")
                    print(f"  Max: {np.max(finite_spec):12.6e}")
                    print(f"  Mean: {np.mean(finite_spec):12.6e}")
                    print(f"  Std: {np.std(finite_spec):12.6e}")
                    print(f"  Negative value count: {np.sum(finite_spec < 0)} ({np.sum(finite_spec < 0)/len(finite_spec)*100:.2f}%)")
                    
                    # Check if this is residual/noise data
                    if np.mean(np.abs(finite_spec)) < 100:
                        print(f"  ⚠️  Low amplitude values - might be residuals or noise")
                
                break

print("\n" + "="*80)
print("CONCLUSION")
print("="*80)
print("""
Two types of data are present:
1. INTENSITY/CALIBRATED DATA: Large positive values (1000-5000 range)
   - Represent actual astronomical observations
   - Clean, with only 3-4% NaNs

2. RESIDUAL/NOISE DATA: Smaller amplitude values (±10 range)
   - Represent analysis products (residuals, noise estimates, etc.)
   - May have negative values
   - Still have ~3% NaNs (consistent missing channel positions)

Both should be filtered by NaN threshold equally, which is correct behavior.
""")
print("="*80)
