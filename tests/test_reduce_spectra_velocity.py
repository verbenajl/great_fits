#!/usr/bin/env python3
"""
Test script to verify VELOCITY_AXIS column is added to reduce_spectra output.
"""

import numpy as np
from astropy.io import fits
from astropy.table import Table
import tempfile
import sys
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent / 'src'))

from oi_zeigt.reduction.core import reduce_spectra


def create_test_fits():
    """Create a simple FITS file with spectral data for testing."""
    # Create test data
    n_spectra = 100
    n_channels = 1264
    
    # Create fake spectra
    spectra = np.random.randn(n_spectra, n_channels).astype(np.float32)
    
    # Create table columns
    columns = [
        fits.Column(name='OBJECT', format='20A', array=np.full(n_spectra, 'TEST', dtype='U20')),
        fits.Column(name='SPECTRUM', format=f'{n_channels}E', array=spectra),
        fits.Column(name='VELOCITY', format='D', array=np.full(n_spectra, 470000.0)),
        fits.Column(name='DELTAV', format='D', array=np.full(n_spectra, 500.0)),
    ]
    
    # Create binary table HDU
    table_hdu = fits.BinTableHDU.from_columns(columns)
    table_hdu.name = 'SPECTRA'
    table_hdu.header['CRPIX1'] = 506
    
    # Create primary HDU
    primary_hdu = fits.PrimaryHDU()
    
    # Create HDUList
    hdul = fits.HDUList([primary_hdu, table_hdu])
    
    return hdul


def test_velocity_axis_addition():
    """Test that VELOCITY_AXIS column is properly added."""
    print("Creating test FITS file...")
    hdul = create_test_fits()
    
    print("Running reduce_spectra with baseline reduction...")
    methods = {
        'baseline': {'order': 2, 'window': (100, 120)}
    }
    
    try:
        reduced_hdul = reduce_spectra(hdul, methods=methods)
        
        print("\n✓ reduce_spectra completed successfully")
        
        # Check if VELOCITY_AXIS was added
        table_data = reduced_hdul[1].data
        
        if 'VELOCITY_AXIS' in table_data.dtype.names:
            print("✓ VELOCITY_AXIS column found in output")
            
            velocity_axis = table_data['VELOCITY_AXIS'][0]
            print(f"  Shape: {velocity_axis.shape}")
            print(f"  First value: {velocity_axis[0]:.1f} m/s")
            print(f"  Middle value (ref): {velocity_axis[505]:.1f} m/s")
            print(f"  Last value: {velocity_axis[-1]:.1f} m/s")
            
            # Verify velocity axis calculation
            # Expected: v[i] = 470000 + (i - 505) * 500
            expected_first = 470000.0 + (0 - 505) * 500.0
            expected_ref = 470000.0 + (505 - 505) * 500.0
            expected_last = 470000.0 + (1263 - 505) * 500.0
            
            print(f"\n  Expected first: {expected_first:.1f} m/s")
            print(f"  Expected ref: {expected_ref:.1f} m/s")
            print(f"  Expected last: {expected_last:.1f} m/s")
            
            tol = 1e-6
            if (abs(velocity_axis[0] - expected_first) < tol and
                abs(velocity_axis[505] - expected_ref) < tol and
                abs(velocity_axis[-1] - expected_last) < tol):
                print("\n✓ Velocity axis values are correct!")
                return True
            else:
                print("\n✗ Velocity axis values don't match expected!")
                return False
        else:
            print("✗ VELOCITY_AXIS column NOT found in output")
            print(f"  Available columns: {table_data.dtype.names}")
            return False
    
    except Exception as e:
        print(f"\n✗ Error during reduce_spectra: {e}")
        import traceback
        traceback.print_exc()
        return False


if __name__ == '__main__':
    success = test_velocity_axis_addition()
    sys.exit(0 if success else 1)
