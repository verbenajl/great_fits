#!/usr/bin/env python3
"""
Test writing reduced FITS file with VELOCITY_AXIS column.
"""

import numpy as np
from astropy.io import fits
import tempfile
import sys
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent / 'src'))

from oi_zeigt.reduction.core import reduce_spectra


def create_test_fits():
    """Create a simple FITS file with spectral data for testing."""
    n_spectra = 50
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


def test_write_and_read():
    """Test writing and reading FITS file with VELOCITY_AXIS."""
    print("=" * 80)
    print("TEST: Write and Read FITS with VELOCITY_AXIS")
    print("=" * 80)
    
    # Create test FITS
    print("\n1. Creating test FITS file...")
    hdul = create_test_fits()
    print("   ✓ Created test FITS with 50 spectra, 1264 channels")
    
    # Apply reduction
    print("\n2. Applying reduce_spectra with baseline correction...")
    methods = {
        'baseline': {'order': 2, 'window': (100, 120)}
    }
    
    try:
        reduced_hdul = reduce_spectra(hdul, methods=methods)
        print("   ✓ Reduction completed")
    except Exception as e:
        print(f"   ✗ Error during reduction: {e}")
        return False
    
    # Write to temporary file
    print("\n3. Writing reduced FITS to temporary file...")
    with tempfile.NamedTemporaryFile(suffix='.fits', delete=False) as tmp:
        tmp_path = tmp.name
    
    try:
        reduced_hdul.writeto(tmp_path, overwrite=True)
        print(f"   ✓ Written to: {tmp_path}")
        print(f"   ✓ File size: {Path(tmp_path).stat().st_size / 1024 / 1024:.1f} MB")
    except Exception as e:
        print(f"   ✗ Error writing file: {e}")
        return False
    
    # Read back and verify
    print("\n4. Reading back FITS file...")
    try:
        with fits.open(tmp_path) as hdul_read:
            print(f"   ✓ Opened successfully")
            print(f"   ✓ Number of HDUs: {len(hdul_read)}")
            
            table_data = hdul_read[1].data
            print(f"   ✓ Table columns: {list(table_data.dtype.names)}")
            
            # Check VELOCITY_AXIS
            if 'VELOCITY_AXIS' not in table_data.dtype.names:
                print("   ✗ VELOCITY_AXIS column not found!")
                return False
            
            print("   ✓ VELOCITY_AXIS column found")
            
            # Check data
            vel_axis = table_data['VELOCITY_AXIS'][0]
            spectrum = table_data['SPECTRUM'][0]
            
            print(f"\n5. Data verification:")
            print(f"   SPECTRUM shape: {spectrum.shape}")
            print(f"   VELOCITY_AXIS shape: {vel_axis.shape}")
            print(f"   SPECTRUM dtype: {spectrum.dtype}")
            print(f"   VELOCITY_AXIS dtype: {vel_axis.dtype}")
            
            print(f"\n6. Velocity axis sample values:")
            print(f"   First channel:  {vel_axis[0]:12.1f} m/s")
            print(f"   Channel 505:    {vel_axis[505]:12.1f} m/s (reference)")
            print(f"   Last channel:   {vel_axis[-1]:12.1f} m/s")
            
            # Verify calculation
            expected_first = 470000.0 + (0 - 505) * 500.0
            expected_ref = 470000.0
            expected_last = 470000.0 + (1263 - 505) * 500.0
            
            print(f"\n7. Verification against expected values:")
            tol = 1e-6
            checks = [
                (vel_axis[0], expected_first, "First channel"),
                (vel_axis[505], expected_ref, "Reference channel"),
                (vel_axis[-1], expected_last, "Last channel"),
            ]
            
            all_correct = True
            for actual, expected, label in checks:
                diff = abs(actual - expected)
                status = "✓" if diff < tol else "✗"
                print(f"   {status} {label:20s}: {actual:12.1f} (expected {expected:12.1f})")
                if diff >= tol:
                    all_correct = False
            
            # Test plotting capability
            print(f"\n8. Test plotting capability:")
            try:
                import matplotlib
                matplotlib.use('Agg')  # Non-interactive backend
                import matplotlib.pyplot as plt
                
                fig, ax = plt.subplots(figsize=(12, 5))
                ax.plot(vel_axis, spectrum, 'b-', linewidth=0.5)
                ax.set_xlabel('Velocity (m/s)')
                ax.set_ylabel('Flux')
                ax.set_title('Reduced Spectrum with Velocity Axis')
                ax.grid(True, alpha=0.3)
                
                plot_file = tmp_path.replace('.fits', '_plot.png')
                fig.savefig(plot_file, dpi=100, bbox_inches='tight')
                plt.close(fig)
                print(f"   ✓ Plot saved to: {plot_file}")
            except ImportError:
                print("   ! Matplotlib not available (skipping plot test)")
            except Exception as e:
                print(f"   ! Could not create plot: {e}")
            
            # Clean up
            Path(tmp_path).unlink()
            print(f"\n✓ Cleaned up temporary file")
            
            return all_correct
            
    except Exception as e:
        print(f"   ✗ Error reading file: {e}")
        import traceback
        traceback.print_exc()
        return False


if __name__ == '__main__':
    success = test_write_and_read()
    
    print("\n" + "=" * 80)
    if success:
        print("RESULT: ✓ ALL TESTS PASSED")
        sys.exit(0)
    else:
        print("RESULT: ✗ TESTS FAILED")
        sys.exit(1)
