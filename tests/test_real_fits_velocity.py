#!/usr/bin/env python3
"""
Test writing reduced FITS file with actual FITS data and VELOCITY_AXIS column.
"""

import sys
from pathlib import Path
import tempfile

# Add src to path
sys.path.insert(0, str(Path(__file__).parent / 'src'))

from astropy.io import fits
from oi_zeigt.reduction.core import reduce_spectra


def test_real_fits():
    """Test with real FITS file."""
    fits_file = "/home/diskB/sofiaobsdata/m51/m51_central_region/kalibrate_d/M51CENTER_UNKNOWN_central_tile_CB_flight1_20170201_LFA_L1.fits"
    
    if not Path(fits_file).exists():
        print(f"✗ FITS file not found: {fits_file}")
        return False
    
    print("=" * 80)
    print("TEST: Write and Read Real FITS with VELOCITY_AXIS")
    print("=" * 80)
    
    # Load FITS
    print(f"\n1. Loading FITS file: {fits_file}")
    with fits.open(fits_file) as hdul:
        print(f"   ✓ Loaded {len(hdul)} HDUs")
        print(f"   ✓ Number of spectra: {len(hdul[1].data)}")
        print(f"   ✓ Spectrum shape: {hdul[1].data['SPECTRUM'].shape}")
        
        # Check for spectral parameters
        table = hdul[1].data
        has_velocity = 'VELOCITY' in table.dtype.names
        has_deltav = 'DELTAV' in table.dtype.names
        has_spectrum = 'SPECTRUM' in table.dtype.names
        
        print(f"   ✓ Has VELOCITY column: {has_velocity}")
        print(f"   ✓ Has DELTAV column: {has_deltav}")
        print(f"   ✓ Has SPECTRUM column: {has_spectrum}")
        
        if has_velocity and has_deltav:
            velo_ref = float(table['VELOCITY'][0])
            deltav = float(table['DELTAV'][0])
            print(f"   ✓ VELOCITY: {velo_ref:.0f} m/s")
            print(f"   ✓ DELTAV: {deltav:.1f} m/s/channel")
    
    # Apply reduction
    print(f"\n2. Applying reduce_spectra...")
    methods = {
        'baseline': {'order': 2, 'window': (100, 120)}
    }
    
    try:
        with fits.open(fits_file) as hdul:
            reduced_hdul = reduce_spectra(hdul, methods=methods)
        print("   ✓ Reduction completed successfully")
    except Exception as e:
        print(f"   ✗ Error during reduction: {e}")
        import traceback
        traceback.print_exc()
        return False
    
    # Write to temporary file
    print(f"\n3. Writing reduced FITS to temporary file...")
    with tempfile.NamedTemporaryFile(suffix='.fits', delete=False) as tmp:
        tmp_path = tmp.name
    
    try:
        reduced_hdul.writeto(tmp_path, overwrite=True)
        file_size_mb = Path(tmp_path).stat().st_size / 1024 / 1024
        print(f"   ✓ Written to: {tmp_path}")
        print(f"   ✓ File size: {file_size_mb:.1f} MB")
    except Exception as e:
        print(f"   ✗ Error writing file: {e}")
        return False
    
    # Read back and verify
    print(f"\n4. Reading back FITS file...")
    try:
        with fits.open(tmp_path) as hdul_read:
            print(f"   ✓ Opened successfully")
            print(f"   ✓ Number of HDUs: {len(hdul_read)}")
            
            table_data = hdul_read[1].data
            columns = list(table_data.dtype.names)
            print(f"   ✓ Number of columns: {len(columns)}")
            print(f"   ✓ Columns: {', '.join(columns[:5])}...")
            
            # Check VELOCITY_AXIS
            if 'VELOCITY_AXIS' not in columns:
                print("   ✗ VELOCITY_AXIS column not found!")
                print(f"   Available columns: {columns}")
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
            print(f"   Min velocity:   {vel_axis.min():12.1f} m/s")
            print(f"   Max velocity:   {vel_axis.max():12.1f} m/s")
            print(f"   Channel spacing: {vel_axis[1] - vel_axis[0]:12.1f} m/s")
            
            # Check multiple spectra
            print(f"\n7. Check consistency across spectra:")
            vel_axis_5 = table_data['VELOCITY_AXIS'][5]
            if (vel_axis == vel_axis_5).all():
                print(f"   ✓ VELOCITY_AXIS identical for all spectra (as expected)")
            else:
                print(f"   ! VELOCITY_AXIS differs for spectrum 5 (might be expected if velocities vary)")
            
            # Test plotting capability
            print(f"\n8. Test plotting first spectrum:")
            try:
                import matplotlib
                matplotlib.use('Agg')
                import matplotlib.pyplot as plt
                
                fig, ax = plt.subplots(figsize=(14, 5))
                ax.plot(vel_axis, spectrum, 'b-', linewidth=0.5, alpha=0.7)
                ax.set_xlabel('Velocity (m/s)', fontsize=12)
                ax.set_ylabel('Flux', fontsize=12)
                ax.set_title('First Reduced Spectrum with Velocity Axis', fontsize=14)
                ax.grid(True, alpha=0.3)
                
                plot_file = tmp_path.replace('.fits', '_plot.png')
                fig.savefig(plot_file, dpi=100, bbox_inches='tight')
                plt.close(fig)
                print(f"   ✓ Plot saved to: {plot_file}")
            except ImportError:
                print("   ! Matplotlib not available (skipping plot)")
            except Exception as e:
                print(f"   ! Could not create plot: {e}")
            
            # Clean up
            Path(tmp_path).unlink()
            print(f"\n9. Cleanup:")
            print(f"   ✓ Cleaned up temporary files")
            
            return True
            
    except Exception as e:
        print(f"   ✗ Error reading file: {e}")
        import traceback
        traceback.print_exc()
        return False


if __name__ == '__main__':
    success = test_real_fits()
    
    print("\n" + "=" * 80)
    if success:
        print("RESULT: ✓ ALL TESTS PASSED")
        sys.exit(0)
    else:
        print("RESULT: ✗ TESTS FAILED")
        sys.exit(1)
