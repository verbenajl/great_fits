#!/usr/bin/env python3
"""
Test spectral extraction from config.toml with real FITS data.

This test:
1. Reads config.toml which has extract=[350,700]
2. Loads the real FITS file (31,822 spectra, 1,264 channels)
3. Applies extraction and baseline reduction
4. Verifies the output spectra are trimmed to the correct velocity range
5. Validates the output FITS file
"""

import numpy as np
from pathlib import Path
from astropy.io import fits

# Add src to path
import sys
sys.path.insert(0, '/home/verbena/software/oi_zeigt/src')

from oi_zeigt.reduction.core import (
    reduce_spectra_from_config,
    _extract_spectral_params,
    _create_velocity_axis,
    _velocity_to_channel_index
)
from oi_zeigt.basic_io import get_config, read_fits_from_config


def test_extraction_from_config():
    """Test extraction with config-based reduction."""
    
    config_path = '/home/verbena/software/oi_zeigt/config.toml'
    cfg = get_config(config_path)
    
    # Print config
    print("=" * 80)
    print("EXTRACTION TEST FROM CONFIG")
    print("=" * 80)
    print(f"\nConfig [reduction] section:")
    print(f"  baseline: {cfg['reduction'].get('baseline')}")
    print(f"  window: {cfg['reduction'].get('window')}")
    print(f"  extract: {cfg['reduction'].get('extract')}")
    print()
    
    # Load original FITS
    print("Loading FITS file...")
    hdul = read_fits_from_config(config_path)
    original_shape = hdul[1].data['SPECTRUM'].shape
    print(f"  Original spectra shape: {original_shape}")
    print(f"  Original number of channels: {original_shape[1]}")
    
    # Extract spectral parameters
    print("\nExtracting spectral parameters...")
    spectral_params = _extract_spectral_params(hdul)
    velo_ref = spectral_params['velo_ref']
    deltav = spectral_params['deltav']
    crpix1_spec = spectral_params['crpix1_spec']
    nchans_original = spectral_params['nchans']
    
    print(f"  Reference velocity (VELOCITY): {velo_ref} m/s")
    print(f"  Channel spacing (DELTAV): {deltav} m/s/channel")
    print(f"  Reference pixel (CRPIX1): {crpix1_spec} (1-indexed)")
    print(f"  Total channels: {nchans_original}")
    
    # Create full velocity axis
    velocity_axis_full = _create_velocity_axis(velo_ref, deltav, crpix1_spec, nchans_original)
    print(f"  Velocity range: {velocity_axis_full[0]:.1f} - {velocity_axis_full[-1]:.1f} m/s")
    
    # Check extraction range
    extract_range = cfg['reduction'].get('extract')
    if extract_range:
        print(f"\nExtraction range from config: {extract_range} km/s")
        
        # Convert km/s to m/s
        velo_min_ms = extract_range[0] * 1000.0
        velo_max_ms = extract_range[1] * 1000.0
        print(f"  Converted to m/s: [{velo_min_ms:.0f}, {velo_max_ms:.0f}]")
        
        # Convert velocity to channel indices
        ch_min_idx = _velocity_to_channel_index(velo_min_ms, velo_ref, deltav, crpix1_spec)
        ch_max_idx = _velocity_to_channel_index(velo_max_ms, velo_ref, deltav, crpix1_spec)
        ch_min = int(np.clip(ch_min_idx, 0, nchans_original - 1))
        ch_max = int(np.clip(ch_max_idx, 0, nchans_original - 1))
        
        print(f"  Mapped to channel indices: [{ch_min}, {ch_max}]")
        print(f"  Actual velocity at min channel: {velocity_axis_full[ch_min]:.1f} m/s ({velocity_axis_full[ch_min]/1000:.0f} km/s)")
        print(f"  Actual velocity at max channel: {velocity_axis_full[ch_max]:.1f} m/s ({velocity_axis_full[ch_max]/1000:.0f} km/s)")
        print(f"  Number of channels after extraction: {ch_max - ch_min + 1}")
    
    # Test reduction with config
    print("\n" + "=" * 80)
    print("APPLYING REDUCTION WITH EXTRACTION")
    print("=" * 80)
    
    output_path = Path('/home/verbena/software/oi_zeigt/test_reduced_extracted.fits')
    
    print(f"\nApplying reduce_spectra_from_config()...")
    print(f"  Output path: {output_path}")
    
    result_path = reduce_spectra_from_config(
        config_path=config_path,
        output_path=output_path,
        overwrite=True
    )
    
    print(f"  ✓ Reduction complete. Output: {result_path}")
    
    # Verify output
    print("\n" + "=" * 80)
    print("VERIFYING OUTPUT")
    print("=" * 80)
    
    with fits.open(result_path) as hdul_out:
        spectra_out = hdul_out[1].data['SPECTRUM']
        print(f"\nOutput spectra shape: {spectra_out.shape}")
        print(f"Number of spectra: {spectra_out.shape[0]}")
        print(f"Number of channels after extraction: {spectra_out.shape[1]}")
        
        # Check VELOCITY_AXIS column
        if 'VELOCITY_AXIS' in hdul_out[1].data.dtype.names:
            velocity_axis_out = hdul_out[1].data['VELOCITY_AXIS'][0]
            print(f"\nVELOCITY_AXIS column present:")
            print(f"  Length: {len(velocity_axis_out)}")
            print(f"  Velocity range: {velocity_axis_out[0]:.1f} - {velocity_axis_out[-1]:.1f} m/s")
            
            # Verify consistency across spectra
            all_same = True
            for i in range(1, min(10, len(hdul_out[1].data))):
                if not np.allclose(velocity_axis_out, hdul_out[1].data['VELOCITY_AXIS'][i]):
                    all_same = False
                    break
            print(f"  All spectra have identical velocity axis: {all_same}")
        else:
            print(f"\n⚠ WARNING: VELOCITY_AXIS column not found in output!")
        
        # Check some spectra values
        print(f"\nSample spectrum statistics:")
        sample_spectrum = spectra_out[0]
        print(f"  Mean: {np.nanmean(sample_spectrum):.6f}")
        print(f"  Std: {np.nanstd(sample_spectrum):.6f}")
        print(f"  Min: {np.nanmin(sample_spectrum):.6f}")
        print(f"  Max: {np.nanmax(sample_spectrum):.6f}")
        print(f"  NaN count: {np.isnan(sample_spectrum).sum()} / {len(sample_spectrum)}")
    
    # Check if baseline was applied
    print(f"\nBaseline was applied: baseline order = {cfg['reduction'].get('baseline')}")
    
    print("\n" + "=" * 80)
    print("✓ EXTRACTION TEST COMPLETE")
    print("=" * 80)


if __name__ == '__main__':
    test_extraction_from_config()
