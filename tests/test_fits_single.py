#!/usr/bin/env python
"""Test on actual FITS spectrum to see if correction makes sense."""

import numpy as np
import sys
sys.path.insert(0, '/home/verbena/software/oi_zeigt/src')

from oi_zeigt.pca_analysis.pca_correct_fits import PCACorrector
from pathlib import Path
from astropy.io import fits

# Load corrector
decomp_file = Path('output/pca_components/decomposition_2017-02-01_GR_F367_2017021_components.pkl')
corrector = PCACorrector(str(decomp_file))

# Load FITS file
fits_file = Path('/home/diskB/sofiaobsdata/m51/m51_central_fits/fits_files_processed/reduced_data.fits')

with fits.open(fits_file) as hdul:
    data = hdul[1].data
    
    # Get first M51 spectrum
    m51_mask = np.char.find(data['OBJECT'], 'M51') >= 0
    m51_idx = np.where(m51_mask)[0][0]
    
    original = data['SPECTRUM'][m51_idx].copy()
    corrected, details = corrector.apply_correction(original)
    difference = original - corrected
    
    print(f"Real FITS Spectrum #{m51_idx} - Object: {data['OBJECT'][m51_idx]}")
    print(f"\nOriginal spectrum:")
    print(f"  Min: {original.min():.2f}, Max: {original.max():.2f}, Mean: {original.mean():.2f}")
    print(f"  Std: {original.std():.2f}")
    
    print(f"\nAfter correction:")
    print(f"  Min: {corrected.min():.2f}, Max: {corrected.max():.2f}, Mean: {corrected.mean():.2f}")
    print(f"  Std: {corrected.std():.2f}")
    
    print(f"\nDifference (Original - Corrected):")
    print(f"  Min: {difference.min():.4f}, Max: {difference.max():.4f}, Mean: {difference.mean():.6f}")
    print(f"  Std: {difference.std():.4f}")
    print(f"  Mean abs: {np.mean(np.abs(difference)):.4f}")
    
    print(f"\nComponents used: {details['n_components_used']}")
    print(f"Coefficients: {details['coefficients']}")
    
    # Check if difference is reasonable
    print(f"\nAnalysis:")
    print(f"  Correction is {np.mean(np.abs(difference)):.1f}% of signal range")
    print(f"  Correction is {(np.mean(np.abs(difference)) / original.std() * 100):.1f}% of signal std")
    
    # Check min/max changes
    min_change = original.min() - corrected.min()
    max_change = original.max() - corrected.max()
    print(f"\n  Min changed by: {min_change:.4f} (original={original.min():.2f}, corrected={corrected.min():.2f})")
    print(f"  Max changed by: {max_change:.4f} (original={original.max():.2f}, corrected={corrected.max():.2f})")
    
    if min_change < 0 and max_change < 0:
        print("  ✓ Both min and max decreased (removed signal)")
    elif min_change > 0 and max_change > 0:
        print("  ✓ Both min and max increased (added signal)")
    else:
        print("  ✗ Min and max changed in opposite directions (asymmetric correction)")
