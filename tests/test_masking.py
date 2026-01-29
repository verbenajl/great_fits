#!/usr/bin/env python
"""Test to see what's happening with channel masking."""

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
    
    # Check for NaN/Inf
    bad_channels = np.isnan(original) | np.isinf(original)
    good_channels = ~bad_channels
    
    print(f"Spectrum #{m51_idx}")
    print(f"  Total channels: {len(original)}")
    print(f"  Bad channels (NaN/Inf): {np.sum(bad_channels)}")
    print(f"  Good channels: {np.sum(good_channels)}")
    
    # Fit coefficients on all channels vs good channels only
    coeff_all = np.dot(corrector.components, original)
    coeff_good = np.dot(corrector.components[:, good_channels], original[good_channels])
    
    print(f"\nCoefficients fitted on ALL channels:")
    print(f"  {coeff_all}")
    
    print(f"\nCoefficients fitted on GOOD channels only:")
    print(f"  {coeff_good}")
    
    # Now correct using both approaches
    corrected_all = original.copy()
    for i in range(len(corrector.components)):
        corrected_all -= coeff_all[i] * corrector.components[i]
    
    corrected_good = original.copy()
    for i in range(len(corrector.components)):
        corrected_good[good_channels] -= coeff_good[i] * corrector.components[i, good_channels]
    
    print(f"\nCorrected using ALL channels:")
    print(f"  Min: {corrected_all.min():.4f}, Max: {corrected_all.max():.4f}, Mean: {corrected_all.mean():.6f}")
    
    print(f"\nCorrected using GOOD channels only (like the code does):")
    print(f"  Min: {corrected_good.min():.4f}, Max: {corrected_good.max():.4f}, Mean: {corrected_good.mean():.6f}")
    
    print(f"\nDifference between approaches:")
    print(f"  Max diff: {np.max(np.abs(corrected_all - corrected_good)):.6f}")
    
    # Check what the actual code does
    corrected_actual, details = corrector.apply_correction(original, good_channels=good_channels)
    print(f"\nActual code output:")
    print(f"  Min: {corrected_actual.min():.4f}, Max: {corrected_actual.max():.4f}, Mean: {corrected_actual.mean():.6f}")
    print(f"  Matches 'good channels only' approach: {np.allclose(corrected_actual, corrected_good)}")
