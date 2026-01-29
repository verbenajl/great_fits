#!/usr/bin/env python
"""Test correction on actual FITS data with limited spectra."""

import numpy as np
import sys
sys.path.insert(0, '/home/verbena/software/oi_zeigt/src')

from oi_zeigt.pca_analysis.pca_correct_fits import PCACorrector
from pathlib import Path
from astropy.io import fits

# Load corrector
decomp_file = Path('output/pca_components/decomposition_2017-02-01_GR_F367_2017021_components.pkl')
print(f"Loading decomposition from: {decomp_file}")
corrector = PCACorrector(str(decomp_file))

# Load FITS file
fits_file = Path('/home/diskB/sofiaobsdata/m51/m51_central_fits/fits_files_processed/reduced_data.fits')
print(f"\nLoading FITS from: {fits_file}")

with fits.open(fits_file) as hdul:
    data = hdul[1].data
    print(f"✓ Loaded {len(data)} spectra")
    
    # Get first 50 M51 spectra
    m51_mask = np.char.find(data['OBJECT'], 'M51') >= 0
    m51_indices = np.where(m51_mask)[0][:50]
    print(f"✓ Found {len(m51_indices)} M51 spectra (showing first 50)")
    
    # Get spectra and correct them
    original = data['SPECTRUM'][m51_indices].copy()
    corrected = np.zeros_like(original)
    
    for i, idx in enumerate(m51_indices):
        spectrum = data['SPECTRUM'][idx]
        corr, details = corrector.apply_correction(spectrum)
        corrected[i] = corr
        if i < 3:
            diff = np.mean(np.abs(spectrum - corr))
            print(f"  Spectrum {i}: correction_magnitude={diff:.4f}, components_used={details['n_components_used']}")
    
    # Statistics
    difference = original - corrected
    print(f"\nOverall Statistics for {len(m51_indices)} spectra:")
    print(f"  Original: min={original.min():.1f}, max={original.max():.1f}, mean={original.mean():.1f}")
    print(f"  Corrected: min={corrected.min():.1f}, max={corrected.max():.1f}, mean={corrected.mean():.1f}")
    print(f"  Difference: min={difference.min():.4f}, max={difference.max():.4f}, mean={difference.mean():.6f}")
    print(f"  Mean abs correction per spectrum: {np.mean(np.abs(difference)):.4f}")
    print(f"  Max abs correction per spectrum: {np.max(np.abs(difference)):.4f}")
    
    # Check if correction is visible
    if np.max(np.abs(difference)) > 0.01:
        print("\n✓ Correction is SIGNIFICANT and should be visible in plots!")
    else:
        print("\n✗ Correction is too small to see in plots")
