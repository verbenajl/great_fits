#!/usr/bin/env python3
"""
Debug script to test PCA correction and see what's being applied.
"""

import pickle
import numpy as np
from pathlib import Path
from astropy.io import fits
import sys
import logging

# Configure logging
logging.basicConfig(
    level=logging.DEBUG,
    format='%(levelname)s:%(name)s:%(message)s'
)
logger = logging.getLogger(__name__)

# Add src to path
sys.path.insert(0, '/home/verbena/software/oi_zeigt/src')

from oi_zeigt.pca_analysis.pca_correct_fits import PCACorrector

# Load decomposition
decomp_file = '/home/verbena/software/oi_zeigt/output/pca_components/decomposition_2017-02-01_GR_F367_LFAH_PX00_S_2017021_components.pkl'

logger.info(f"Loading decomposition from {decomp_file}")
with open(decomp_file, 'rb') as f:
    decomp = pickle.load(f)

logger.info(f"Decomposition keys: {decomp.keys()}")
logger.info(f"Components shape: {decomp['components'].shape}")
logger.info(f"Explained variance: {decomp['explained_variance']}")
logger.info(f"Explained variance ratio: {decomp['explained_variance_ratio']}")
logger.info(f"Mean spectrum shape: {decomp['mean_spectrum'].shape}")

# Create corrector
corrector = PCACorrector(decomp_file)

logger.info(f"\nCorrectorLoaded:")
logger.info(f"  Components shape: {corrector.components.shape}")
logger.info(f"  Variance ratio: {corrector.explained_variance_ratio}")

print("\n" + "="*80)
print("COMPONENT STATISTICS")
print("="*80)
for i, comp in enumerate(corrector.components):
    print(f"Component {i}: range=[{comp.min():.6f}, {comp.max():.6f}], "
          f"norm={np.linalg.norm(comp):.6f}, "
          f"mean={comp.mean():.6f}, std={comp.std():.6f}")

# Load a FITS file
fits_file = '/home/diskB/sofiaobsdata/m51/m51_central_region/kalibrate_d/M51CENTER_UNKNOWN_central_tile_CB_flight1_20170201_LFA_L1.fits'
logger.info(f"\nLoading FITS file: {fits_file}")

with fits.open(fits_file) as hdul:
    data = hdul[1].data
    logger.info(f"FITS table shape: {data.shape}")
    
    # Get first 5 spectra
    print("\n" + "="*80)
    print("TESTING FIRST 5 SPECTRA")
    print("="*80)
    
    for test_idx in range(min(5, len(data))):
        spec = data['SPECTRUM'][test_idx]
        
        print(f"\nSpectrum {test_idx}:")
        print(f"  Range: [{spec.min():.6f}, {spec.max():.6f}]")
        print(f"  Mean: {spec.mean():.6f}, Std: {spec.std():.6f}")
        
        # Check for bad channels
        bad_channels = np.isnan(spec) | np.isinf(spec) | (spec == 0)
        good_channels = ~bad_channels
        print(f"  Bad channels: {np.sum(bad_channels)}/{len(spec)}")
        
        # Fit coefficients
        coeff = corrector.fit_coefficients(spec, good_channels)
        print(f"  Fitted coefficients: {coeff}")
        print(f"  Coeff magnitudes: {np.abs(coeff)}")
        
        # Check noise ratios
        print(f"  Noise ratios:")
        for i, comp in enumerate(corrector.components):
            noise_ratio = corrector.get_noise_ratio(spec, comp)
            scaled_comp = coeff[i] * comp
            print(f"    Comp {i}: noise_ratio={noise_ratio:.4f}, coeff={coeff[i]:.6f}, "
                  f"scaled_range=[{scaled_comp.min():.6f}, {scaled_comp.max():.6f}]")
        
        # Apply correction WITHOUT noise ratio cutoff
        print(f"  Applying correction (no cutoff):")
        corrected, details = corrector.apply_correction(
            spec,
            good_channels=good_channels,
            cutoff_variance=None,
            cutoff_noise_ratio=None,
            verbose=False
        )
        
        diff = spec - corrected
        print(f"    Components used: {details['n_components_used']}")
        print(f"    Correction range: [{diff.min():.6f}, {diff.max():.6f}]")
        print(f"    Correction mean: {np.mean(np.abs(diff)):.6f}")
        
        # Apply correction WITH noise ratio cutoff=15
        print(f"  Applying correction (cutoff_noise_ratio=15):")
        corrected2, details2 = corrector.apply_correction(
            spec,
            good_channels=good_channels,
            cutoff_variance=None,
            cutoff_noise_ratio=15.0,
            verbose=False
        )
        
        diff2 = spec - corrected2
        print(f"    Components used: {details2['n_components_used']}")
        print(f"    Correction range: [{diff2.min():.6f}, {diff2.max():.6f}]")
        print(f"    Correction mean: {np.mean(np.abs(diff2)):.6f}")
        
        if details2['n_components_used'] == 0:
            print(f"    WARNING: No components passed the noise_ratio cutoff!")
