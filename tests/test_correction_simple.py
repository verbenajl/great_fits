#!/usr/bin/env python
"""Simple test to verify PCA correction is working."""

import numpy as np
import sys
sys.path.insert(0, '/home/verbena/software/oi_zeigt/src')

from oi_zeigt.pca_analysis.pca_correct_fits import PCACorrector
from pathlib import Path

# Load corrector
decomp_file = Path('output/pca_components/decomposition_2017-02-01_GR_F367_2017021_components.pkl')
print(f"Loading decomposition from: {decomp_file}")

corrector = PCACorrector(str(decomp_file))
print(f"✓ Loaded {len(corrector.components)} components")
print(f"  Variance explained: {sum(corrector.explained_variance_ratio)*100:.2f}%")

# Create a synthetic test spectrum with strong signal
test_spectrum = np.ones(700) * 100.0  # Base level
test_spectrum[300:350] = 150.0  # Add some feature
test_spectrum[400:450] = 120.0

print(f"\nTest spectrum:")
print(f"  Min: {test_spectrum.min():.2f}, Max: {test_spectrum.max():.2f}, Mean: {test_spectrum.mean():.2f}")

# Apply correction
corrected, details = corrector.apply_correction(test_spectrum)

print(f"\nAfter correction:")
print(f"  Min: {corrected.min():.2f}, Max: {corrected.max():.2f}, Mean: {corrected.mean():.2f}")

# Calculate difference
difference = test_spectrum - corrected
print(f"\nDifference (Original - Corrected):")
print(f"  Min: {difference.min():.2f}, Max: {difference.max():.2f}, Mean: {difference.mean():.4f}")
print(f"  Total correction: {np.sum(np.abs(difference)):.2f}")

# Check correction details
print(f"\nCorrection details:")
print(f"  Components used: {details['n_components_used']}")
print(f"  Coefficients: {details['coefficients']}")

# Verify actual subtraction
total_correction_manual = np.zeros_like(test_spectrum)
for i in range(len(corrector.components)):
    if i < len(details['coefficients']):
        coeff = details['coefficients'][i]
        comp = corrector.components[i]
        total_correction_manual += coeff * comp
        print(f"  Component {i}: coeff={coeff:.4f}, max_abs_scaled={np.max(np.abs(coeff*comp)):.4f}")

manual_corrected = test_spectrum - total_correction_manual
print(f"\nManual correction verification:")
print(f"  Max difference between manual and actual: {np.max(np.abs(manual_corrected - corrected)):.6f}")

if np.max(np.abs(difference)) > 0.001:
    print("\n✓ CORRECTION IS BEING APPLIED!")
else:
    print("\n✗ WARNING: No significant correction detected!")
