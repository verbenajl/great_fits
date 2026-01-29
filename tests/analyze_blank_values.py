#!/usr/bin/env python3
"""
Analyze spectra for unusual values that might be blank/missing markers.
Different astronomical software (GILDAS, CLASS, etc.) uses different blank values.
"""
import sys
from pathlib import Path
import numpy as np
from astropy.io import fits
from collections import Counter

# Add src to path
sys.path.insert(0, str(Path(__file__).parent / 'src'))

from oi_zeigt.reduction.core import detect_nan_channels

# Common blank value markers in astronomical software
COMMON_BLANK_VALUES = {
    'NaN (IEEE)': np.nan,
    'GILDAS/CLASS blank': -9.99e30,  # Common in GILDAS
    'ALMA blank': 0.0,
    'IRAC blank': -9999.0,
    'Spitzer blank': -1.0,
    'Very small values': 1e-30,
    'Zero': 0.0,
    'Negative infinity': -np.inf,
    'Positive infinity': np.inf,
}

def analyze_spectrum_for_blanks(spectrum, tolerance=1e-10):
    """Analyze a single spectrum for potential blank values."""
    analysis = {
        'total_values': len(spectrum),
        'nan_count': np.sum(np.isnan(spectrum)),
        'inf_count': np.sum(np.isinf(spectrum)),
        'negative_inf_count': np.sum(np.isneginf(spectrum)),
        'positive_inf_count': np.sum(np.isposinf(spectrum)),
        'zero_count': np.sum(spectrum == 0.0),
        'near_zero_count': np.sum(np.abs(spectrum) < tolerance),
        'negative_count': np.sum(spectrum < 0),
        'gildas_blank_count': np.sum(np.isclose(spectrum, -9.99e30, rtol=1e-5)),
        'very_small_count': np.sum((np.abs(spectrum) < 1e-20) & (spectrum != 0.0) & ~np.isnan(spectrum)),
        'min_finite_value': np.min(spectrum[np.isfinite(spectrum)]) if np.any(np.isfinite(spectrum)) else None,
        'max_value': np.max(spectrum[np.isfinite(spectrum)]) if np.any(np.isfinite(spectrum)) else None,
        'mean_value': np.mean(spectrum[np.isfinite(spectrum)]) if np.any(np.isfinite(spectrum)) else None,
        'std_value': np.std(spectrum[np.isfinite(spectrum)]) if np.any(np.isfinite(spectrum)) else None,
    }
    
    return analysis

def find_suspicious_values(spectrum, num_samples=20):
    """Find the most common values in a spectrum."""
    finite_values = spectrum[np.isfinite(spectrum)]
    
    if len(finite_values) == 0:
        return []
    
    # Get statistics
    value_counts = Counter(finite_values.round(6))  # Round to avoid float precision issues
    most_common = value_counts.most_common(num_samples)
    
    return most_common

def check_gildas_blank_pattern(spectrum):
    """Check if spectrum contains GILDAS-style blanks (-9.99e30 or similar)."""
    # GILDAS uses values close to -1e31
    suspicious = np.abs(spectrum) > 1e30
    if np.any(suspicious):
        return {
            'has_large_negative': True,
            'count': np.sum(suspicious),
            'values': np.unique(spectrum[suspicious])[:10]  # Show first 10 unique values
        }
    return {'has_large_negative': False, 'count': 0}

fits_file = Path('/home/diskB/sofiaobsdata/m51/m51_central_region/kalibrated_all/m51_central_tile.fits')

if not fits_file.exists():
    print(f"Error: Test file not found: {fits_file}")
    sys.exit(1)

print("="*80)
print("SPECTRAL VALUE ANALYSIS - Looking for Blank/Missing Value Markers")
print("="*80)
print(f"\nFile: {fits_file.name}")
print(f"Size: {fits_file.stat().st_size / (1024**3):.2f} GB\n")

with fits.open(fits_file) as hdul:
    # Collect statistics from multiple spectra
    all_stats = {
        'nan_count': 0,
        'zero_count': 0,
        'inf_count': 0,
        'gildas_blank_count': 0,
        'very_small_count': 0,
        'negative_count': 0,
        'spectra_with_nans': 0,
        'spectra_with_zeros': 0,
        'spectra_with_infs': 0,
        'spectra_with_gildas_blanks': 0,
    }
    
    hdu_count = 0
    spectrum_count = 0
    
    # Find first spectrum HDU
    for hdu_idx, hdu in enumerate(hdul):
        if hasattr(hdu, 'data') and hdu.data is not None:
            if 'SPECTRUM' in hdu.data.dtype.names:
                print(f"\nHDU {hdu_idx} ({hdu.name}): {len(hdu.data)} spectra")
                print("-" * 80)
                
                # Analyze first 20 spectra for detailed breakdown
                for i in range(min(20, len(hdu.data))):
                    spectrum = hdu.data[i]['SPECTRUM']
                    
                    # Get analysis
                    analysis = analyze_spectrum_for_blanks(spectrum)
                    
                    # Update global stats
                    all_stats['nan_count'] += analysis['nan_count']
                    all_stats['zero_count'] += analysis['zero_count']
                    all_stats['inf_count'] += analysis['inf_count']
                    all_stats['gildas_blank_count'] += analysis['gildas_blank_count']
                    all_stats['very_small_count'] += analysis['very_small_count']
                    all_stats['negative_count'] += analysis['negative_count']
                    
                    if analysis['nan_count'] > 0:
                        all_stats['spectra_with_nans'] += 1
                    if analysis['zero_count'] > 0:
                        all_stats['spectra_with_zeros'] += 1
                    if analysis['inf_count'] > 0:
                        all_stats['spectra_with_infs'] += 1
                    if analysis['gildas_blank_count'] > 0:
                        all_stats['spectra_with_gildas_blanks'] += 1
                    
                    spectrum_count += 1
                    
                    # Print details
                    print(f"\n  Row {i}:")
                    print(f"    Shape: {spectrum.shape}")
                    print(f"    NaN channels: {analysis['nan_count']} ({analysis['nan_count']/len(spectrum)*100:.2f}%)")
                    print(f"    Zero channels: {analysis['zero_count']} ({analysis['zero_count']/len(spectrum)*100:.2f}%)")
                    print(f"    Infinity channels: {analysis['inf_count']}")
                    print(f"    GILDAS blanks: {analysis['gildas_blank_count']}")
                    print(f"    Very small values (<1e-20): {analysis['very_small_count']}")
                    print(f"    Negative values: {analysis['negative_count']} ({analysis['negative_count']/len(spectrum)*100:.2f}%)")
                    
                    if analysis['min_finite_value'] is not None:
                        print(f"    Value range: {analysis['min_finite_value']:.6e} to {analysis['max_value']:.6e}")
                        print(f"    Mean: {analysis['mean_value']:.6e}, Std: {analysis['std_value']:.6e}")
                    
                    # Check for GILDAS-style blanks
                    gildas_check = check_gildas_blank_pattern(spectrum)
                    if gildas_check['has_large_negative']:
                        print(f"    ⚠️  GILDAS-STYLE BLANKS DETECTED: {gildas_check['count']} channels")
                        print(f"       Sample values: {gildas_check['values'][:3]}")
                    
                    # Show most common suspicious values
                    if analysis['zero_count'] > 0 or analysis['negative_count'] > 0:
                        common = find_suspicious_values(spectrum, num_samples=5)
                        if common:
                            print(f"    Most common values in spectrum:")
                            for value, count in common[:5]:
                                freq = count / len(spectrum) * 100
                                print(f"      {value:12.6e}: {count:6d} times ({freq:5.2f}%)")
                
                hdu_count += 1
                if hdu_count >= 1:  # Only analyze first HDU for speed
                    break

print("\n" + "="*80)
print("SUMMARY STATISTICS")
print("="*80)
print(f"Total spectra analyzed: {spectrum_count}")
print(f"\nBlank/Missing Value Markers Found:")
print(f"  NaN channels: {all_stats['nan_count']} total")
print(f"    In {all_stats['spectra_with_nans']} spectra")
print(f"  Zero channels: {all_stats['zero_count']} total")
print(f"    In {all_stats['spectra_with_zeros']} spectra")
print(f"  Infinity channels: {all_stats['inf_count']} total")
print(f"    In {all_stats['spectra_with_infs']} spectra")
print(f"  GILDAS blanks (≈-9.99e30): {all_stats['gildas_blank_count']} total")
print(f"    In {all_stats['spectra_with_gildas_blanks']} spectra")
print(f"  Very small values (<1e-20): {all_stats['very_small_count']} total")
print(f"  Negative values: {all_stats['negative_count']} total")

print("\n" + "="*80)
print("INTERPRETATION")
print("="*80)

if all_stats['spectra_with_nans'] > 0:
    print("✓ NaNs are being used as blank markers")
else:
    print("✗ No standard NaNs found - check for other blank markers")

if all_stats['spectra_with_gildas_blanks'] > 0:
    print("⚠️  GILDAS-style blanks detected! Need special handling")

if all_stats['spectra_with_zeros'] > spectrum_count * 0.5:
    print("⚠️  Many zero values - may be used as blanks or actual data")

print("\nNote: GILDAS uses -9.99e30 or similar for blanks")
print("      CLASS may use different conventions")
print("      Check your data source documentation!")

print("\n" + "="*80)
