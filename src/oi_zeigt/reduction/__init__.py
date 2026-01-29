"""
Data reduction module for OI ZEIGT.

This package contains functions for reducing and processing OI GREAT data.
"""

"""Convenience imports for the reduction subpackage."""

from .core import (
    detect_nan_channels,
    detect_missing_channels,
    baseline_subtract,
    apply_baseline_to_hdul,
    apply_baseline_from_config,
    reduce_spectra,
    reduce_spectra_from_config,
    average_spectra,
    average_spectra_from_config,
    filter_and_save_fits,
)

__all__ = [
    'detect_nan_channels',
    'detect_missing_channels',
    'baseline_subtract',
    'apply_baseline_to_hdul',
    'apply_baseline_from_config',
    'reduce_spectra',
    'reduce_spectra_from_config',
    'average_spectra',
    'average_spectra_from_config',
    'filter_and_save_fits',
]