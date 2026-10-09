"""
Statistics and quality analysis module for OI ZEIGT.

This module provides tools for analyzing spectral data quality,
computing statistical metrics, and generating quality reports.
"""

from .quality import (
    get_spechistogram,
    get_rmsratio_histogram,
    rmsratio_statistics,
    ripple_ratio,
)

__all__ = [
    'get_spechistogram',
    'get_rmsratio_histogram',
    'rmsratio_statistics',
    'ripple_ratio',
]
