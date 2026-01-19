"""
Mapping and gridding functions for OI ZEIGT spectral data.

This module provides tools for creating proper WCS-based spatial maps from GREAT observations,
including cygrid-based gridding for optimal Gaussian kernel weighting.

Available functions:
- grid_spectra_with_cygrid: Core WCS-based gridding with cygrid
- create_wcs_header: Create WCS-compliant FITS headers
- create_map_from_column: Create spatial map from a scalar column
- create_integrated_map: Create spatial map from integrated intensity
"""

from .gridding import (
    grid_spectra_with_cygrid,
    create_wcs_header,
    create_map_from_column,
    create_integrated_map,
)

__all__ = [
    'grid_spectra_with_cygrid',
    'create_wcs_header',
    'create_map_from_column',
    'create_integrated_map',
]
