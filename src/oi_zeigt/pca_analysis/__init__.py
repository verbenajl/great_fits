"""
Principal Component Analysis (PCA) module for spectral datacube analysis.

This module provides tools for analyzing spectral datacubes using PCA techniques,
including:

Phase 1 (Foundation):
- Configuration management
- Utilities and helpers
- FITS file indexing and loading
- Science line detection
- Error handling

Phase 2 (Decomposition):
- Spectrum preparation and normalization
- PCA decomposition using scikit-learn
- Data loading and quality filtering
- Results storage and serialization

Phase 3 (Correction):
- Spectrum projection onto components
- Correction calculation via least-squares
- Receiver variation removal
- Science line preservation
- Corrected spectra output

Additional capabilities:
- Spectral PCA decomposition
- Spatial correlation analysis
- Noise characterization
- Component visualization
"""

from .core import (
    perform_spectral_pca,
    reconstruct_from_pca,
    get_component_stats,
    filter_datacube_by_pca,
)

# Phase 1 foundation modules
from .errors import (
    PCAError,
    ConfigurationError,
    MissingConfigurationError,
    InvalidConfigurationError,
    NoDataFoundError,
    NoSpectraFoundError,
    InsufficientDataError,
    InsufficientSpectraError,
    LineDetectionError,
    FITSError,
    DecompositionError,
    CorrectionError,
    BaselineError,
)

from .config import ConfigLoader, create_default_config_template

from .utilities import (
    fit_baseline_poly,
    create_channel_mask,
    combine_masks,
    apply_mask_to_spectrum,
    normalize_spectrum,
    denormalize_spectrum,
    check_spectrum_validity,
    get_spectrum_stats,
    rolling_window,
    ensure_directory,
    safe_pickle_save,
    safe_pickle_load,
)

from .fits_indexing import (
    FITSIndexer,
    load_spectral_data,
    save_spectral_data,
)

from .line_detection import (
    detect_science_line_waterfall,
    detect_science_line_max_intensity,
    detect_science_line_centroid,
    refine_line_width,
    LineDetectionConfig,
    create_science_line_mask,
    create_artifact_mask,
    validate_line_detection,
)

# Phase 2 decomposition module
from .decompose import (
    DecompositionConfig,
    SpectrumPreparator,
    DecompositionDataLoader,
    PCADecomposer,
    DecompositionResult,
    decompose_spectra,
)

# Phase 3 correction module
from .correct import (
    CorrectionConfig,
    SpectrumProjector,
    CorrectionCalculator,
    ScienceDataLoader,
    CorrectionResult,
    correct_spectra,
)

__all__ = [
    # Core functions
    'perform_spectral_pca',
    'reconstruct_from_pca',
    'get_component_stats',
    'filter_datacube_by_pca',
    
    # Error classes
    'PCAError',
    'ConfigurationError',
    'MissingConfigurationError',
    'InvalidConfigurationError',
    'NoDataFoundError',
    'NoSpectraFoundError',
    'InsufficientDataError',
    'InsufficientSpectraError',
    'LineDetectionError',
    'FITSError',
    'DecompositionError',
    'CorrectionError',
    'BaselineError',
    
    # Configuration
    'ConfigLoader',
    'create_default_config_template',
    
    # Utilities
    'fit_baseline_poly',
    'create_channel_mask',
    'combine_masks',
    'apply_mask_to_spectrum',
    'normalize_spectrum',
    'denormalize_spectrum',
    'check_spectrum_validity',
    'get_spectrum_stats',
    'rolling_window',
    'ensure_directory',
    'safe_pickle_save',
    'safe_pickle_load',
    
    # FITS indexing
    'FITSIndexer',
    'load_spectral_data',
    'save_spectral_data',
    
    # Line detection
    'detect_science_line_waterfall',
    'detect_science_line_max_intensity',
    'detect_science_line_centroid',
    'refine_line_width',
    'LineDetectionConfig',
    'create_science_line_mask',
    'create_artifact_mask',
    'validate_line_detection',
    
    # Phase 2 Decomposition
    'DecompositionConfig',
    'SpectrumPreparator',
    'DecompositionDataLoader',
    'PCADecomposer',
    'DecompositionResult',
    'decompose_spectra',
    
    # Phase 3 Correction
    'CorrectionConfig',
    'SpectrumProjector',
    'CorrectionCalculator',
    'ScienceDataLoader',
    'CorrectionResult',
    'correct_spectra',
]
