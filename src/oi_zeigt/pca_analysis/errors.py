"""
Custom exception classes for PCA analysis module.

Provides specific error types for different failure modes in PCA decomposition
and correction workflows.
"""


class PCAError(Exception):
    """Base exception for all PCA analysis errors."""
    pass


class ConfigurationError(PCAError):
    """Raised when configuration is missing or invalid."""
    pass


class MissingConfigurationError(ConfigurationError):
    """Raised when a required configuration parameter is missing."""
    pass


class InvalidConfigurationError(ConfigurationError):
    """Raised when a configuration parameter has an invalid value."""
    pass


class NoDataFoundError(PCAError):
    """Raised when no matching data is found for the given filters."""
    pass


class NoSpectraFoundError(NoDataFoundError):
    """Raised when no spectra match the filtering criteria."""
    pass


class InsufficientDataError(PCAError):
    """Raised when there is not enough data to perform the operation."""
    pass


class LineDetectionError(PCAError):
    """Raised when science line detection fails."""
    pass


class FITSError(PCAError):
    """Raised when there is an error reading/writing FITS files."""
    pass


class DecompositionError(PCAError):
    """Raised when PCA decomposition fails."""
    pass


class InsufficientSpectraError(DecompositionError):
    """Raised when there aren't enough spectra for the requested components."""
    pass


class CorrectionError(PCAError):
    """Raised when PCA correction fails."""
    pass


class BaselineError(PCAError):
    """Raised when baseline fitting fails."""
    pass
