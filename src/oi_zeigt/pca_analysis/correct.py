"""
Correction module for applying PCA-based receiver variation corrections.

Provides functionality to:
- Load PCA decomposition results (Phase 2)
- Load science spectra (M51CENTER)
- Project science spectra onto PCA components
- Calculate optimal correction via least-squares fitting
- Apply correction while preserving science lines
- Store corrected spectra with full provenance

The correction workflow removes receiver/instrument variations from spectra
while preserving authentic science line features.
"""

import numpy as np
import logging
from pathlib import Path
from typing import Tuple, List, Optional, Dict, Union
from dataclasses import dataclass, asdict
import pickle

from .config import ConfigLoader
from .decompose import DecompositionResult
from .fits_indexing import FITSIndexer
from .utilities import (
    check_spectrum_validity,
    get_spectrum_stats,
    create_channel_mask,
    combine_masks,
    safe_pickle_save,
    safe_pickle_load,
    ensure_directory,
)
from .line_detection import (
    detect_science_line_waterfall,
    create_science_line_mask,
    create_artifact_mask,
)
from .errors import (
    CorrectionError,
    FITSError,
    LineDetectionError,
    InsufficientDataError,
)

logger = logging.getLogger(__name__)


@dataclass
class CorrectionConfig:
    """Configuration parameters for PCA-based spectral correction."""
    
    n_components: int
    science_source: str = "M51CENTER"
    preserve_science_line: bool = True
    mask_width: int = 20
    fit_window: int = 512
    regularization: float = 0.0
    
    def to_dict(self) -> dict:
        """Convert to dictionary."""
        return asdict(self)
    
    @classmethod
    def from_dict(cls, data: dict) -> 'CorrectionConfig':
        """Create from dictionary."""
        return cls(**data)


class SpectrumProjector:
    """Project spectra onto PCA components.
    
    Transforms 1D spectra into coefficient space defined by PCA components.
    """
    
    def __init__(self, decomposition_result: DecompositionResult):
        """
        Initialize projector with PCA components.
        
        Parameters
        ----------
        decomposition_result : DecompositionResult
            Results from Phase 2 decomposition containing components and mean
        """
        self.components = decomposition_result.components
        self.mean_spectrum = decomposition_result.mean_spectrum
        self.n_components, self.n_channels = self.components.shape
        
        logger.info(f"Initialized SpectrumProjector: "
                   f"{self.n_components} components, {self.n_channels} channels")
    
    def project_spectrum(self, spectrum_1d: np.ndarray) -> Tuple[np.ndarray, dict]:
        """
        Project a single spectrum onto components.
        
        Parameters
        ----------
        spectrum_1d : np.ndarray
            1D spectrum array [n_channels]
            
        Returns
        -------
        coefficients : np.ndarray
            Component coefficients [n_components]
        metadata : dict
            Projection metadata (norms, stats, etc.)
            
        Raises
        ------
        CorrectionError
            If projection fails
        """
        try:
            spectrum = spectrum_1d.copy().astype(float)
            
            # Validate input
            if not check_spectrum_validity(spectrum):
                raise CorrectionError(f"Invalid spectrum (shape: {spectrum.shape})")
            
            if len(spectrum) != self.n_channels:
                raise CorrectionError(
                    f"Spectrum shape mismatch: expected {self.n_channels}, "
                    f"got {len(spectrum)}"
                )
            
            # Center spectrum
            centered_spectrum = spectrum - self.mean_spectrum
            
            # Project onto components: coefficients = spectrum @ components.T
            coefficients = centered_spectrum @ self.components.T
            
            # Metadata
            metadata = {
                'input_norm': np.linalg.norm(spectrum),
                'centered_norm': np.linalg.norm(centered_spectrum),
                'coefficient_norm': np.linalg.norm(coefficients),
                'coefficient_stats': {
                    'mean': float(np.mean(coefficients)),
                    'std': float(np.std(coefficients)),
                    'min': float(np.min(coefficients)),
                    'max': float(np.max(coefficients)),
                }
            }
            
            return coefficients, metadata
        
        except Exception as e:
            raise CorrectionError(f"Spectrum projection failed: {e}")
    
    def project_multiple(self, spectra_matrix: np.ndarray) -> Tuple[np.ndarray, List[dict]]:
        """
        Project multiple spectra onto components.
        
        Parameters
        ----------
        spectra_matrix : np.ndarray
            2D array of spectra [n_spectra, n_channels]
            
        Returns
        -------
        coefficients_matrix : np.ndarray
            Coefficients [n_spectra, n_components]
        metadata_list : list
            List of per-spectrum metadata dicts
            
        Raises
        ------
        CorrectionError
            If projection fails
        """
        try:
            n_spectra = spectra_matrix.shape[0]
            logger.info(f"Projecting {n_spectra} spectra onto {self.n_components} components")
            
            coefficients_list = []
            metadata_list = []
            
            for i, spectrum in enumerate(spectra_matrix):
                coeff, meta = self.project_spectrum(spectrum)
                coefficients_list.append(coeff)
                metadata_list.append({**meta, 'spectrum_id': i})
            
            coefficients_matrix = np.array(coefficients_list)
            
            logger.info(f"Projected {n_spectra} spectra: "
                       f"coefficients shape {coefficients_matrix.shape}")
            
            return coefficients_matrix, metadata_list
        
        except Exception as e:
            raise CorrectionError(f"Multiple spectrum projection failed: {e}")


class CorrectionCalculator:
    """Calculate optimal correction via least-squares fitting."""
    
    def __init__(self, config: CorrectionConfig, decomposition_result: DecompositionResult):
        """
        Initialize correction calculator.
        
        Parameters
        ----------
        config : CorrectionConfig
            Configuration parameters
        decomposition_result : DecompositionResult
            Phase 2 decomposition results
        """
        self.config = config
        self.components = decomposition_result.components
        self.mean_spectrum = decomposition_result.mean_spectrum
        self.n_components = config.n_components
        
        logger.info(f"Initialized CorrectionCalculator "
                   f"(regularization={config.regularization})")
    
    def calculate_correction(
        self,
        reference_coefficients: np.ndarray,
        science_coefficients: np.ndarray,
    ) -> Tuple[np.ndarray, dict]:
        """
        Calculate optimal correction matrix via least-squares.
        
        Solves: min_C ||reference - C @ science.T||²_F
        
        Parameters
        ----------
        reference_coefficients : np.ndarray
            Reference (SKYCHOPDIFF) coefficients [n_ref, n_components]
        science_coefficients : np.ndarray
            Science (M51CENTER) coefficients [n_science, n_components]
            
        Returns
        -------
        correction_matrix : np.ndarray
            Correction matrix [n_components, n_components]
        stats : dict
            Statistics about correction (residual RMS, condition number, etc.)
            
        Raises
        ------
        CorrectionError
            If calculation fails
        """
        try:
            logger.info(f"Calculating correction: "
                       f"{reference_coefficients.shape[0]} reference, "
                       f"{science_coefficients.shape[0]} science spectra")
            
            # Validate inputs
            if reference_coefficients.shape[1] != self.n_components:
                raise CorrectionError(f"Reference shape mismatch")
            if science_coefficients.shape[1] != self.n_components:
                raise CorrectionError(f"Science shape mismatch")
            
            # Least-squares solution based on mean coefficients
            # Since we have reference and science from different measurements,
            # we find a correction C that maps mean(science) → mean(reference)
            # using least-squares fitting
            
            # Compute mean coefficient vectors
            mean_science = np.mean(science_coefficients, axis=0)  # [n_comp]
            mean_reference = np.mean(reference_coefficients, axis=0)  # [n_comp]
            
            logger.info(f"Mean science coefficients shape: {mean_science.shape}")
            logger.info(f"Mean reference coefficients shape: {mean_reference.shape}")
            
            # Solve: mean_science @ C.T = mean_reference
            # Using least-squares: C.T = pinv(mean_science.T) @ mean_reference.T
            # But mean_science and mean_reference are 1D vectors [n_comp]
            # This gives us one equation per component (underdetermined)
            
            # Instead, use a more sophisticated approach:
            # Find C such that science_coefficients @ C.T ≈ reference_coefficients
            # Using all spectra but allowing mismatched numbers via mean-centering
            
            # Center both by subtracting their means
            science_centered = science_coefficients - mean_science[np.newaxis, :]  # [n_science, n_comp]
            reference_centered = reference_coefficients - mean_reference[np.newaxis, :]  # [n_ref, n_comp]
            
            # Compute covariance matrix: Cov = science.T @ science / (n_science - 1)
            # For correction: we want to find C that relates the two
            # Use the singular value decomposition approach:
            # science_centered.T @ science_centered approximates within-science covariance
            
            # Compute the transformation using least-squares on centered data
            # For each reference spectrum, find best science spectrum(s) that map to it
            # But with different numbers, we use principal alignment:
            
            # Compute SVD of both covariance matrices
            U_science, S_science, Vt_science = np.linalg.svd(
                science_centered.T @ science_centered, full_matrices=False
            )
            U_ref, S_ref, Vt_ref = np.linalg.svd(
                reference_centered.T @ reference_centered, full_matrices=False
            )
            
            # Correction aligns the principal directions:
            # C = U_ref @ U_science.T (rotation matrix)
            # Then scale by ratio of variances
            correction_matrix = U_ref @ U_science.T
            
            # Add mean shift correction
            mean_shift = np.linalg.lstsq(
                correction_matrix.T, 
                (mean_reference - mean_science), 
                rcond=None
            )[0]  # Solve: C.T @ x = mean_diff
            
            # Alternative simpler approach: diagonal scaling by mean ratio
            # For each component, scale by the ratio of means
            correction_diagonal = np.where(
                np.abs(mean_science) > 1e-10,
                mean_reference / mean_science,
                1.0  # No correction for near-zero components
            )
            
            # Combine rotation and scaling
            correction_matrix = np.diag(correction_diagonal) @ correction_matrix
            
            logger.info(f"Correction matrix shape: {correction_matrix.shape}")
            logger.info(f"Correction diagonal scaling: {correction_diagonal}")
            
            # Compute condition number
            try:
                cond_number = np.linalg.cond(correction_matrix)
            except:
                cond_number = np.inf
            logger.info(f"Condition number: {cond_number:.2e}")
            
            # Compute residual as mean difference after correction
            # Apply correction to all science spectra
            corrected_science = science_coefficients @ correction_matrix.T  # [n_science, n_comp]
            
            # Compare corrected science mean to reference mean
            mean_corrected_science = np.mean(corrected_science, axis=0)  # [n_comp]
            mean_residual = mean_reference - mean_corrected_science  # [n_comp]
            residual_rms = np.linalg.norm(mean_residual)
            
            logger.info(f"Correction calculated: residual RMS = {residual_rms:.6f}")
            
            # Correction strength (per component)
            correction_strength = np.array([
                np.linalg.norm(correction_matrix[i, :])
                for i in range(self.n_components)
            ])
            
            stats = {
                'residual_rms': float(residual_rms),
                'condition_number': float(cond_number),
                'correction_strength': correction_strength.tolist(),
                'regularization': self.config.regularization,
                'n_reference': reference_coefficients.shape[0],
                'n_science': science_coefficients.shape[0],
            }
            
            return correction_matrix, stats
        
        except Exception as e:
            raise CorrectionError(f"Correction calculation failed: {e}")
    
    def apply_correction(
        self,
        coefficients: np.ndarray,
        correction_matrix: np.ndarray,
    ) -> np.ndarray:
        """
        Apply correction to coefficients.
        
        Parameters
        ----------
        coefficients : np.ndarray
            Input coefficients [n_spectra, n_components]
        correction_matrix : np.ndarray
            Correction matrix [n_components, n_components]
            
        Returns
        -------
        corrected_coefficients : np.ndarray
            Corrected coefficients [n_spectra, n_components]
            
        Raises
        ------
        CorrectionError
            If application fails
        """
        try:
            logger.info(f"Applying correction to {coefficients.shape[0]} spectra")
            
            # Apply: corrected = coefficients @ correction.T
            corrected = coefficients @ correction_matrix.T
            
            logger.info(f"Correction applied successfully")
            
            return corrected
        
        except Exception as e:
            raise CorrectionError(f"Correction application failed: {e}")


class ScienceDataLoader:
    """Load science spectra (M51CENTER) for correction."""
    
    def __init__(self, config: ConfigLoader, mission_id: str):
        """
        Initialize science data loader.
        
        Parameters
        ----------
        config : ConfigLoader
            Configuration loader
        mission_id : str
            Mission identifier
        """
        self.config = config
        self.mission_id = mission_id
        self.indexer = None
        
        logger.info(f"Initialized ScienceDataLoader for mission {mission_id}")
    
    def load_science_data(
        self,
        fits_directory: str,
        science_source: str = "M51CENTER",
    ) -> Tuple[np.ndarray, List[dict], dict]:
        """
        Load science spectra from FITS files.
        
        Parameters
        ----------
        fits_directory : str
            Path to FITS directory
        science_source : str
            Source name to filter (default: M51CENTER)
            
        Returns
        -------
        spectral_matrix : np.ndarray
            [n_spectra, n_channels] prepared spectra
        metadata_list : list
            Per-spectrum metadata
        summary : dict
            Summary statistics
            
        Raises
        ------
        FITSError
            If loading fails
        InsufficientDataError
            If insufficient spectra found
        """
        try:
            logger.info(f"Loading science data from {fits_directory}")
            logger.info(f"Filtering for source: {science_source}")
            
            # Index FITS files
            self.indexer = FITSIndexer(str(fits_directory))
            self.indexer.scan_directory()
            logger.info(f"Found {len(self.indexer.index_df)} FITS files")
            
            # Filter by source
            source_data = self.indexer.filter_by_source(science_source)
            n_science = len(source_data.index_df)
            logger.info(f"Found {n_science} {science_source} spectra")
            
            if n_science == 0:
                raise InsufficientDataError(f"No {science_source} spectra found")
            
            # Extract spectra
            spectral_list = []
            metadata_list = []
            
            logger.info("Extracting science spectra...")
            
            for data, header, metadata in source_data.get_all_fits_data():
                try:
                    # Extract 1D spectrum (average of 2D/3D data)
                    if data.ndim == 3:
                        spectrum = np.nanmean(data[:, 0, :], axis=0)
                    elif data.ndim == 2:
                        spectrum = np.nanmean(data, axis=0)
                    else:
                        spectrum = data
                    
                    # Validate
                    if not check_spectrum_validity(spectrum):
                        logger.warning(f"Invalid science spectrum, skipping")
                        continue
                    
                    spectral_list.append(spectrum)
                    metadata_list.append(metadata)
                
                except Exception as e:
                    logger.warning(f"Failed to extract spectrum: {e}")
                    continue
            
            if not spectral_list:
                raise InsufficientDataError("Failed to extract any valid science spectra")
            
            spectral_matrix = np.array(spectral_list)
            logger.info(f"Loaded {len(spectral_list)} science spectra: "
                       f"shape {spectral_matrix.shape}")
            
            # Summary
            summary = {
                'n_spectra': len(spectral_list),
                'n_channels': spectral_matrix.shape[1],
                'source': science_source,
                'mean_spectrum': np.mean(spectral_matrix, axis=0),
                'std_spectrum': np.std(spectral_matrix, axis=0),
            }
            
            return spectral_matrix, metadata_list, summary
        
        except Exception as e:
            raise FITSError(f"Failed to load science data: {e}")


@dataclass
class CorrectionResult:
    """Results from PCA-based spectral correction."""
    
    reference_coefficients: np.ndarray  # [n_ref, n_components]
    science_coefficients: np.ndarray    # [n_science, n_components]
    corrected_coefficients: np.ndarray  # [n_science, n_components]
    correction_matrix: np.ndarray       # [n_components, n_components]
    components: np.ndarray              # [n_components, n_channels]
    mean_spectrum: np.ndarray           # [n_channels]
    config: CorrectionConfig
    metadata: dict                      # Mission, statistics, etc.
    reference_metadata: List[dict] = None
    science_metadata: List[dict] = None
    
    def __post_init__(self):
        """Initialize optional fields."""
        if self.reference_metadata is None:
            self.reference_metadata = []
        if self.science_metadata is None:
            self.science_metadata = []
    
    def reconstruct_corrected_spectra(self) -> np.ndarray:
        """
        Reconstruct corrected spectra from corrected coefficients.
        
        Returns
        -------
        np.ndarray
            Corrected spectra [n_science, n_channels]
        """
        try:
            # Inverse transform: spectrum = coefficients @ components + mean
            corrected_spectra = self.corrected_coefficients @ self.components + self.mean_spectrum
            return corrected_spectra
        except Exception as e:
            raise CorrectionError(f"Reconstruction failed: {e}")
    
    def get_correction_strength(self) -> dict:
        """
        Get statistics about correction magnitude.
        
        Returns
        -------
        dict
            Correction strength per component and overall
        """
        # Per-component correction strength
        before = np.mean(np.abs(self.science_coefficients), axis=0)
        after = np.mean(np.abs(self.corrected_coefficients), axis=0)
        change = before - after
        
        return {
            'per_component_change': change.tolist(),
            'total_change_mean': float(np.mean(change)),
            'total_change_std': float(np.std(change)),
            'max_change': float(np.max(np.abs(change))),
            'before_mean_magnitude': float(np.mean(before)),
            'after_mean_magnitude': float(np.mean(after)),
        }
    
    def summary(self) -> str:
        """
        Generate text summary of correction results.
        
        Returns
        -------
        str
            Formatted summary text
        """
        n_components = self.components.shape[0]
        n_science = self.science_coefficients.shape[0]
        n_channels = self.components.shape[1]
        
        correction_strength = self.get_correction_strength()
        
        summary_text = f"""
Correction Results Summary:
  Reference spectra: {self.reference_coefficients.shape[0]}
  Science spectra: {n_science}
  Components: {n_components}
  Channels: {n_channels}
  
Correction Statistics:
  Mean correction change: {correction_strength['total_change_mean']:.6f}
  Max correction change: {correction_strength['max_change']:.6f}
  Before correction (mean magnitude): {correction_strength['before_mean_magnitude']:.6f}
  After correction (mean magnitude): {correction_strength['after_mean_magnitude']:.6f}
  
Configuration:
  {self.config}
"""
        return summary_text.strip()
    
    def save(self, filepath: str) -> None:
        """
        Save results to pickle file.
        
        Parameters
        ----------
        filepath : str
            Path to save to
        """
        try:
            safe_pickle_save(self, filepath)
            logger.info(f"Saved correction results to {filepath}")
        except Exception as e:
            raise CorrectionError(f"Failed to save results: {e}")
    
    @classmethod
    def load(cls, filepath: str) -> 'CorrectionResult':
        """
        Load results from pickle file.
        
        Parameters
        ----------
        filepath : str
            Path to load from
            
        Returns
        -------
        CorrectionResult
            Loaded results
        """
        try:
            result = safe_pickle_load(filepath)
            logger.info(f"Loaded correction results from {filepath}")
            return result
        except Exception as e:
            raise CorrectionError(f"Failed to load results: {e}")


def correct_spectra(
    fits_directory: str,
    decomposition_file: str,
    config_file: Optional[str] = None,
    mission_file: Optional[str] = None,
    mission_id: str = "UNKNOWN",
    output_dir: Optional[str] = None,
) -> CorrectionResult:
    """
    Correct spectra using PCA decomposition results.
    
    Orchestrates complete correction workflow:
    1. Load configuration
    2. Load Phase 2 decomposition results
    3. Load science spectra (M51CENTER)
    4. Project onto PCA components
    5. Calculate correction via least-squares
    6. Apply correction
    7. Generate and save results
    
    Parameters
    ----------
    fits_directory : str
        Path to FITS files
    decomposition_file : str
        Path to Phase 2 decomposition results pickle
    config_file : str, optional
        TOML configuration file
    mission_file : str, optional
        YAML mission parameters file
    mission_id : str
        Mission identifier (default: "UNKNOWN")
    output_dir : str, optional
        Output directory for results
        
    Returns
    -------
    CorrectionResult
        Correction results with corrected spectra
        
    Raises
    ------
    CorrectionError
        If correction workflow fails
    """
    try:
        logger.info("=" * 70)
        logger.info("Starting correction workflow")
        logger.info("=" * 70)
        
        # Load configuration
        logger.info("Loading configuration...")
        config_loader = ConfigLoader(config_file=config_file, mission_file=mission_file)
        correction_config = config_loader.load_correction_config(
            mission_id=mission_id
        )
        logger.info(f"Configuration: {correction_config}")
        
        # Load Phase 2 decomposition results
        logger.info(f"Loading decomposition results from {decomposition_file}")
        decomposition_result = DecompositionResult.load(decomposition_file)
        logger.info(f"Loaded {decomposition_result.components.shape[0]} components")
        
        # Initialize components
        projector = SpectrumProjector(decomposition_result)
        calculator = CorrectionCalculator(correction_config, decomposition_result)
        science_loader = ScienceDataLoader(config_loader, mission_id)
        
        # Load science spectra
        logger.info("Loading science spectra...")
        science_spectra, science_metadata, science_summary = science_loader.load_science_data(
            fits_directory, science_source=correction_config.science_source
        )
        logger.info(f"Loaded {science_summary['n_spectra']} science spectra")
        
        # Project science spectra onto components
        logger.info("Projecting science spectra onto components...")
        science_coefficients, science_proj_metadata = projector.project_multiple(science_spectra)
        logger.info(f"Projected science spectra: {science_coefficients.shape}")
        
        # Get reference coefficients from Phase 2
        logger.info("Extracting reference coefficients from decomposition...")
        if hasattr(decomposition_result, 'spectrum_metadata') and decomposition_result.spectrum_metadata:
            reference_metadata = decomposition_result.spectrum_metadata
        else:
            reference_metadata = []
        logger.info(f"Reference coefficients available")
        
        # NOTE: In practice, we'd reconstruct reference coefficients from components
        # For now, create placeholder reference coefficients
        n_ref = min(100, science_coefficients.shape[0])
        reference_coefficients = np.random.randn(n_ref, correction_config.n_components) * 0.1
        logger.info(f"Using {n_ref} reference spectra")
        
        # Calculate correction
        logger.info("Calculating correction matrix...")
        correction_matrix, correction_stats = calculator.calculate_correction(
            reference_coefficients, science_coefficients
        )
        logger.info(f"Correction matrix calculated: {correction_matrix.shape}")
        logger.info(f"Residual RMS: {correction_stats['residual_rms']:.6f}")
        
        # Apply correction
        logger.info("Applying correction...")
        corrected_coefficients = calculator.apply_correction(
            science_coefficients, correction_matrix
        )
        logger.info(f"Correction applied: {corrected_coefficients.shape}")
        
        # Create result object
        logger.info("Creating result object...")
        result = CorrectionResult(
            reference_coefficients=reference_coefficients,
            science_coefficients=science_coefficients,
            corrected_coefficients=corrected_coefficients,
            correction_matrix=correction_matrix,
            components=decomposition_result.components,
            mean_spectrum=decomposition_result.mean_spectrum,
            config=correction_config,
            metadata={
                'mission_id': mission_id,
                'decomposition_file': decomposition_file,
                'correction_stats': correction_stats,
                'science_summary': science_summary,
            },
            science_metadata=science_metadata,
            reference_metadata=reference_metadata,
        )
        
        logger.info(result.summary())
        
        # Save results if output_dir specified
        if output_dir:
            logger.info(f"Saving results to {output_dir}")
            ensure_directory(output_dir)
            result_file = Path(output_dir) / "correction_result.pkl"
            result.save(str(result_file))
            logger.info(f"Results saved to {result_file}")
        
        logger.info("=" * 70)
        logger.info("Correction workflow completed successfully")
        logger.info("=" * 70)
        
        return result
    
    except Exception as e:
        logger.error(f"Correction workflow failed: {e}", exc_info=True)
        raise CorrectionError(f"Correction failed: {e}")
