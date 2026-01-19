"""
PCA Decomposition module for spectral analysis.

Provides functionality to:
- Load reference spectra (SKYCHOPDIFF)
- Prepare spectra (baseline, mask, normalize)
- Perform PCA decomposition
- Store and manage decomposition results
"""

import numpy as np
import logging
from pathlib import Path
from typing import Tuple, List, Optional, Dict, Union
from dataclasses import dataclass, asdict
import pickle

from .config import ConfigLoader
from .fits_indexing import FITSIndexer, load_spectral_data
from .utilities import (
    fit_baseline_poly,
    create_channel_mask,
    combine_masks,
    normalize_spectrum,
    check_spectrum_validity,
    get_spectrum_stats,
    safe_pickle_save,
    safe_pickle_load,
    ensure_directory,
)
from .line_detection import (
    detect_science_line_waterfall,
    create_science_line_mask,
    create_artifact_mask,
    validate_line_detection,
)
from .errors import (
    DecompositionError,
    FITSError,
    LineDetectionError,
    BaselineError,
    InsufficientDataError,
)

logger = logging.getLogger(__name__)


@dataclass
class DecompositionConfig:
    """Configuration parameters for PCA decomposition."""
    
    n_components: int
    pca_source: str = "SKYCHOPDIFF"
    normalize: bool = True
    scale: bool = False
    add_sky_diff: bool = False
    noise_cutoff: bool = False
    scramble: bool = False
    
    def to_dict(self) -> dict:
        """Convert to dictionary."""
        return asdict(self)
    
    @classmethod
    def from_dict(cls, data: dict) -> 'DecompositionConfig':
        """Create from dictionary."""
        return cls(**data)


class SpectrumPreparator:
    """Prepare individual spectra for PCA decomposition.
    
    Performs per-spectrum processing:
    1. Detect science line from 2D data
    2. Create artifact and science line masks
    3. Fit baseline to masked regions
    4. Subtract baseline
    5. Normalize
    """
    
    def __init__(self, config: ConfigLoader, mission_id: str):
        """
        Initialize spectrum preparator.
        
        Parameters
        ----------
        config : ConfigLoader
            Configuration loader with mission parameters
        mission_id : str
            Mission ID for accessing mission-specific parameters
        """
        self.config = config
        self.mission_id = mission_id
        
        # Get mission-specific line parameters
        try:
            line_params = config.get_line_parameters(mission_id)
            self.artifact_center = line_params.get('telluric_line_center', 50)
            self.artifact_width = line_params.get('telluric_line_width', 10)
        except Exception:
            logger.warning(f"Could not load line parameters for {mission_id}, using defaults")
            self.artifact_center = 50
            self.artifact_width = 10
        
        # Get baseline fitting parameters
        self.baseline_order = config.get('pca.common.baseline_order', default=3)
        self.smoothing_kernel = config.get('pca.common.smoothing_kernel_size', default=3)
        self.rolling_window = config.get('pca.common.rolling_noise_window', default=11)
    
    def prepare(self, spectrum_1d: np.ndarray, spectrum_2d: np.ndarray) -> Tuple[np.ndarray, dict]:
        """
        Prepare a single spectrum for decomposition.
        
        Parameters
        ----------
        spectrum_1d : np.ndarray
            1D spectrum array [n_channels]
        spectrum_2d : np.ndarray
            2D spectral data [n_readouts, n_channels] for line detection
            
        Returns
        -------
        prepared_spectrum : np.ndarray
            Prepared spectrum (baseline subtracted, masked, normalized)
        metadata : dict
            Processing metadata (line position, baseline stats, etc.)
            
        Raises
        ------
        DecompositionError
            If spectrum preparation fails
        BaselineError
            If baseline fitting fails
        LineDetectionError
            If line detection fails
        """
        try:
            metadata = {}
            spectrum = spectrum_1d.copy().astype(float)
            n_channels = len(spectrum)
            
            # Validate input spectrum
            if not check_spectrum_validity(spectrum, min_valid=n_channels // 2):
                raise DecompositionError(f"Spectrum has insufficient valid data ({n_channels})")
            
            metadata['input_stats'] = get_spectrum_stats(spectrum)
            
            # Step 1: Detect science line from 2D data
            try:
                line_center, line_confidence = detect_science_line_waterfall(
                    spectrum_2d,
                    kernel_size=self.smoothing_kernel
                )
                metadata['line_center'] = line_center
                metadata['line_confidence'] = float(line_confidence)
                
                # Validate detection
                if not validate_line_detection(line_center, n_channels, width=20):
                    logger.warning(f"Line detection at {line_center} failed validation")
                    line_center = np.argmax(np.nanmean(spectrum_2d, axis=0))
                    metadata['line_center'] = line_center
                    metadata['line_confidence'] = 0.0
            
            except LineDetectionError as e:
                logger.warning(f"Line detection failed: {e}, using peak position")
                line_center = int(np.nanargmax(spectrum))
                metadata['line_center'] = line_center
                metadata['line_confidence'] = 0.0
            
            # Step 2: Create masks
            # Artifact mask (from config - fixed hardware issue)
            artifact_mask = create_artifact_mask(
                n_channels,
                artifact_center=self.artifact_center,
                artifact_width=self.artifact_width
            )
            
            # Science line mask (auto-detected)
            line_mask = create_science_line_mask(spectrum, center=line_center, width=20)
            
            # Combined mask (exclude both artifact and science line during baseline fitting)
            combined_mask = combine_masks(artifact_mask, line_mask)
            metadata['masked_channels'] = int(np.sum(combined_mask))
            
            # Step 3: Fit baseline to unmasked regions
            valid_mask = ~combined_mask
            if np.sum(valid_mask) < self.baseline_order + 10:
                raise BaselineError(
                    f"Not enough valid channels for baseline fitting: {np.sum(valid_mask)}"
                )
            
            try:
                baseline = fit_baseline_poly(spectrum, order=self.baseline_order, mask=valid_mask)
                metadata['baseline_order'] = self.baseline_order
                metadata['baseline_stats'] = get_spectrum_stats(baseline)
            except BaselineError as e:
                logger.warning(f"Baseline fitting failed: {e}, using zero baseline")
                baseline = np.zeros_like(spectrum)
                metadata['baseline_order'] = -1
            
            # Step 4: Subtract baseline
            spectrum = spectrum - baseline
            
            # Step 5: Apply mask (set masked regions to zero)
            spectrum[combined_mask] = 0.0
            
            # Step 6: Normalize
            if np.std(spectrum) > 0:
                normalized, mean, std = normalize_spectrum(spectrum)
                metadata['normalization_mean'] = float(mean)
                metadata['normalization_std'] = float(std)
            else:
                normalized = spectrum.copy()
                metadata['normalization_mean'] = 0.0
                metadata['normalization_std'] = 1.0
            
            metadata['output_stats'] = get_spectrum_stats(normalized)
            metadata['processing_status'] = 'success'
            
            return normalized, metadata
        
        except Exception as e:
            raise DecompositionError(f"Spectrum preparation failed: {e}")


class DecompositionDataLoader:
    """Load and organize FITS data for PCA decomposition.
    
    Handles:
    - FITS file discovery and filtering
    - Data quality filtering
    - Spectrum extraction and organization
    """
    
    def __init__(self, fits_directory: Union[str, Path], config: ConfigLoader):
        """
        Initialize data loader.
        
        Parameters
        ----------
        fits_directory : str or Path
            Directory containing FITS files
        config : ConfigLoader
            Configuration loader
        """
        self.fits_directory = Path(fits_directory)
        self.config = config
        self.indexer = None
        self.mission_id = None
    
    def load_decomposition_data(self, mission_id: str) -> Tuple[np.ndarray, List[dict], dict]:
        """
        Load and prepare decomposition data.
        
        Loads all SKYCHOPDIFF spectra, applies quality filters, and prepares them.
        
        Parameters
        ----------
        mission_id : str
            Mission ID for accessing mission-specific parameters
            
        Returns
        -------
        spectral_matrix : np.ndarray
            2D array [n_spectra, n_channels] of prepared spectra
        metadata_list : list
            List of metadata dictionaries for each spectrum
        summary : dict
            Summary statistics of loaded data
            
        Raises
        ------
        FITSError
            If FITS files can't be read
        InsufficientDataError
            If not enough spectra are available
        """
        try:
            self.mission_id = mission_id
            
            # Index FITS files
            logger.info(f"Scanning FITS directory: {self.fits_directory}")
            self.indexer = FITSIndexer(str(self.fits_directory))
            self.indexer.scan_directory()
            logger.info(f"Found {len(self.indexer.index_df)} FITS files")
            
            # Filter by source (SKYCHOPDIFF)
            pca_source = self.config.get('pca.common.pca_source', default='SKYCHOPDIFF')
            logger.info(f"Filtering for source: {pca_source}")
            source_data = self.indexer.filter_by_source(pca_source)
            logger.info(f"Found {len(source_data.index_df)} {pca_source} spectra")
            
            # Apply drop filters (exclude bad data)
            try:
                drop_filters = self.config.get_drop_filters(mission_id)
                filtered_df = source_data.index_df.copy()
                
                if drop_filters.get('exclude_telescopes'):
                    for telescope in drop_filters['exclude_telescopes']:
                        filtered_df = filtered_df[filtered_df['telescope'] != telescope]
                
                if drop_filters.get('exclude_scans'):
                    for scan in drop_filters['exclude_scans']:
                        filtered_df = filtered_df[filtered_df['scan'] != scan]
                
                logger.info(f"After filtering: {len(filtered_df)} spectra")
                
                if len(filtered_df) == 0:
                    raise InsufficientDataError("No spectra remain after applying drop filters")
                
                source_data.index_df = filtered_df.reset_index(drop=True)
            
            except Exception as e:
                logger.warning(f"Could not apply drop filters: {e}")
            
            # Load and prepare spectra
            preparator = SpectrumPreparator(self.config, mission_id)
            spectral_list = []
            metadata_list = []
            
            logger.info("Loading and preparing spectra...")
            
            for data, header, metadata in source_data.get_all_fits_data():
                try:
                    # Extract 1D and 2D spectra
                    if data.ndim == 3:
                        # [n_readouts, n_subscan, n_channels]
                        spectrum_2d = data[:, 0, :]  # Use first subscan
                        spectrum_1d = np.nanmean(spectrum_2d, axis=0)
                    elif data.ndim == 2:
                        # [n_readouts, n_channels]
                        spectrum_2d = data
                        spectrum_1d = np.nanmean(spectrum_2d, axis=0)
                    else:
                        # [n_channels]
                        spectrum_1d = data
                        spectrum_2d = spectrum_1d.reshape(1, -1)
                    
                    # Prepare spectrum
                    prepared, prep_metadata = preparator.prepare(spectrum_1d, spectrum_2d)
                    
                    spectral_list.append(prepared)
                    
                    # Combine metadata
                    combined_metadata = {**metadata, **prep_metadata}
                    metadata_list.append(combined_metadata)
                
                except Exception as e:
                    logger.warning(f"Failed to prepare spectrum: {e}")
                    continue
            
            if not spectral_list:
                raise InsufficientDataError("Failed to prepare any valid spectra")
            
            logger.info(f"Successfully prepared {len(spectral_list)} spectra")
            
            # Create spectral matrix
            spectral_matrix = np.array(spectral_list)
            
            # Summary statistics
            summary = {
                'n_spectra': len(spectral_list),
                'n_channels': spectral_matrix.shape[1],
                'mission_id': mission_id,
                'pca_source': pca_source,
                'mean_spectrum': np.mean(spectral_matrix, axis=0),
                'std_spectrum': np.std(spectral_matrix, axis=0),
            }
            
            logger.info(f"Decomposition data shape: {spectral_matrix.shape}")
            
            return spectral_matrix, metadata_list, summary
        
        except Exception as e:
            raise FITSError(f"Failed to load decomposition data: {e}")


class PCADecomposer:
    """Perform PCA decomposition on spectral data.
    
    Uses sklearn's PCA for eigenvalue decomposition.
    """
    
    def __init__(self, n_components: int, scale: bool = False):
        """
        Initialize PCA decomposer.
        
        Parameters
        ----------
        n_components : int
            Number of PCA components to extract
        scale : bool
            Whether to standardize features (default: False)
        """
        self.n_components = n_components
        self.scale = scale
        self.pca_model = None
        self.fitted = False
    
    def fit(self, spectral_matrix: np.ndarray) -> 'PCADecomposer':
        """
        Fit PCA to spectral data.
        
        Parameters
        ----------
        spectral_matrix : np.ndarray
            2D array [n_spectra, n_channels] of spectral data
            
        Returns
        -------
        self
            Returns self for method chaining
            
        Raises
        ------
        DecompositionError
            If PCA fitting fails
        """
        try:
            from sklearn.decomposition import PCA
            
            logger.info(f"Fitting PCA with {self.n_components} components")
            
            # Validate input
            if spectral_matrix.shape[0] < self.n_components:
                raise DecompositionError(
                    f"Need at least {self.n_components} spectra, "
                    f"but got {spectral_matrix.shape[0]}"
                )
            
            # Fit PCA
            self.pca_model = PCA(
                n_components=self.n_components,
                whiten=self.scale,
                random_state=42
            )
            
            self.pca_model.fit(spectral_matrix)
            self.fitted = True
            
            logger.info(f"PCA fitting complete")
            logger.info(f"Explained variance ratio: {self.pca_model.explained_variance_ratio_}")
            logger.info(f"Total variance explained: {np.sum(self.pca_model.explained_variance_ratio_):.1%}")
            
            return self
        
        except ImportError:
            raise DecompositionError("scikit-learn is required. Install with: pip install scikit-learn")
        except Exception as e:
            raise DecompositionError(f"PCA fitting failed: {e}")
    
    def get_components(self) -> np.ndarray:
        """
        Get PCA components.
        
        Returns
        -------
        np.ndarray
            Components array [n_components, n_channels]
        """
        if not self.fitted or self.pca_model is None:
            raise DecompositionError("PCA model not fitted. Call fit() first.")
        
        return self.pca_model.components_
    
    def get_explained_variance_ratio(self) -> np.ndarray:
        """
        Get explained variance ratio per component.
        
        Returns
        -------
        np.ndarray
            Variance ratio [n_components]
        """
        if not self.fitted or self.pca_model is None:
            raise DecompositionError("PCA model not fitted. Call fit() first.")
        
        return self.pca_model.explained_variance_ratio_
    
    def get_explained_variance(self) -> np.ndarray:
        """
        Get explained variance per component.
        
        Returns
        -------
        np.ndarray
            Variance [n_components]
        """
        if not self.fitted or self.pca_model is None:
            raise DecompositionError("PCA model not fitted. Call fit() first.")
        
        return self.pca_model.explained_variance_
    
    def transform(self, spectral_matrix: np.ndarray) -> np.ndarray:
        """
        Project data onto PCA components.
        
        Parameters
        ----------
        spectral_matrix : np.ndarray
            Input data [n_spectra, n_channels]
            
        Returns
        -------
        np.ndarray
            Coefficients [n_spectra, n_components]
        """
        if not self.fitted or self.pca_model is None:
            raise DecompositionError("PCA model not fitted. Call fit() first.")
        
        return self.pca_model.transform(spectral_matrix)
    
    def reconstruct(self, coefficients: np.ndarray) -> np.ndarray:
        """
        Reconstruct from PCA coefficients.
        
        Parameters
        ----------
        coefficients : np.ndarray
            PCA coefficients [n_spectra, n_components]
            
        Returns
        -------
        np.ndarray
            Reconstructed spectra [n_spectra, n_channels]
        """
        if not self.fitted or self.pca_model is None:
            raise DecompositionError("PCA model not fitted. Call fit() first.")
        
        return self.pca_model.inverse_transform(coefficients)


@dataclass
class DecompositionResult:
    """Results from PCA decomposition."""
    
    components: np.ndarray  # [n_components, n_channels]
    explained_variance_ratio: np.ndarray  # [n_components]
    explained_variance: np.ndarray  # [n_components]
    mean_spectrum: np.ndarray  # [n_channels] - mean of training data
    config: DecompositionConfig
    metadata: dict  # Mission, date, tag, etc.
    spectrum_metadata: List[dict] = None  # Per-spectrum metadata
    
    def __post_init__(self):
        """Validate after initialization."""
        if self.spectrum_metadata is None:
            self.spectrum_metadata = []
    
    def summary(self) -> str:
        """Get summary of decomposition results."""
        lines = [
            "Decomposition Results Summary:",
            f"  Components: {self.components.shape}",
            f"  Total variance explained: {np.sum(self.explained_variance_ratio):.1%}",
            f"  Variance per component: {self.explained_variance_ratio}",
            f"  Configuration: {self.config}",
        ]
        return "\n".join(lines)
    
    def save(self, filepath: Union[str, Path]) -> None:
        """
        Save decomposition results to pickle file.
        
        Parameters
        ----------
        filepath : str or Path
            Output pickle file path
        """
        filepath = Path(filepath)
        logger.info(f"Saving decomposition results to {filepath}")
        
        safe_pickle_save(asdict(self), filepath)
        logger.info(f"Successfully saved to {filepath}")
    
    @classmethod
    def load(cls, filepath: Union[str, Path]) -> 'DecompositionResult':
        """
        Load decomposition results from pickle file.
        
        Parameters
        ----------
        filepath : str or Path
            Pickle file path
            
        Returns
        -------
        DecompositionResult
            Loaded results
        """
        filepath = Path(filepath)
        logger.info(f"Loading decomposition results from {filepath}")
        
        data = safe_pickle_load(filepath)
        
        # Reconstruct config from dict
        config_dict = data['config']
        config = DecompositionConfig.from_dict(config_dict)
        
        result = cls(
            components=data['components'],
            explained_variance_ratio=data['explained_variance_ratio'],
            explained_variance=data['explained_variance'],
            mean_spectrum=data['mean_spectrum'],
            config=config,
            metadata=data['metadata'],
            spectrum_metadata=data.get('spectrum_metadata', [])
        )
        
        logger.info(f"Successfully loaded from {filepath}")
        return result


def decompose_spectra(fits_directory: Union[str, Path],
                     config_file: Union[str, Path],
                     mission_file: Union[str, Path],
                     mission_id: str,
                     output_dir: Optional[Union[str, Path]] = None) -> DecompositionResult:
    """
    Perform complete PCA decomposition on FITS spectra.
    
    Main entry point for Phase 2. Orchestrates the full decomposition workflow:
    1. Load configuration
    2. Index FITS files
    3. Load SKYCHOPDIFF spectra
    4. Prepare spectra (baseline, mask, normalize)
    5. Fit PCA
    6. Return results
    
    Parameters
    ----------
    fits_directory : str or Path
        Directory containing FITS files
    config_file : str or Path
        Path to TOML configuration file
    mission_file : str or Path
        Path to YAML mission parameters file
    mission_id : str
        Mission ID (e.g., 'SOFIA_MISSION')
    output_dir : str or Path, optional
        Directory to save results (default: None = don't save)
        
    Returns
    -------
    DecompositionResult
        PCA decomposition results with components and variance
        
    Raises
    ------
    DecompositionError
        If decomposition fails
    FITSError
        If FITS files can't be read
    ConfigurationError
        If configuration is invalid
    """
    try:
        logger.info("Starting PCA decomposition workflow...")
        
        # Step 1: Load configuration
        logger.info("Loading configuration...")
        config = ConfigLoader(config_file=str(config_file), 
                             mission_file=str(mission_file))
        
        # Step 2: Load decomposition parameters
        n_components = config.get('pca.decompose.number_components')
        scale = config.get('pca.common.do_scale', default=False)
        tag = config.get('pca.common.tag', default='decomposition')
        
        decomp_config = DecompositionConfig(
            n_components=n_components,
            scale=scale,
            add_sky_diff=config.get('pca.decompose.add_sky_diff', default=False),
            noise_cutoff=config.get('pca.decompose.noise_cutoff', default=False),
        )
        
        logger.info(f"Configuration loaded: {n_components} components, scale={scale}")
        
        # Step 3: Load data
        logger.info("Loading spectral data...")
        data_loader = DecompositionDataLoader(fits_directory, config)
        spectral_matrix, spectrum_metadata, summary = data_loader.load_decomposition_data(mission_id)
        
        logger.info(f"Loaded spectral data: shape {spectral_matrix.shape}")
        logger.info(f"  Mean spectrum: {np.mean(spectral_matrix, axis=0)[:10]}...")
        
        # Step 4: Fit PCA
        logger.info("Fitting PCA decomposition...")
        decomposer = PCADecomposer(n_components=n_components, scale=scale)
        decomposer.fit(spectral_matrix)
        
        # Step 5: Create result object
        logger.info("Creating decomposition result...")
        result = DecompositionResult(
            components=decomposer.get_components(),
            explained_variance_ratio=decomposer.get_explained_variance_ratio(),
            explained_variance=decomposer.get_explained_variance(),
            mean_spectrum=summary['mean_spectrum'],
            config=decomp_config,
            metadata={
                'mission_id': mission_id,
                'tag': tag,
                'n_spectra': summary['n_spectra'],
                'n_channels': summary['n_channels'],
                'pca_source': summary['pca_source'],
                'total_variance_explained': float(np.sum(decomposer.get_explained_variance_ratio())),
            },
            spectrum_metadata=spectrum_metadata
        )
        
        logger.info(result.summary())
        
        # Step 6: Save results if output directory specified
        if output_dir:
            output_dir = ensure_directory(output_dir)
            output_file = output_dir / f"decomposition_{tag}.pkl"
            result.save(output_file)
            logger.info(f"Results saved to {output_file}")
        
        logger.info("PCA decomposition workflow complete!")
        return result
    
    except Exception as e:
        logger.error(f"Decomposition workflow failed: {e}")
        raise DecompositionError(f"Decomposition failed: {e}")
