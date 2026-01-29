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
        
        # Get mission-specific line parameters (in km/s from YAML)
        try:
            line_params = config.get_line_parameters(mission_id)
            self.telluric_center = line_params.get('center', 50)
            self.telluric_width = line_params.get('width', 10)
            logger.info(
                f"Loaded telluric parameters for {mission_id}: "
                f"center={self.telluric_center} km/s, width={self.telluric_width} km/s"
            )
        except Exception:
            logger.warning(f"Could not load line parameters for {mission_id}, using defaults")
            self.telluric_center = 50
            self.telluric_width = 10
        
        # Get baseline fitting parameters
        self.baseline_order = config.get('pca.common.baseline_order', default=3)
        self.smoothing_kernel = config.get('pca.common.smoothing_kernel_size', default=3)
        self.rolling_window = config.get('pca.common.rolling_noise_window', default=11)
        
        # Store for header-based WCS conversion and velocity axis
        self.header = None
        self.velocity_axis = None
    
    def prepare(self, spectrum_1d: np.ndarray, spectrum_2d: np.ndarray,
                header: Optional[dict] = None, velocity_axis: Optional[np.ndarray] = None) -> Tuple[np.ndarray, dict]:
        """
        Prepare a single spectrum for decomposition.
        
        Parameters
        ----------
        spectrum_1d : np.ndarray
            1D spectrum array [n_channels]
        spectrum_2d : np.ndarray
            2D spectral data [n_readouts, n_channels] for line detection
        header : dict, optional
            FITS header with WCS information (CRVAL1, CDELT1, CRPIX1)
            Used for velocity-to-channel conversion for telluric masking
        velocity_axis : np.ndarray, optional
            1D array of velocity values (from VELOCITY_AXIS column in FITS file)
            If provided, this is used preferentially over WCS parameters
            
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
            
            # Store header and velocity axis for telluric masking
            if header is not None:
                self.header = header
            if velocity_axis is not None:
                self.velocity_axis = np.asarray(velocity_axis).flatten()
                metadata['has_velocity_axis'] = True
            else:
                self.velocity_axis = None
                metadata['has_velocity_axis'] = False
            
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
            # Artifact mask (telluric lines - velocity-based with WCS or VELOCITY_AXIS)
            artifact_mask = create_artifact_mask(
                n_channels,
                artifact_center=self.telluric_center,
                artifact_width=self.telluric_width,
                header=self.header,
                velocity_axis=self.velocity_axis,
                use_velocity=True
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
                    # Extract VELOCITY_AXIS column if available (preferred for telluric masking)
                    velocity_axis = None
                    
                    # Check if data is a structured array with columns (from FITS binary table)
                    if isinstance(data, np.ndarray) and data.dtype.names is not None:
                        if 'VELOCITY_AXIS' in data.dtype.names:
                            # Extract first row's velocity axis (should be same for all)
                            velocity_axis = data['VELOCITY_AXIS'][0] if len(data) > 0 else None
                            logger.debug("Using VELOCITY_AXIS from FITS file")
                    
                    # Extract 1D and 2D spectra from SPECTRUM column (if binary table)
                    if isinstance(data, np.ndarray) and data.dtype.names is not None and 'SPECTRUM' in data.dtype.names:
                        # Binary table with SPECTRUM column
                        # Use first spectrum for decomposition (all should be similar)
                        spectrum_1d = data['SPECTRUM'][0].flatten() if len(data) > 0 else np.array([])
                        spectrum_2d = data['SPECTRUM'].reshape(len(data), -1) if len(data) > 0 else spectrum_1d.reshape(1, -1)
                    else:
                        # Regular array data (e.g., from single HDU)
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
                    
                    # Prepare spectrum (pass header and velocity_axis for telluric masking)
                    prepared, prep_metadata = preparator.prepare(
                        spectrum_1d, spectrum_2d, 
                        header=header, 
                        velocity_axis=velocity_axis
                    )
                    
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
    reference_spectrum_indices: np.ndarray = None  # [n_reference_spectra] - indices of SKYCHOPDIFF spectra in original FITS file
    
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
        
        # Explicitly include all fields to ensure reference_spectrum_indices is saved
        save_dict = {
            'components': self.components,
            'explained_variance_ratio': self.explained_variance_ratio,
            'explained_variance': self.explained_variance,
            'mean_spectrum': self.mean_spectrum,
            'config': asdict(self.config),
            'metadata': self.metadata,
            'spectrum_metadata': self.spectrum_metadata,
            'reference_spectrum_indices': self.reference_spectrum_indices
        }
        
        if self.reference_spectrum_indices is not None:
            logger.info(f"  Including {len(self.reference_spectrum_indices)} spectrum indices in save")
        else:
            logger.info(f"  reference_spectrum_indices is None - will not include in save")
        
        safe_pickle_save(save_dict, filepath)
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
            spectrum_metadata=data.get('spectrum_metadata', []),
            reference_spectrum_indices=data.get('reference_spectrum_indices', None)
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


def _plot_decomposition_results(result: DecompositionResult, velocity_axis: np.ndarray, 
                                 mission_id: str, flight_date: str, output_dir: Path) -> Tuple[Path, Path]:
    """
    Create visualization plots of PCA components.
    
    Parameters
    ----------
    result : DecompositionResult
        Decomposition result object
    velocity_axis : np.ndarray
        Velocity axis in km/s for x-axis labels
    mission_id : str
        Mission ID for plot filenames
    flight_date : str
        Flight date for plot filenames
    output_dir : Path
        Output directory for plots
        
    Returns
    -------
    tuple of Path
        Paths to the two created plot files
    """
    try:
        import matplotlib.pyplot as plt
        import matplotlib
        # Use non-interactive backend
        matplotlib.use('Agg')
    except ImportError:
        logger.warning("matplotlib not available, skipping plots")
        return None, None
    
    components = result.components
    explained_variance = result.explained_variance_ratio
    mean_spectrum = result.mean_spectrum
    
    n_components = components.shape[0]
    
    # Create component plots with velocity axis
    fig, axes = plt.subplots(n_components + 1, 1, figsize=(14, 3*(n_components + 1)))
    
    # Plot mean spectrum with velocity axis
    axes[0].plot(velocity_axis, mean_spectrum, 'b-', linewidth=1.5)
    axes[0].axvspan(592-15, 592+15, alpha=0.3, color='red', label='Telluric mask (592±15 km/s)')
    axes[0].set_title('Mean Spectrum (All Spectra) - Telluric Region Masked', fontsize=12, fontweight='bold')
    axes[0].set_ylabel('Intensity')
    axes[0].set_xlim(velocity_axis.min(), velocity_axis.max())
    axes[0].grid(True, alpha=0.3)
    axes[0].legend()
    
    # Plot each component
    for i in range(n_components):
        axes[i+1].plot(velocity_axis, components[i], 'r-', linewidth=1.5)
        axes[i+1].axvspan(592-15, 592+15, alpha=0.3, color='red')
        axes[i+1].set_title(
            f'Component {i+1} (Variance: {100*explained_variance[i]:.2f}%)',
            fontsize=12, fontweight='bold'
        )
        axes[i+1].set_ylabel('Loadings')
        axes[i+1].set_xlim(velocity_axis.min(), velocity_axis.max())
        axes[i+1].grid(True, alpha=0.3)
        if i == n_components - 1:
            axes[i+1].set_xlabel('Velocity (km/s)')
    
    plt.tight_layout()
    components_path = output_dir / f"pca_components_{mission_id}_{flight_date}.png"
    plt.savefig(components_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    logger.info(f"  Saved component plots to {components_path}")
    
    # Create variance explained plot
    fig2, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
    
    # Individual variance
    ax1.bar(range(1, n_components+1), 100*explained_variance, color='steelblue', alpha=0.8, edgecolor='black')
    ax1.set_xlabel('Component', fontsize=12, fontweight='bold')
    ax1.set_ylabel('Variance Explained (%)', fontsize=12, fontweight='bold')
    ax1.set_title('Individual Component Variance (With Telluric Masking)', fontsize=12, fontweight='bold')
    ax1.grid(True, alpha=0.3, axis='y')
    for i, v in enumerate(explained_variance):
        ax1.text(i+1, 100*v + 0.1, f'{100*v:.2f}%', ha='center', va='bottom', fontsize=10, fontweight='bold')
    
    # Cumulative variance
    cumsum_variance = np.cumsum(explained_variance)
    ax2.plot(range(1, n_components+1), 100*cumsum_variance, 'o-', linewidth=2.5, markersize=8, color='darkred')
    ax2.fill_between(range(1, n_components+1), 0, 100*cumsum_variance, alpha=0.3, color='darkred')
    ax2.set_xlabel('Number of Components', fontsize=12, fontweight='bold')
    ax2.set_ylabel('Cumulative Variance Explained (%)', fontsize=12, fontweight='bold')
    ax2.set_title('Cumulative Variance Explained (With Telluric Masking)', fontsize=12, fontweight='bold')
    ax2.grid(True, alpha=0.3)
    ax2.set_ylim([0, np.max(100*cumsum_variance) * 1.1])
    for i, v in enumerate(cumsum_variance):
        ax2.text(i+1, 100*v + 0.3, f'{100*v:.2f}%', ha='center', va='bottom', fontsize=10, fontweight='bold')
    
    plt.tight_layout()
    variance_path = output_dir / f"pca_variance_{mission_id}_{flight_date}.png"
    plt.savefig(variance_path, dpi=150, bbox_inches='tight')
    plt.close(fig2)
    logger.info(f"  Saved variance plot to {variance_path}")
    
    return components_path, variance_path


def main_cli():
    """Command-line interface for PCA decomposition."""
    import sys
    import argparse
    from astropy.io import fits as fits_io
    
    parser = argparse.ArgumentParser(
        description="PCA decomposition of spectral reference data",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Example:\n  pca_decompose --config config.toml\n"
               "  pca_decompose --config config.toml --n-components 10 --pca-source M51\n"
               "  pca_decompose --config config.toml --fits /path/to/custom.fits"
    )
    parser.add_argument(
        "--config",
        type=str,
        default="config.toml",
        help="Path to configuration file (default: config.toml)"
    )
    parser.add_argument(
        "--fits",
        type=str,
        default=None,
        help="FITS file to analyze (overrides config file [output][reduced_fits])"
    )
    parser.add_argument(
        "--n-components",
        type=int,
        default=None,
        help="Number of PCA components to extract (overrides config file, default: 5)"
    )
    parser.add_argument(
        "--pca-source",
        type=str,
        default=None,
        help="PCA source object name (overrides config file, default: SKYCHOPDIFF)"
    )
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Verbose output"
    )
    parser.add_argument(
        "--plot-components",
        action="store_true",
        help="Generate visualization plots of PCA components with velocity axis"
    )
    
    args = parser.parse_args()
    
    # Set logging level
    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)
    
    try:
        # Load config
        logger.info(f"Loading configuration from {args.config}")
        config = ConfigLoader(args.config)
        
        # Load mission-specific parameters from mission_id_parameters.yml
        mission_yml = Path(__file__).parent / "mission_id_parameters.yml"
        if mission_yml.exists():
            logger.info(f"Loading mission parameters from {mission_yml}")
            config.load_missions(mission_yml)
            logger.debug(f"Loaded mission parameters. Available missions: {len(config.mission_params)}")
        else:
            logger.warning(f"Mission parameters file not found: {mission_yml}")
        
        # Determine FITS file (command line > config [output][reduced_fits] > config [input][fits_file])
        if args.fits is not None:
            fits_file = args.fits
            logger.info(f"✓ FITS file = {fits_file} (from command line)")
        else:
            # Try [output][reduced_fits] first
            output_config = config.get('output', {})
            fits_file = output_config.get('reduced_fits')
            
            if fits_file:
                logger.info(f"✓ FITS file = {fits_file} (from config [output][reduced_fits])")
            else:
                # Fall back to [input][fits_file]
                input_config = config.get('input', {})
                fits_file = input_config.get('fits_file')
                if fits_file:
                    logger.info(f"✓ FITS file = {fits_file} (from config [input][fits_file])")
                else:
                    logger.error("No FITS file specified in configuration or command line")
                    logger.error("Specify via:")
                    logger.error("  1. Command line: pca_decompose --fits /path/to/file.fits")
                    logger.error("  2. Config [output][reduced_fits]")
                    logger.error("  3. Config [input][fits_file]")
                    sys.exit(1)
        
        # Verify file exists
        if not Path(fits_file).exists():
            logger.error(f"FITS file not found: {fits_file}")
            sys.exit(1)
        
        logger.info(f"Analyzing file: {fits_file}")

        
        # Determine n_components (command line > config file > default)
        if args.n_components is not None:
            n_components = args.n_components
            logger.info(f"✓ n_components = {n_components} (from command line)")
        else:
            pca_config = config.get('pca', {})
            n_components = pca_config.get('n_components', 5)
            if 'n_components' in pca_config:
                logger.info(f"✓ n_components = {n_components} (from config file)")
            else:
                logger.info(f"✓ n_components = {n_components} (default)")
        
        # Determine pca_source (command line > config file > default)
        if args.pca_source is not None:
            pca_source = args.pca_source
            logger.info(f"✓ pca_source = {pca_source} (from command line)")
        else:
            pca_config = config.get('pca', {})
            pca_source = pca_config.get('pca_source', 'SKYCHOPDIFF')
            if 'pca_source' in pca_config:
                logger.info(f"✓ pca_source = {pca_source} (from config file)")
            else:
                logger.info(f"✓ pca_source = {pca_source} (default)")
        
        # Load spectra and extract mission_id from data
        logger.info(f"Loading {pca_source} spectra from {fits_file}")
        with fits_io.open(fits_file) as hdul:
            # Determine which HDU contains the matrix data
            matrix_hdu = None
            if 'MATRIX' in hdul:
                matrix_hdu = hdul['MATRIX']
            elif len(hdul) > 1:
                # Try the first binary table HDU after primary
                for hdu in hdul[1:]:
                    if isinstance(hdu, fits_io.BinTableHDU):
                        matrix_hdu = hdu
                        break
            
            if matrix_hdu is None:
                logger.error(f"Could not find data HDU in {fits_file}")
                sys.exit(1)
            
            # Extract spectra for the requested source
            # Handle both string and bytes comparison for OBJECT column
            if isinstance(matrix_hdu.data['OBJECT'][0], bytes):
                mask = matrix_hdu.data['OBJECT'] == pca_source.encode() if isinstance(pca_source, str) else matrix_hdu.data['OBJECT'] == pca_source
            else:
                mask = matrix_hdu.data['OBJECT'] == pca_source
            spectra_raw = np.array([row['SPECTRUM'] for row in matrix_hdu.data[mask]])
            
            if len(spectra_raw) == 0:
                logger.error(f"No spectra found with OBJECT = '{pca_source}'")
                logger.info("Available objects:")
                unique_objects = np.unique(matrix_hdu.data['OBJECT'])
                for obj in unique_objects:
                    obj_str = obj.decode('utf-8') if isinstance(obj, bytes) else str(obj)
                    count = np.sum(matrix_hdu.data['OBJECT'] == obj)
                    logger.info(f"  {obj_str}: {count} spectra")
                sys.exit(1)
            
            # Get flight date and mission_id for filename
            skychopdiff_data = matrix_hdu.data[mask]
            if len(skychopdiff_data) > 0:
                date_obs = skychopdiff_data[0]['DATE-OBS']
                if isinstance(date_obs, bytes):
                    date_obs = date_obs.decode().strip()
                flight_date = date_obs.split('T')[0].replace('-', '')
                
                # Get mission_id from MISSION_ID column if it exists
                if 'MISSION_ID' in matrix_hdu.data.dtype.names:
                    mission_id = skychopdiff_data[0]['MISSION_ID']
                    if isinstance(mission_id, bytes):
                        mission_id = mission_id.decode().strip()
                else:
                    mission_id = 'unknown'
            else:
                flight_date = 'unknown'
                mission_id = 'unknown'
            
            # Extract velocity axis if available (for telluric masking)
            velocity_axis = None
            if 'VELOCITY_AXIS' in matrix_hdu.data.dtype.names:
                velocity_axis = matrix_hdu.data['VELOCITY_AXIS'][0] if len(matrix_hdu.data) > 0 else None
                logger.debug("Using VELOCITY_AXIS from FITS file for telluric masking")
            
            # Get header for WCS information
            matrix_header = matrix_hdu.header
        
        logger.info(f"  Loaded {len(spectra_raw)} {pca_source} spectra")
        if mission_id != 'unknown':
            logger.info(f"  Mission ID: {mission_id}")
        
        # Apply spectrum preparation with telluric masking
        logger.debug(f"Creating SpectrumPreparator for mission {mission_id}")
        preparator = SpectrumPreparator(config, mission_id)
        spectra = []
        failed_count = 0
        
        for i, spectrum_row in enumerate(spectra_raw):
            try:
                # For single spectra, use as both 1D and 2D
                spectrum_1d = spectrum_row.flatten()
                spectrum_2d = spectrum_row.reshape(1, -1)
                
                # Apply preparation (includes telluric masking)
                prepared, _ = preparator.prepare(
                    spectrum_1d, spectrum_2d,
                    header=matrix_header,
                    velocity_axis=velocity_axis
                )
                
                # Ensure no NaNs in prepared spectrum
                if np.any(np.isnan(prepared)):
                    logger.warning(f"Spectrum {i} contains NaN after preparation, filling with zeros")
                    prepared = np.nan_to_num(prepared, nan=0.0)
                
                spectra.append(prepared)
            except Exception as e:
                logger.warning(f"Could not prepare spectrum {i}: {e}")
                failed_count += 1
                # Don't fall back to raw spectrum - skip it instead
                continue
        
        if failed_count > 0:
            logger.warning(f"Failed to prepare {failed_count} spectra")
        
        spectra = np.array(spectra)
        logger.info(f"Prepared spectra: {spectra.shape[0]} spectra × {spectra.shape[1]} channels")
        
        logger.info(f"Prepared spectra: {spectra.shape[0]} spectra × {spectra.shape[1]} channels")
        
        # Perform PCA decomposition
        logger.info(f"Performing PCA decomposition ({n_components} components)...")
        decomposer = PCADecomposer(
            n_components=n_components,
            scale=False
        )
        decomposer.fit(spectra)
        
        # Create result
        result = DecompositionResult(
            mean_spectrum=decomposer.pca_model.mean_,
            components=decomposer.get_components(),
            explained_variance=decomposer.get_explained_variance(),
            explained_variance_ratio=decomposer.get_explained_variance_ratio(),
            config=DecompositionConfig(n_components=n_components),
            metadata={
                'n_reference_spectra': len(spectra),
                'n_channels': spectra.shape[1],
                'source': pca_source,
                'mission_id': mission_id,
            }
        )
        
        # Save results
        output_dir = Path("output/pca_components")
        output_dir.mkdir(parents=True, exist_ok=True)
        output_file = output_dir / f"decomposition_{mission_id}_{flight_date}_components.pkl"
        result.save(str(output_file))
        
        logger.info(f"✓ Results saved to {output_file}")
        logger.info(f"  Total variance explained: {result.explained_variance_ratio.sum():.2%}")
        logger.info(f"  Parameters used:")
        logger.info(f"    n_components: {n_components}")
        logger.info(f"    pca_source: {pca_source}")
        
        # Generate plots if requested
        if args.plot_components:
            logger.info("Generating component visualization plots...")
            try:
                # Load velocity axis from FITS for plotting
                with fits_io.open(fits_file) as hdul:
                    matrix_hdu = None
                    if 'MATRIX' in hdul:
                        matrix_hdu = hdul['MATRIX']
                    elif len(hdul) > 1:
                        for hdu in hdul[1:]:
                            if isinstance(hdu, fits_io.BinTableHDU):
                                matrix_hdu = hdu
                                break
                    
                    if matrix_hdu is not None and 'VELOCITY_AXIS' in matrix_hdu.data.dtype.names:
                        velocity_axis_ms = matrix_hdu.data['VELOCITY_AXIS'][0]
                        velocity_axis_kms = velocity_axis_ms / 1000.0
                        
                        _plot_decomposition_results(
                            result, velocity_axis_kms, mission_id, flight_date, output_dir
                        )
                        logger.info("✓ Plots generated successfully")
                    else:
                        logger.warning("Could not find VELOCITY_AXIS in FITS file, skipping plots")
            except Exception as e:
                logger.warning(f"Could not generate plots: {e}")
        
    except Exception as e:
        logger.error(f"Error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)