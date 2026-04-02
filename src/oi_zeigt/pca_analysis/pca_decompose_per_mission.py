#!/usr/bin/env python3
"""
PCA Decomposition per MISSION_ID/TELESCOPE combination

This script performs PCA decomposition on SKYCHOPDIFF reference spectra,
grouped by MISSION_ID and TELESCOP. Each mission/telescope combination gets 
its own set of PCA components, which is useful for correcting mission- and
instrument-specific sky contamination patterns.

VELOCITY AXIS HANDLING:
- Mission-specific telluric line parameters are in km/s (from mission_id_parameters.yml)
- The velocity axis is extracted from the FITS VELOCITY_AXIS column (in m/s)
- During decomposition: velocity ranges are converted to channel indices for masking
- The components themselves are eigenvectors stored in CHANNEL space (0-699)
- When applying correction: spectra from the SAME FITS file have the SAME velocity axis
  This ensures automatic alignment by velocity without manual conversion needed

Usage:
    python pca_decompose_per_mission.py --config config.toml [--n-components 5]

Output:
    Creates one pickle file per mission/telescope combination:
    output/pca_components/decomposition_<MISSION_ID>_<TELESCOP>_components.pkl
"""

import sys
import argparse
import logging
import pickle
from pathlib import Path
from typing import Dict, Tuple

import numpy as np
from astropy.io import fits

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

try:
    from oi_zeigt.pca_analysis.decompose import PCADecomposer, DecompositionResult, DecompositionConfig
except ImportError:
    logger.error("oi_zeigt module not found. Make sure it's installed with: pip install -e .")
    sys.exit(1)


def load_spectra_by_mission(fits_file: str, mission_params: dict = None) -> Dict[str, Dict]:
    """
    Load spectra from FITS file, grouped by MISSION_ID and TELESCOP.

    Drop rules from mission_id_parameters.yml are applied when ``mission_params``
    is provided:
      - ``drop.telescope`` — skip every combination whose TELESCOP matches.
      - ``drop.scans.complete`` — remove spectra with those SCAN numbers
        regardless of telescope.
      - ``drop.scans.telescope.<TELE>`` — remove those SCAN numbers only for
        the named telescope.

    Parameters
    ----------
    fits_file : str
        Path to FITS file with SPECTRA table containing MISSION_ID and TELESCOP columns.
    mission_params : dict, optional
        Parsed content of mission_id_parameters.yml.  When given, the ``drop``
        sub-section of each mission entry is applied.

    Returns
    -------
    dict
        Dictionary with structure:
        {
            'mission_id_1_telescope_1': {
                'spectra': array [n_spectra, n_channels],
                'mission_id': str,
                'telescop': str,
                'date': str (from first DATE-OBS),
                'metadata': {
                    'n_spectra': int,
                    'n_channels': int,
                    'source': str,
                }
            },
            'mission_id_1_telescope_2': {...},
            ...
        }
    """
    logger.info(f"Loading FITS file: {fits_file}")
    
    with fits.open(fits_file) as hdul:
        # Find the binary table HDU with MISSION_ID and TELESCOP columns
        # Accept HDUs with standard names OR empty names (from prepare_for_pca)
        spectra_hdu = None
        for hdu in hdul:
            # Check if this HDU has the required columns first
            if hasattr(hdu, 'data') and hdu.data is not None:
                if 'MISSION_ID' in hdu.data.dtype.names and 'TELESCOP' in hdu.data.dtype.names:
                    # Accept either standard names or empty name (from prepare_for_pca)
                    if hdu.name in ['SPECTRA', 'MATRIX', 'AVERAGE', '']:
                        spectra_hdu = hdu
                        break
        
        if spectra_hdu is None:
            raise ValueError("No HDU with both MISSION_ID and TELESCOP columns found in FITS file")
        
        data = spectra_hdu.data
        logger.info(f"Loaded {len(data)} total spectra from {spectra_hdu.name} HDU")
        
        # Extract mission IDs and telescopes
        def _col_to_str(x):
            if isinstance(x, bytes):
                return x.decode().strip()
            if isinstance(x, str):
                return x.strip()
            if isinstance(x, np.ndarray):
                return ''  # 18C complex array (CLASS artifact)
            return str(x).strip()

        mission_id_col = np.array([_col_to_str(x) for x in data['MISSION_ID']])
        telescop_col = np.array([_col_to_str(x) for x in data['TELESCOP']])
        
        # Get unique mission/telescope combinations
        unique_pairs = np.unique(list(zip(mission_id_col, telescop_col)), axis=0)
        logger.info(f"Found {len(unique_pairs)} unique MISSION_ID/TELESCOP combinations:")
        
        mission_data = {}
        
        for mission_id, telescop in unique_pairs:
            mask = (mission_id_col == mission_id) & (telescop_col == telescop)
            mission_spectra = data[mask]
            n_spectra = np.sum(mask)
            
            logger.info(f"  {mission_id}/{telescop}: {n_spectra} spectra")
            
            # Filter for SKYCHOPDIFF if OBJECT column exists
            if 'OBJECT' in data.dtype.names:
                objects = np.array([str(x).strip() if isinstance(x, bytes) else str(x).strip()
                                   for x in mission_spectra['OBJECT']])
                sky_mask = objects == 'SKYCHOPDIFF'
                sky_spectra = mission_spectra[sky_mask]['SPECTRUM']
                # Get the original indices of SKYCHOPDIFF spectra in the full data array
                mission_indices = np.where(mask)[0]
                sky_indices = mission_indices[sky_mask]
                logger.info(f"    └─ SKYCHOPDIFF: {np.sum(sky_mask)} spectra")
            else:
                sky_spectra = mission_spectra['SPECTRUM']
                sky_indices = np.where(mask)[0]
                logger.info(f"    └─ (no OBJECT filter, using all {n_spectra} spectra)")
            
            # Get date from first spectrum (YYYYMMDD format)
            if 'DATE-OBS' in mission_spectra.dtype.names:
                date_obs = str(mission_spectra['DATE-OBS'][0])
                if isinstance(date_obs, bytes):
                    date_obs = date_obs.decode().strip()
                # Extract just the date part (YYYY-MM-DD format or variations)
                date_part = date_obs.split('T')[0]
                # Split by dash to get year, month, day parts
                date_parts = date_part.split('-')
                if len(date_parts) == 3:
                    try:
                        # Reconstruct as YYYYMMDD with proper zero-padding
                        year = str(date_parts[0]).zfill(4)
                        month = str(date_parts[1]).zfill(2)
                        day = str(date_parts[2]).zfill(2)
                        flight_date = f"{year}{month}{day}"
                        if len(flight_date) != 8:
                            logger.warning(f"Unexpected date format: DATE-OBS={date_obs}")
                            flight_date = 'unknown'
                        else:
                            logger.debug(f"Extracted flight_date={flight_date} from DATE-OBS={date_obs}")
                    except (ValueError, AttributeError):
                        logger.warning(f"Could not parse date: DATE-OBS={date_obs}")
                        flight_date = 'unknown'
                else:
                    logger.warning(f"Unexpected date format (not YYYY-MM-DD): DATE-OBS={date_obs}")
                    flight_date = 'unknown'
            else:
                flight_date = 'unknown'
            
            # Create a unique key for this mission/telescope combination
            # Replace special characters in telescope name for filesystem safety
            safe_telescop = telescop.replace('/', '_').replace(' ', '_')
            key = f"{mission_id}_{safe_telescop}"
            
            # Store mission data
            mission_data[key] = {
                'spectra': np.array(sky_spectra),
                'indices': np.array(sky_indices),  # Track original indices in FITS file
                'mission_id': mission_id,
                'telescop': telescop,
                'date': flight_date,
                'metadata': {
                    'n_spectra': len(sky_spectra),
                    'n_channels': sky_spectra[0].shape[0] if len(sky_spectra) > 0 else 0,
                    'source': 'SKYCHOPDIFF' if 'OBJECT' in data.dtype.names else 'all',
                }
            }
        
        return mission_data


def preprocess_spectra(spectra: np.ndarray,
                       velocity_axis: np.ndarray = None):
    """
    Preprocess spectra: filter bad spectra and fill NaNs.

    Parameters
    ----------
    spectra : np.ndarray
        Spectra array [n_spectra, n_channels]
    velocity_axis : np.ndarray, optional
        Velocity axis in km/s, length n_channels. If given, trimmed to the
        same valid channel range as the spectra and returned as second value.

    Returns
    -------
    spectra : np.ndarray
        Preprocessed spectra with valid channels and NaNs filled.
    velocity_axis : np.ndarray or None
        Velocity axis trimmed to the valid channel range, or None if not given.
    """
    # Filter out spectra that are mostly NaN (>50% NaN)
    nan_fraction_per_spectrum = np.sum(np.isnan(spectra), axis=1) / spectra.shape[1]
    good_spectra_mask = nan_fraction_per_spectrum < 0.5
    original_count = len(spectra)
    spectra = spectra[good_spectra_mask]
    removed_count = original_count - len(spectra)
    if removed_count > 0:
        logger.info(f"Filtered out {removed_count} spectra with >50% NaN values")
    
    if len(spectra) == 0:
        logger.error("No spectra remain after NaN filtering!")
        raise ValueError("All spectra have >50% NaN values")
    
    # Find valid channels (not all NaN)
    valid_channels = ~np.all(np.isnan(spectra), axis=0)
    if not np.any(valid_channels):
        logger.error("No valid channels found!")
        raise ValueError("All channels are NaN in all spectra")
    
    original_n_channels = spectra.shape[1]
    first_valid = np.where(valid_channels)[0][0]
    last_valid = np.where(valid_channels)[0][-1]
    spectra = spectra[:, first_valid:last_valid+1]
    if velocity_axis is not None and len(velocity_axis) == original_n_channels:
        velocity_axis = velocity_axis[first_valid:last_valid + 1]

    n_valid_channels = last_valid - first_valid + 1
    logger.info(f"Trimmed to {n_valid_channels} valid channels (out of {original_n_channels} original)")
    
    # Fill remaining NaNs with per-channel mean
    n_nans_filled = 0
    for ch in range(spectra.shape[1]):
        channel_data = spectra[:, ch]
        if np.any(np.isnan(channel_data)):
            n_nans = np.sum(np.isnan(channel_data))
            n_nans_filled += n_nans
            mean_val = np.nanmean(channel_data)
            if np.isnan(mean_val):
                # All values are NaN in this channel, use 0
                mean_val = 0.0
            spectra[np.isnan(channel_data), ch] = mean_val
    
    if n_nans_filled > 0:
        logger.info(f"Filled {n_nans_filled} NaN values with per-channel means")
    
    # Final safety check: replace any remaining NaNs with 0
    remaining_nans = np.sum(np.isnan(spectra))
    if remaining_nans > 0:
        logger.warning(f"Found {remaining_nans} remaining NaNs after preprocessing, replacing with 0")
        spectra = np.nan_to_num(spectra, nan=0.0)
    
    # Final safety check: replace any inf or -inf
    spectra = np.nan_to_num(spectra, nan=0.0, posinf=0.0, neginf=0.0)

    return spectra, velocity_axis


def decompose_mission_spectra(mission_id: str, telescope: str, spectra: np.ndarray,
                              flight_date: str, n_components: int = 5,
                              velocity_axis: np.ndarray = None,
                              line_window_kms: tuple = None,
                              spectrum_indices: np.ndarray = None,
                              smoothing_kernel_size: int = None) -> DecompositionResult:
    """
    Perform PCA decomposition on spectra from a single mission/telescope combination.
    
    Parameters
    ----------
    mission_id : str
        Mission identifier
    telescope : str
        Telescope identifier
    spectra : np.ndarray
        Spectra array [n_spectra, n_channels]
    flight_date : str
        Flight date (YYYYMMDD format)
    n_components : int
        Number of PCA components to extract
    velocity_axis : np.ndarray, optional
        Velocity axis in km/s for masking lines
    line_window_kms : tuple, optional
        (v_min, v_max) in km/s to mask emission lines during decomposition
    spectrum_indices : np.ndarray, optional
        Indices of these spectra in the original FITS file
    smoothing_kernel_size : int, optional
        If set, apply a boxcar convolution of this length to each PCA component
        after fitting. Matches legacy smooth_pca_components() in pca_decompose.py.

    Returns
    -------
    DecompositionResult
        PCA decomposition result with components, explained variance, etc.
    """
    logger.info(f"\n{'='*80}")
    logger.info(f"DECOMPOSING: {mission_id} / {telescope}")
    logger.info(f"{'='*80}")
    logger.info(f"Spectra shape: {spectra.shape}")
    logger.info(f"N components: {n_components}")
    
    # Skip if no spectra
    if len(spectra) == 0:
        logger.warning(f"Skipping {mission_id}/{telescope}: no SKYCHOPDIFF spectra found")
        return None
    
    # Preprocess — also trims velocity_axis to match valid channel range
    logger.info(f"Preprocessing spectra...")
    spectra, velocity_axis = preprocess_spectra(spectra, velocity_axis)
    logger.info(f"After preprocessing: {spectra.shape}")
    
    # Baseline-subtract each SKYCHOPDIFF spectrum before decomposition,
    # matching the original pyclass behaviour where prepare_spectrum() calls
    # baseline() prior to adding a spectrum to the PCA input set.
    # The telluric line window is excluded from the baseline fit so that
    # the atmospheric emission feature does not bias the polynomial.
    from oi_zeigt.reduction.core import baseline_subtract
    baseline_window = None
    if line_window_kms is not None and velocity_axis is not None:
        v_min_bl, v_max_bl = line_window_kms
        line_mask_bl = (velocity_axis >= v_min_bl) & (velocity_axis <= v_max_bl)
        ch_indices = np.where(line_mask_bl)[0]
        if len(ch_indices) > 0:
            baseline_window = (int(ch_indices[0]), int(ch_indices[-1]) + 1)
    baselined = np.empty_like(spectra)
    for s_idx in range(len(spectra)):
        baselined[s_idx] = baseline_subtract(spectra[s_idx], order=1,
                                             window=baseline_window)
    spectra = baselined
    logger.info(f"  ✓ Baseline-subtracted {len(spectra)} spectra "
                f"(order=1, excluded channels {baseline_window})")

    # Apply telluric line masking (from mission-specific parameters)
    if line_window_kms is not None and velocity_axis is not None:
        logger.info(f"Masking telluric line window {line_window_kms[0]:.1f}-{line_window_kms[1]:.1f} km/s for decomposition")
        logger.info(f"  Velocity axis range: {velocity_axis[0]:.1f} to {velocity_axis[-1]:.1f} km/s ({len(velocity_axis)} channels)")
        # Find channels in the telluric line window
        v_min, v_max = line_window_kms
        line_mask = (velocity_axis >= v_min) & (velocity_axis <= v_max)
        line_channels = np.where(line_mask)[0]
        
        if len(line_channels) > 0:
            logger.info(f"  Mapped to {len(line_channels)} channels: indices {line_channels[0]}-{line_channels[-1]}")
            logger.info(f"  Actual velocity range in data: {velocity_axis[line_channels[0]]:.1f} to {velocity_axis[line_channels[-1]]:.1f} km/s")
            # Mask out line regions by setting to mean of continuum
            for ch in line_channels:
                if ch < spectra.shape[1]:
                    # Set masked channels to mean of unmasked channels
                    continuum = spectra[:, ~np.isin(np.arange(spectra.shape[1]), line_channels)]
                    spectra[:, ch] = np.mean(continuum, axis=1) if continuum.size > 0 else 0
            logger.info(f"  ✓ Masked {len(line_channels)} channels in line region")
        else:
            logger.warning(f"  ⚠ No channels found in line window {v_min:.1f}-{v_max:.1f} km/s!")
    
    # Fit PCA
    logger.info(f"Fitting PCA with {n_components} components...")
    decomposer = PCADecomposer(n_components=n_components, scale=False)
    decomposer.fit(spectra)

    # Smooth PCA components with a boxcar kernel (matches legacy smooth_pca_components()
    # in pca_decompose.py, which used np.convolve with mode='same').
    if smoothing_kernel_size:
        kernel_size = int(smoothing_kernel_size)
        kernel = np.ones(kernel_size) / kernel_size
        smoothed = []
        for component in decomposer.pca_model.components_:
            smoothed.append(np.convolve(component.copy(), kernel, mode='same'))
        decomposer.pca_model.components_ = np.array(smoothed)
        logger.info(f"Smoothed {len(smoothed)} components with boxcar kernel size {kernel_size}")

    # Create result
    if spectrum_indices is not None:
        logger.info(f"Creating DecompositionResult with {len(spectrum_indices)} spectrum indices")
    else:
        logger.info(f"Creating DecompositionResult with NO spectrum indices (spectrum_indices is None)")
    
    result = DecompositionResult(
        mean_spectrum=decomposer.pca_model.mean_,
        components=decomposer.get_components(),
        explained_variance=decomposer.get_explained_variance(),
        explained_variance_ratio=decomposer.get_explained_variance_ratio(),
        config=DecompositionConfig(n_components=n_components),
        metadata={
            'n_reference_spectra': len(spectra),
            'n_channels': spectra.shape[1],
            'source': 'SKYCHOPDIFF',
            'mission_id': mission_id,
            'velocity_axis_kms': velocity_axis,
        },
        reference_spectrum_indices=spectrum_indices  # Store indices to retrieve spectra from original FITS
    )
    
    # Log results
    logger.info(f"\nResults:")
    logger.info(f"  Components shape: {result.components.shape}")
    logger.info(f"  Explained variance ratio: {result.explained_variance_ratio}")
    logger.info(f"  Total variance explained: {result.explained_variance_ratio.sum():.2%}")
    
    for i, ratio in enumerate(result.explained_variance_ratio, 1):
        logger.info(f"    Component {i}: {ratio*100:6.2f}%")
    
    return result


def save_results(result: DecompositionResult, mission_id: str, telescop: str) -> Path:
    """
    Save decomposition result to pickle file.
    
    Parameters
    ----------
    result : DecompositionResult
        Decomposition result to save
    mission_id : str
        Mission identifier for filename (includes date like 2017-02-01_...)
    telescop : str
        Telescope identifier for filename
    
    Returns
    -------
    Path
        Path to saved pickle file
    """
    output_dir = Path("output/pca_components")
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Replace special characters in telescope name for filesystem safety
    safe_telescop = telescop.replace('/', '_').replace(' ', '_')
    output_file = output_dir / f"decomposition_{mission_id}_{safe_telescop}_components.pkl"
    
    with open(output_file, 'wb') as f:
        pickle.dump({
            'components': result.components,
            'mean_spectrum': result.mean_spectrum,
            'explained_variance': result.explained_variance,
            'explained_variance_ratio': result.explained_variance_ratio,
            'config': result.config.to_dict(),
            'metadata': result.metadata,
            'spectrum_metadata': result.spectrum_metadata,
            'reference_spectrum_indices': result.reference_spectrum_indices,  # Include spectrum indices
        }, f)
    
    logger.info(f"✓ Saved: {output_file}")
    return output_file


def main_cli():
    """Command-line interface for per-mission/telescope PCA decomposition."""
    parser = argparse.ArgumentParser(
        description="PCA decomposition of spectral reference data per MISSION_ID/TELESCOPE combination",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Input FITS file resolution order (first found is used):\n"
            "  1. --fits argument\n"
            "  2. [output].prepared_for_pca from config.toml  (default)\n"
            "  3. [output].reduced_fits from config.toml\n"
            "  4. [input].fits_file from config.toml\n\n"
            "Example:\n"
            "  pca_decompose --config config.toml --n-components 5"
        )
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
        help="FITS file to decompose (overrides config file [output][prepared_for_pca] or [output][reduced_fits])"
    )
    parser.add_argument(
        "--n-components",
        type=int,
        default=None,
        help="Number of PCA components to extract (default: from config [pca][n_components], or 5)"
    )
    parser.add_argument(
        "--mission-id",
        type=str,
        default=None,
        help="Limit decomposition to this MISSION_ID (e.g. 2016-05-18_GR_F298). "
             "Can be combined with --telescope."
    )
    parser.add_argument(
        "--telescope",
        type=str,
        default=None,
        help="Limit decomposition to this TELESCOP value (e.g. LFAH_PX00_S). "
             "Can be combined with --mission-id."
    )
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Verbose output"
    )
    parser.add_argument(
        "--plot-components",
        action="store_true",
        help="Generate visualization plots of PCA components"
    )
    
    args = parser.parse_args()
    
    # Set logging level
    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)
    
    try:
        # Load configuration
        from oi_zeigt.basic_io import get_config
        logger.info(f"Loading configuration from {args.config}")
        config = get_config(args.config)
        
        # Resolve n_components: CLI > config [pca][n_components] > default 5
        pca_config = config.get('pca', {})
        if args.n_components is not None:
            n_components = args.n_components
            logger.info(f"✓ N components = {n_components} (from command line)")
        elif pca_config.get('n_components') is not None:
            n_components = int(pca_config['n_components'])
            logger.info(f"✓ N components = {n_components} (from config [pca][n_components])")
        else:
            n_components = 5
            logger.info(f"✓ N components = {n_components} (default)")

        # Determine FITS file (command line > config [output][prepared_for_pca] > config [output][reduced_fits] > config [input][fits_file])
        if args.fits is not None:
            fits_file = args.fits
            logger.info(f"✓ FITS file = {fits_file} (from command line)")
        else:
            # Try [output][prepared_for_pca] first
            output_config = config.get('output', {})
            fits_file = output_config.get('prepared_for_pca')
            
            if fits_file:
                logger.info(f"✓ FITS file = {fits_file} (from config [output][prepared_for_pca])")
            else:
                # Fall back to [output][reduced_fits]
                fits_file = output_config.get('reduced_fits')
                
                if fits_file:
                    logger.info(f"✓ FITS file = {fits_file} (from config [output][reduced_fits])")
                else:
                    # Final fall back to [input][fits_file]
                    input_config = config.get('input', {})
                    fits_file = input_config.get('fits_file')
                    if fits_file:
                        logger.info(f"✓ FITS file = {fits_file} (from config [input][fits_file])")
                    else:
                        logger.error("No FITS file specified in configuration or command line")
                        logger.error("Specify via:")
                        logger.error("  1. Command line: pca_decompose --fits /path/to/file.fits")
                        logger.error("  2. Config [output][prepared_for_pca]")
                        logger.error("  3. Config [output][reduced_fits]")
                        logger.error("  4. Config [input][fits_file]")
                        sys.exit(1)
        
        # Verify file exists
        if not Path(fits_file).exists():
            logger.error(f"FITS file not found: {fits_file}")
            sys.exit(1)
        
        # Load spectra by mission
        mission_data = load_spectra_by_mission(fits_file)
        
        # Extract velocity axis from FITS file
        velocity_axis = None
        try:
            from oi_zeigt.basic_io import reconstruct_velocity_axis
            with fits.open(fits_file) as hdul:
                for hdu in hdul:
                    if hasattr(hdu, 'data') and hdu.data is not None and hasattr(hdu.data, 'dtype'):
                        names = hdu.data.dtype.names or []
                        if 'VELOCITY' in names and 'DELTAV' in names and 'SPECTRUM' in names:
                            vel_axis_ms = reconstruct_velocity_axis(hdu)
                            velocity_axis = vel_axis_ms / 1000.0  # Convert m/s to km/s
                            logger.info(f"✓ Reconstructed velocity axis: {len(velocity_axis)} channels, "
                                      f"range: {velocity_axis[0]:.1f} to {velocity_axis[-1]:.1f} km/s")
                            break
        except Exception as e:
            logger.warning(f"Could not reconstruct velocity axis from FITS: {e}")
            logger.warning(f"Proceeding without velocity axis (line masking will be skipped)")
        
        # Load mission-specific parameters from mission_id_parameters.yml
        mission_line_windows = {}
        try:
            mission_yml = Path(__file__).parent / "mission_id_parameters.yml"
            if mission_yml.exists():
                import yaml
                with open(mission_yml, 'r') as f:
                    mission_params = yaml.safe_load(f) or {}
                
                # Extract line window parameters for each mission
                for mission_id in set(d['mission_id'] for d in mission_data.values()):
                    if mission_id in mission_params:
                        params = mission_params[mission_id]
                        if 'telluric_line_center' in params and 'telluric_line_width' in params:
                            center = params['telluric_line_center']
                            width = params['telluric_line_width']
                            v_min = center - width / 2.0
                            v_max = center + width / 2.0
                            mission_line_windows[mission_id] = (v_min, v_max)
                            logger.info(f"✓ {mission_id}: telluric line {v_min:.1f}-{v_max:.1f} km/s")
        except Exception as e:
            logger.warning(f"Could not load mission parameters: {e}")
            logger.warning(f"Proceeding without line masking")
        
        if not mission_data:
            logger.error("No MISSION_ID data found in FITS file")
            sys.exit(1)

        # Filter to a single mission_id / telescope if requested
        if args.mission_id or args.telescope:
            before = len(mission_data)
            mission_data_full = mission_data  # keep for error reporting
            mission_data = {
                k: v for k, v in mission_data.items()
                if (args.mission_id is None or v['mission_id'] == args.mission_id)
                and (args.telescope is None or v['telescop'] == args.telescope)
            }
            if not mission_data:
                logger.error(
                    f"No combinations matched mission_id={args.mission_id!r} "
                    f"telescope={args.telescope!r} (had {before} total)."
                )
                # Show available values to help the user correct the typo
                all_missions  = sorted({v['mission_id'] for v in mission_data_full.values()})
                all_telescopes = sorted({v['telescop']   for v in mission_data_full.values()})
                if args.mission_id:
                    import difflib
                    close = difflib.get_close_matches(args.mission_id, all_missions, n=5, cutoff=0.6)
                    if close:
                        logger.error(f"  Did you mean one of: {', '.join(close)}")
                    else:
                        logger.error(f"  Available mission_ids (first 10): {', '.join(all_missions[:10])}")
                if args.telescope:
                    import difflib
                    close = difflib.get_close_matches(args.telescope, all_telescopes, n=5, cutoff=0.6)
                    if close:
                        logger.error(f"  Did you mean telescope: {', '.join(close)}")
                    else:
                        logger.error(f"  Available telescopes: {', '.join(all_telescopes)}")
                sys.exit(1)
            logger.info(
                f"Filtered to {len(mission_data)}/{before} combination(s): "
                f"mission_id={args.mission_id!r}  telescope={args.telescope!r}"
            )

        # Check for missions with no SKYCHOPDIFF spectra
        missions_no_sky = [m for m, d in mission_data.items() if d['metadata']['n_spectra'] == 0]
        if missions_no_sky:
            logger.warning(f"Missions with NO SKYCHOPDIFF spectra: {', '.join(str(m) for m in missions_no_sky)}")

        logger.info(f"\n{'='*80}")
        logger.info(f"PROCESSING {len(mission_data)} MISSION/TELESCOPE COMBINATIONS")
        logger.info(f"{'='*80}\n")
        
        # Decompose each mission/telescope combination
        results = {}
        output_files = []
        
        for key in sorted(mission_data.keys()):
            data = mission_data[key]
            mission_id = data['mission_id']
            telescop = data['telescop']
            spectrum_indices = data.get('indices')
            if spectrum_indices is not None:
                logger.info(f"Using {len(spectrum_indices)} spectrum indices for {mission_id}/{telescop}")
            
            # Decompose with velocity axis and line masking if available
            result = decompose_mission_spectra(
                mission_id,
                telescop,
                data['spectra'],
                data['date'],
                n_components=n_components,
                velocity_axis=velocity_axis,
                line_window_kms=mission_line_windows.get(mission_id),
                spectrum_indices=spectrum_indices,
                smoothing_kernel_size=pca_config.get('smoothing_kernel_size')
            )
            
            if result is None:
                continue
            
            results[key] = result
            
            # Save
            output_file = save_results(result, mission_id, telescop)
            output_files.append(output_file)
        
        # Summary
        logger.info(f"\n{'='*80}")
        logger.info(f"SUMMARY")
        logger.info(f"{'='*80}")
        logger.info(f"Processed {len(results)} mission/telescope combinations")
        logger.info(f"Output files:")
        for output_file in output_files:
            logger.info(f"  ✓ {output_file}")
        
        # Generate plots if requested
        if args.plot_components and results:
            logger.info(f"\nGenerating component visualization plots...")
            try:
                import matplotlib.pyplot as plt
                import matplotlib
                # Use non-interactive backend
                matplotlib.use('Agg')
                
                output_dir = Path("output/pca_components")
                
                for key in sorted(results.keys()):
                    result = results[key]
                    if key in mission_data:
                        mission_id = mission_data[key]['mission_id']
                        telescop = mission_data[key]['telescop']
                        flight_date = mission_data[key]['date']
                        
                        try:
                            components = result.components
                            explained_variance_ratio = result.explained_variance_ratio
                            mean_spectrum = result.mean_spectrum
                            n_components = components.shape[0]
                            n_channels = components.shape[1]

                            # Use the trimmed velocity axis stored in result metadata
                            vax = result.metadata.get('velocity_axis_kms')
                            if vax is not None and len(vax) == n_channels:
                                x_axis = vax
                                x_label = 'Velocity (km/s)'
                            else:
                                x_axis = np.arange(n_channels)
                                x_label = 'Channel'

                            # Create component plots
                            fig, axes = plt.subplots(n_components + 1, 1, figsize=(14, 3*(n_components + 1)))

                            # Plot mean spectrum
                            axes[0].plot(x_axis, mean_spectrum, 'b-', linewidth=1.5)
                            axes[0].set_title(f'Mean Spectrum - {mission_id}/{telescop}', fontsize=12, fontweight='bold')
                            axes[0].set_ylabel('Intensity')
                            axes[0].set_xlim(x_axis[0], x_axis[-1])
                            axes[0].grid(True, alpha=0.3)

                            # Plot each component
                            for i in range(n_components):
                                axes[i+1].plot(x_axis, components[i], 'r-', linewidth=1.5)
                                axes[i+1].set_title(
                                    f'Component {i+1} (Variance: {100*explained_variance_ratio[i]:.2f}%) - {mission_id}/{telescop}',
                                    fontsize=12, fontweight='bold'
                                )
                                axes[i+1].set_ylabel('Loadings')
                                axes[i+1].set_xlim(x_axis[0], x_axis[-1])
                                axes[i+1].grid(True, alpha=0.3)
                                if i == n_components - 1:
                                    axes[i+1].set_xlabel(x_label)
                            
                            plt.tight_layout()
                            
                            # Safe telescope name for filename
                            safe_telescop = telescop.replace('/', '_').replace(' ', '_')
                            components_path = output_dir / f"pca_components_{mission_id}_{safe_telescop}.png"
                            plt.savefig(components_path, dpi=150, bbox_inches='tight')
                            plt.close(fig)
                            logger.info(f"✓ Saved component plots to {components_path}")
                            
                            # Create variance explained plot
                            fig2, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
                            
                            # Individual variance
                            ax1.bar(range(1, n_components + 1), 100 * explained_variance_ratio, 
                                   color='steelblue', alpha=0.7, edgecolor='black')
                            ax1.set_xlabel('Component', fontsize=12)
                            ax1.set_ylabel('Explained Variance (%)', fontsize=12)
                            ax1.set_title(f'Explained Variance per Component - {mission_id}/{telescop}', 
                                         fontsize=12, fontweight='bold')
                            ax1.grid(True, alpha=0.3)
                            for i, v in enumerate(100 * explained_variance_ratio, 1):
                                ax1.text(i, v + 0.5, f'{v:.1f}%', ha='center', fontsize=10)
                            
                            # Cumulative variance
                            cumsum = np.cumsum(100 * explained_variance_ratio)
                            ax2.plot(range(1, n_components + 1), cumsum, 'o-', linewidth=2, markersize=8, color='darkgreen')
                            ax2.axhline(y=90, color='red', linestyle='--', linewidth=2, label='90% threshold')
                            ax2.set_xlabel('Number of Components', fontsize=12)
                            ax2.set_ylabel('Cumulative Explained Variance (%)', fontsize=12)
                            ax2.set_title(f'Cumulative Explained Variance - {mission_id}/{telescop}', 
                                         fontsize=12, fontweight='bold')
                            ax2.grid(True, alpha=0.3)
                            ax2.legend()
                            ax2.set_ylim(0, 105)
                            
                            plt.tight_layout()
                            variance_path = output_dir / f"pca_variance_{mission_id}_{safe_telescop}.png"
                            plt.savefig(variance_path, dpi=150, bbox_inches='tight')
                            plt.close(fig2)
                            logger.info(f"✓ Saved variance plots to {variance_path}")
                            
                        except Exception as e:
                            logger.warning(f"Could not generate plots for {mission_id}/{telescop}: {e}")
                    
            except Exception as e:
                logger.warning(f"Could not generate plots: {e}")
        
        logger.info(f"\n✓ All decompositions completed successfully!")
        
    except Exception as e:
        logger.error(f"Error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main_cli()
