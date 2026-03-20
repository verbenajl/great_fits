"""
Prepare FITS data for PCA analysis.

This module filters FITS data to include only specific sources and objects,
and fills telluric line regions with Gaussian noise.
"""

import logging
import os
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Optional

import numpy as np
import yaml
from astropy.io import fits

logger = logging.getLogger(__name__)


def load_mission_parameters(mission_id: str) -> dict:
    """
    Load mission-specific parameters from YAML file.
    
    Parameters
    ----------
    mission_id : str
        Mission identifier
    
    Returns
    -------
    dict
        Mission parameters including telluric_line_center and telluric_line_width
    """
    yaml_file = Path(__file__).parent / 'mission_id_parameters.yml'
    
    if not yaml_file.exists():
        logger.warning(f"Mission parameters file not found: {yaml_file}")
        return {}
    
    try:
        with open(yaml_file, 'r') as f:
            config = yaml.safe_load(f)
        
        # Find matching mission ID (may be a substring)
        for key in config:
            if key and mission_id in key:
                logger.debug(f"Found mission config for '{mission_id}': {key}")
                return config[key] or {}
        
        logger.warning(f"Mission '{mission_id}' not found in parameters file")
        return {}
    except Exception as e:
        logger.warning(f"Error loading mission parameters: {e}")
        return {}


def get_telluric_indices(mission_id: str, velocity_axis_kms: np.ndarray) -> Optional[np.ndarray]:
    """
    Get channel indices for telluric line regions based on mission parameters.
    
    Parameters
    ----------
    mission_id : str
        Mission identifier
    velocity_axis_kms : ndarray
        Velocity axis in km/s
    
    Returns
    -------
    ndarray or None
        Boolean array of channel indices to fill with noise
    """
    params = load_mission_parameters(mission_id)
    
    if not params or 'telluric_line_center' not in params:
        return None
    
    center_km_s = params['telluric_line_center']
    width_km_s = params.get('telluric_line_width', 30)
    
    # Create mask for velocities within the line region
    v_min = center_km_s - width_km_s / 2.0
    v_max = center_km_s + width_km_s / 2.0
    
    telluric_mask = (velocity_axis_kms >= v_min) & (velocity_axis_kms <= v_max)
    
    if np.any(telluric_mask):
        indices = np.where(telluric_mask)[0]
        logger.info(f"Telluric mask for {mission_id}: center={center_km_s} km/s, width={width_km_s} km/s")
        logger.info(f"  Channels {indices[0]}-{indices[-1]} to be filled with noise")
        return telluric_mask
    
    return None


def _fill_telluric_chunk(args: tuple) -> np.ndarray:
    """
    Worker function: fill telluric channels with Gaussian noise for a chunk of spectra.

    Must be defined at module level so ProcessPoolExecutor can pickle it.

    Parameters
    ----------
    args : tuple
        (spectra_chunk, velocity_axes_chunk, mission_ids_chunk, objects_chunk, mission_params_cache)
        - spectra_chunk       : (N, nchan) float array
        - velocity_axes_chunk : (N, nchan) float array in km/s
        - mission_ids_chunk   : list of str
        - objects_chunk       : list of str
        - mission_params_cache: dict {mission_id -> params dict}

    Returns
    -------
    np.ndarray
        (N, nchan) array with telluric channels filled.
    """
    spectra_chunk, velocity_axes_chunk, mission_ids_chunk, objects_chunk, mission_params_cache = args
    result = spectra_chunk.copy()

    for i in range(len(result)):
        if objects_chunk[i] in ('TSYS', 'TAU_SIG'):
            continue

        params = mission_params_cache.get(mission_ids_chunk[i], {})
        if not params or 'telluric_line_center' not in params:
            continue

        center_km_s = params['telluric_line_center']
        width_km_s = params.get('telluric_line_width', 30)
        v_min = center_km_s - width_km_s / 2.0
        v_max = center_km_s + width_km_s / 2.0

        velocity_kms = velocity_axes_chunk[i]
        telluric_mask = (velocity_kms >= v_min) & (velocity_kms <= v_max)

        if not np.any(telluric_mask):
            continue

        spectrum = result[i]
        non_telluric_values = spectrum[~telluric_mask]
        n_valid = np.sum(~np.isnan(non_telluric_values))
        noise_level = np.nanstd(non_telluric_values) if n_valid >= 2 else np.nan

        if np.isnan(noise_level):
            n_valid_full = np.sum(~np.isnan(spectrum))
            noise_level = np.nanstd(spectrum) if n_valid_full >= 2 else np.nan

        if not np.isnan(noise_level) and noise_level > 0:
            result[i][telluric_mask] = np.random.normal(0, noise_level, np.sum(telluric_mask))

    return result


def fill_telluric_with_noise(fits_file: str, output_fits: str,
                            pca_source: str = "SKYCHOPDIFF",
                            object_filter: str = "M51CENTER",
                            mission_id: Optional[str] = None,
                            scan: Optional[int] = None,
                            fill_noise: bool = False) -> None:
    """
    Filter FITS file and optionally fill telluric lines with Gaussian noise.

    Also associates each M51CENTER spectrum with its corresponding TSYS and TAU_SIG
    measurements, storing the indices as new columns in the output.

    Parameters
    ----------
    fits_file : str
        Input FITS file path
    output_fits : str
        Output FITS file path
    pca_source : str
        Source to filter for (default: SKYCHOPDIFF)
    object_filter : str
        Object substring to filter for (default: M51CENTER)
    mission_id : str, optional
        Filter to specific mission_id (e.g., "2017-02-01_GR_F367")
    scan : int, optional
        Filter to specific SCAN number
    fill_noise : bool, optional
        If True, fill telluric channels with Gaussian noise.
        If False (default), leave telluric channels unchanged.
    """
    logger.info(f"Opening {fits_file}...")
    
    with fits.open(fits_file) as hdul:
        header = hdul[1].header
        data = hdul[1].data
        
        # Get velocity axis from data columns
        n_channels = len(data['SPECTRUM'][0])
        
        try:
            # Try to get VELOCITY_AXIS column first (velocity per channel)
            if 'VELOCITY_AXIS' in data.dtype.names:
                velocity_axis_ms = data['VELOCITY_AXIS'][0]  # m/s per channel
                velocity_axis_kms = velocity_axis_ms / 1000.0  # Convert to km/s
                logger.info(f"Velocity axis from VELOCITY_AXIS column: {len(velocity_axis_kms)} channels, {velocity_axis_kms[0]:.1f} to {velocity_axis_kms[-1]:.1f} km/s")
            # Fall back to constructing from VELOCITY and DELTAV
            elif ('VELOCITY' in data.dtype.names and 'DELTAV' in data.dtype.names):
                velo_ref = float(data['VELOCITY'][0])
                deltav = float(data['DELTAV'][0])
                
                # Get reference pixel from header
                crpix1_spec = float(header.get('CRPIX1', 1.0))
                
                # Create velocity axis
                channel_indices = np.arange(n_channels, dtype=np.float64)
                velocity_axis = velo_ref + (channel_indices - (crpix1_spec - 1.0)) * deltav
                velocity_axis_kms = velocity_axis / 1000.0  # Convert m/s to km/s
                logger.info(f"Velocity axis from VELOCITY/DELTAV: {len(velocity_axis_kms)} channels, {velocity_axis_kms[0]:.1f} to {velocity_axis_kms[-1]:.1f} km/s")
            else:
                raise KeyError("VELOCITY_AXIS, VELOCITY, or DELTAV columns missing")
        except (AttributeError, KeyError, TypeError, ValueError) as e:
            logger.warning(f"Could not construct velocity axis: {e}")
            velocity_axis_kms = np.arange(n_channels)
        
        # Filter data
        logger.info(f"Filtering to OBJECT = '{pca_source}' or '{object_filter}' or TSYS/TAU_SIG...")
        
        # Get masks for filtering
        def _to_str(s):
            if isinstance(s, (bytes, str)):
                return s.decode().strip() if isinstance(s, bytes) else s.strip()
            return ''  # e.g. 18C complex array (CLASS artifact) — treat as empty

        objects = np.array([s.strip() if isinstance(s, str) else
                            s.decode().strip() if isinstance(s, bytes) else ''
                            for s in data['OBJECT']])
        mission_ids = np.array([_to_str(s) for s in data['MISSION_ID']])
        
        # Keep spectra where OBJECT equals either pca_source or object_filter
        # Also include TSYS and TAU_SIG rows
        combined_mask = (objects == pca_source) | (objects == object_filter) | \
                       (objects == 'TSYS') | (objects == 'TAU_SIG')
        
        # Apply mission_id filter if specified
        if mission_id:
            logger.info(f"Filtering to MISSION_ID = '{mission_id}'...")
            mission_filter = mission_ids == mission_id
            before_filter = np.sum(combined_mask)
            combined_mask = combined_mask & mission_filter
            after_filter = np.sum(combined_mask)
            logger.info(f"  After mission_id filter: {after_filter}/{before_filter} rows kept")
        
        # Apply scan filter if specified
        if scan is not None:
            logger.info(f"Filtering to SCAN = {scan}...")
            scans = np.array(data['SCAN'])
            scan_filter = scans == scan
            before_filter = np.sum(combined_mask)
            combined_mask = combined_mask & scan_filter
            after_filter = np.sum(combined_mask)
            logger.info(f"  After scan filter: {after_filter}/{before_filter} rows kept")
        
        n_original = len(data)
        n_filtered = np.sum(combined_mask)
        n_tsys = np.sum(objects[combined_mask] == 'TSYS')
        n_tau_sig = np.sum(objects[combined_mask] == 'TAU_SIG')
        logger.info(f"Filtered from {n_original} to {n_filtered} spectra")
        logger.info(f"  - {np.sum((objects[combined_mask] == pca_source) | (objects[combined_mask] == object_filter))} science spectra")
        logger.info(f"  - {n_tsys} TSYS spectra")
        logger.info(f"  - {n_tau_sig} TAU_SIG spectra")
        
        # Filter all data
        filtered_data = data[combined_mask]
        filtered_mission_ids = mission_ids[combined_mask]
        filtered_objects = objects[combined_mask]
        
        if fill_noise:
            # Fill telluric lines with Gaussian noise (parallelized)
            logger.info("Filling telluric lines with Gaussian noise...")

            # Pre-cache mission parameters — avoids re-parsing the YAML file for every spectrum
            unique_mission_ids = list(set(filtered_mission_ids))
            mission_params_cache = {mid: load_mission_parameters(mid) for mid in unique_mission_ids}
            logger.info(f"  Cached parameters for {len(mission_params_cache)} unique mission IDs")

            # Build per-spectrum velocity axes as a plain numpy array for parallel workers
            spectra_array = np.array(filtered_data['SPECTRUM'], dtype=np.float64)
            n_spectra, n_chan = spectra_array.shape

            if 'VELOCITY_AXIS' in filtered_data.dtype.names:
                velocity_axes = np.array(filtered_data['VELOCITY_AXIS'], dtype=np.float64) / 1000.0
            else:
                velocity_axes = np.tile(velocity_axis_kms, (n_spectra, 1))

            # Split work into chunks — one per CPU
            n_workers = os.cpu_count() or 1
            chunk_size = max(1, n_spectra // n_workers)
            chunks = []
            for start in range(0, n_spectra, chunk_size):
                end = min(start + chunk_size, n_spectra)
                chunks.append((
                    spectra_array[start:end],
                    velocity_axes[start:end],
                    list(filtered_mission_ids[start:end]),
                    list(filtered_objects[start:end]),
                    mission_params_cache,
                ))

            logger.info(f"  Processing {n_spectra} spectra across {len(chunks)} chunks ({n_workers} workers)...")

            filled_chunks = []
            with ProcessPoolExecutor(max_workers=n_workers) as executor:
                for chunk_result in executor.map(_fill_telluric_chunk, chunks):
                    filled_chunks.append(chunk_result)

            # Reassemble and count filled spectra
            filled_spectra = np.concatenate(filled_chunks, axis=0)
            n_filled = int(np.sum(
                np.any(filled_spectra != spectra_array, axis=1) &
                np.array([obj not in ('TSYS', 'TAU_SIG') for obj in filtered_objects])
            ))

            # Write back into filtered_data in-place
            filtered_data['SPECTRUM'][:] = filled_spectra.astype(filtered_data['SPECTRUM'].dtype)

            logger.info(f"Filled telluric lines in {n_filled} spectra")
        else:
            logger.info("Skipping telluric noise filling (--fill-telluric-with-noise not set)")
        
        # Create index columns for linking M51CENTER to TSYS and TAU_SIG
        logger.info("Creating TSYS and TAU_SIG index columns...")
        
        # Initialize index columns with -1 (indicating no match)
        tsys_indices = np.full(len(filtered_data), -1, dtype=np.int32)
        tau_sig_indices = np.full(len(filtered_data), -1, dtype=np.int32)
        
        # Build structured lookup for TSYS and TAU_SIG with subscan info
        # Structure: {(mission_id, scan) -> [(subscan, idx, ut_time), ...]}
        tsys_lookup = {}  # (mission_id, scan) -> [(subscan, idx, ut_time), ...]
        tau_sig_lookup = {}  # (mission_id, scan) -> [(subscan, idx, ut_time), ...]
        
        # Check if UT or LST columns exist for temporal fallback
        has_ut_col = 'UT' in filtered_data.dtype.names
        has_lst_col = 'LST' in filtered_data.dtype.names
        time_col = 'UT' if has_ut_col else ('LST' if has_lst_col else None)
        
        for idx in range(len(filtered_data)):
            obj = filtered_objects[idx]
            if obj in ('TSYS', 'TAU_SIG'):
                key = (filtered_mission_ids[idx], filtered_data['SCAN'][idx])
                subscan = filtered_data['SUBSCAN'][idx]
                
                # Build timestamp for temporal fallback (UT or LST in fractional days)
                timestamp = None
                if time_col:
                    try:
                        timestamp = float(filtered_data[time_col][idx])
                    except:
                        pass
                
                entry = (subscan, idx, timestamp)
                
                if obj == 'TSYS':
                    if key not in tsys_lookup:
                        tsys_lookup[key] = []
                    tsys_lookup[key].append(entry)
                else:  # TAU_SIG
                    if key not in tau_sig_lookup:
                        tau_sig_lookup[key] = []
                    tau_sig_lookup[key].append(entry)
        
        logger.info(f"Found {len(tsys_lookup)} TSYS groups and {len(tau_sig_lookup)} TAU_SIG groups")
        
        # Link M51CENTER spectra to TSYS and TAU_SIG
        n_tsys_linked = 0
        n_tsys_fallback = 0
        n_tau_sig_linked = 0
        n_tau_sig_fallback = 0
        
        for idx in range(len(filtered_data)):
            if filtered_objects[idx] == object_filter:  # M51CENTER
                key = (filtered_mission_ids[idx], filtered_data['SCAN'][idx])
                m51_subscan = filtered_data['SUBSCAN'][idx]
                m51_timestamp = None
                
                if time_col:
                    try:
                        m51_timestamp = float(filtered_data[time_col][idx])
                    except:
                        pass
                
                # Link to TSYS using pattern-aware matching
                if key in tsys_lookup and len(tsys_lookup[key]) > 0:
                    # Pattern-aware: look for TSYS with subscan = M51CENTER_subscan - 1
                    tsys_entries = tsys_lookup[key]
                    
                    # First try: exact pattern match (TSYS subscan = M51CENTER subscan - 1)
                    pattern_matches = [entry for entry in tsys_entries if entry[0] == m51_subscan - 1]
                    
                    if pattern_matches:
                        # Use the first pattern match
                        tsys_indices[idx] = pattern_matches[0][1]
                        n_tsys_linked += 1
                        logger.debug(f"M51CENTER idx={idx} subscan={m51_subscan} -> TSYS idx={pattern_matches[0][1]} subscan={pattern_matches[0][0]} (pattern match)")
                    else:
                        # Fallback: use nearest TSYS by time if available, else use first TSYS
                        if m51_timestamp is not None and any(entry[2] is not None for entry in tsys_entries):
                            # Find TSYS with closest timestamp
                            valid_entries = [e for e in tsys_entries if e[2] is not None]
                            nearest = min(valid_entries, key=lambda e: abs(e[2] - m51_timestamp))
                            tsys_indices[idx] = nearest[1]
                            n_tsys_fallback += 1
                            logger.debug(f"M51CENTER idx={idx} subscan={m51_subscan} -> TSYS idx={nearest[1]} subscan={nearest[0]} (temporal fallback)")
                        else:
                            # Last resort: use first TSYS in group
                            tsys_indices[idx] = tsys_entries[0][1]
                            n_tsys_fallback += 1
                            logger.debug(f"M51CENTER idx={idx} subscan={m51_subscan} -> TSYS idx={tsys_entries[0][1]} subscan={tsys_entries[0][0]} (first match fallback)")
                
                # Link to TAU_SIG using pattern-aware matching
                if key in tau_sig_lookup and len(tau_sig_lookup[key]) > 0:
                    # Pattern-aware: look for TAU_SIG with subscan = M51CENTER_subscan - 1
                    tau_sig_entries = tau_sig_lookup[key]
                    
                    # First try: exact pattern match (TAU_SIG subscan = M51CENTER subscan - 1)
                    pattern_matches = [entry for entry in tau_sig_entries if entry[0] == m51_subscan - 1]
                    
                    if pattern_matches:
                        # Use the first pattern match
                        tau_sig_indices[idx] = pattern_matches[0][1]
                        n_tau_sig_linked += 1
                        logger.debug(f"M51CENTER idx={idx} subscan={m51_subscan} -> TAU_SIG idx={pattern_matches[0][1]} subscan={pattern_matches[0][0]} (pattern match)")
                    else:
                        # Fallback: use nearest TAU_SIG by time if available, else use first TAU_SIG
                        if m51_timestamp is not None and any(entry[2] is not None for entry in tau_sig_entries):
                            # Find TAU_SIG with closest timestamp
                            valid_entries = [e for e in tau_sig_entries if e[2] is not None]
                            nearest = min(valid_entries, key=lambda e: abs(e[2] - m51_timestamp))
                            tau_sig_indices[idx] = nearest[1]
                            n_tau_sig_fallback += 1
                            logger.debug(f"M51CENTER idx={idx} subscan={m51_subscan} -> TAU_SIG idx={nearest[1]} subscan={nearest[0]} (temporal fallback)")
                        else:
                            # Last resort: use first TAU_SIG in group
                            tau_sig_indices[idx] = tau_sig_entries[0][1]
                            n_tau_sig_fallback += 1
                            logger.debug(f"M51CENTER idx={idx} subscan={m51_subscan} -> TAU_SIG idx={tau_sig_entries[0][1]} subscan={tau_sig_entries[0][0]} (first match fallback)")
        
        m51_total = np.sum(filtered_objects == object_filter)
        logger.info(f"Linked {n_tsys_linked}/{m51_total} M51CENTER spectra to TSYS (pattern-aware)")
        logger.info(f"Linked {n_tsys_fallback}/{m51_total} M51CENTER spectra to TSYS (fallback)")
        logger.info(f"Linked {n_tau_sig_linked}/{m51_total} M51CENTER spectra to TAU_SIG (pattern-aware)")
        logger.info(f"Linked {n_tau_sig_fallback}/{m51_total} M51CENTER spectra to TAU_SIG (fallback)")
        
        # Add the index columns to filtered_data by extending the table
        # Build column list from filtered data arrays
        col_list = []
        for name in filtered_data.columns.names:
            col_array = filtered_data[name]  # Get sliced data
            col_format = filtered_data.columns[name].format
            col_list.append(fits.Column(name=name, format=col_format, array=col_array))
        
        # Add TSYS_INDEX and TAU_SIG_INDEX columns as simple int32 arrays
        col_list.append(fits.Column(name='TSYS_INDEX', format='J', array=tsys_indices.astype(np.int32)))
        col_list.append(fits.Column(name='TAU_SIG_INDEX', format='J', array=tau_sig_indices.astype(np.int32)))
        
        # Create new table with all columns
        # CRITICAL: Preserve all header keywords from original filtered_data
        # BinTableHDU.from_columns() creates a new table but may not preserve all keywords
        header_keywords_to_preserve = dict(header)
        filtered_data_with_indices = fits.BinTableHDU.from_columns(col_list, header=header)
        
        # Restore any keywords that might have been lost during table creation
        for key, value in header_keywords_to_preserve.items():
            if key not in filtered_data_with_indices.header:
                logger.debug(f"Restoring header keyword {key}")
                filtered_data_with_indices.header[key] = value
        
        # Create new HDU with filtered data
        primary_hdu = fits.PrimaryHDU()
        table_hdu = filtered_data_with_indices
        
        # Create output FITS
        logger.info(f"Writing output to {output_fits}...")
        os.makedirs(os.path.dirname(output_fits) or '.', exist_ok=True)
        
        hdul_out = fits.HDUList([primary_hdu, table_hdu])
        # Use overwrite=True and ensure checksum is not used to avoid encoding issues
        hdul_out.writeto(output_fits, overwrite=True, checksum=False)
        
        logger.info(f"✓ Wrote {n_filtered} spectra to {output_fits}")


def prepare_for_pca(fits_file: Optional[str] = None,
                   output_fits: Optional[str] = None,
                   config: Optional[str] = None,
                   pca_source: Optional[str] = None,
                   object_filter: Optional[str] = None,
                   mission_id: Optional[str] = None,
                   scan: Optional[int] = None,
                   fill_noise: bool = False) -> None:
    """
    Main entry point for prepare_for_pca functionality.
    
    Loads configuration from config.toml if provided, using defaults for missing parameters.
    
    Parameters
    ----------
    fits_file : str, optional
        Input FITS file (overrides config)
    output_fits : str, optional
        Output FITS file (overrides config)
    config : str, optional
        Path to config.toml file
    pca_source : str, optional
        PCA source to filter for (overrides config)
    object_filter : str, optional
        Object substring to filter for (overrides config)
    mission_id : str, optional
        Filter to specific mission_id (e.g., "2017-02-01_GR_F367")
    scan : int, optional
        Filter to specific SCAN number
    """
    # Configure logging to display output
    if not logger.handlers:
        handler = logging.StreamHandler()
        handler.setLevel(logging.INFO)
        formatter = logging.Formatter('%(message)s')
        handler.setFormatter(formatter)
        logger.addHandler(handler)
        logger.setLevel(logging.INFO)
    
    # Load configuration
    if config or (not fits_file and not output_fits):
        try:
            import tomllib
        except ModuleNotFoundError:
            import tomli as tomllib
        
        config_path = config or "config.toml"
        logger.info(f"Loading configuration from {config_path}")
        
        with open(config_path, 'rb') as f:
            cfg = tomllib.load(f)
        
        # Get input file
        if not fits_file:
            output_cfg = cfg.get('output', {})
            if 'reduced_fits' in output_cfg:
                fits_file = output_cfg['reduced_fits']
                logger.info(f"Using reduced_fits from config: {fits_file}")
            else:
                raise ValueError("Could not find input file. Specify --fits or ensure output.reduced_fits in config.toml")
        
        # Get output file
        if not output_fits:
            output_cfg = cfg.get('output', {})
            if 'prepared_for_pca' in output_cfg:
                output_fits = output_cfg['prepared_for_pca']
                logger.info(f"Using prepared_for_pca from config: {output_fits}")
            else:
                raise ValueError("Could not find output file. Specify --output or ensure output.prepared_for_pca in config.toml")
        
        # Get PCA source
        if not pca_source:
            pca_cfg = cfg.get('pca', {})
            pca_source = pca_cfg.get('pca_source', 'SKYCHOPDIFF')
        
        # Get object filter
        if not object_filter:
            parameters_cfg = cfg.get('parameters', {})
            object_filter = parameters_cfg.get('object', 'M51CENTER')
    
    # Validate required parameters
    if not fits_file:
        raise ValueError("Input FITS file not specified")
    if not output_fits:
        raise ValueError("Output FITS file not specified")
    if not pca_source:
        pca_source = 'SKYCHOPDIFF'
    if not object_filter:
        object_filter = 'M51CENTER'
    
    logger.info(f"PCA source filter: {pca_source}")
    logger.info(f"Object filter: {object_filter}")
    if mission_id:
        logger.info(f"Mission ID filter: {mission_id}")
    if scan is not None:
        logger.info(f"Scan filter: {scan}")
    
    # Process the file
    fill_telluric_with_noise(fits_file, output_fits, pca_source, object_filter,
                            mission_id=mission_id, scan=scan, fill_noise=fill_noise)
