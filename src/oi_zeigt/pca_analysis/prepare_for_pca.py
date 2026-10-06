"""
Prepare FITS data for PCA analysis.

This module filters FITS data to include only specific sources and objects,
and fills telluric line regions with Gaussian noise.
"""

import logging
import os
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Optional, Union

import numpy as np
import yaml
from astropy.io import fits

logger = logging.getLogger(__name__)


def build_tau_tsys_indices(data, object_filter: str):
    """
    Build TSYS_INDEX and TAU_SIG_INDEX arrays for a FITS BinTable data array.

    For each science spectrum (rows whose OBJECT contains *object_filter*) finds
    the best-matching TSYS / TAU_SIG row within the same (MISSION_ID, SCAN) group:
      1. Pattern match  : calibration SUBSCAN == science SUBSCAN - 1
      2. Temporal nearest: closest by UT or LST timestamp
      3. First-in-group : fallback of last resort

    Returns int32 arrays of length len(data); unmatched rows have value -1.
    """
    def _to_str(s):
        return s.decode().strip() if isinstance(s, bytes) else str(s).strip()

    objects     = np.array([_to_str(s) for s in data['OBJECT']])
    mission_ids = np.array([_to_str(s) for s in data['MISSION_ID']])

    tsys_indices    = np.full(len(data), -1, dtype=np.int32)
    tau_sig_indices = np.full(len(data), -1, dtype=np.int32)

    has_ut  = 'UT'  in data.dtype.names
    has_lst = 'LST' in data.dtype.names
    time_col = 'UT' if has_ut else ('LST' if has_lst else None)

    tsys_lookup    = {}
    tau_sig_lookup = {}

    for idx in range(len(data)):
        obj = objects[idx]
        if obj in ('TSYS', 'TAU_SIG'):
            key     = (mission_ids[idx], data['SCAN'][idx])
            subscan = data['SUBSCAN'][idx]
            timestamp = None
            if time_col:
                try:
                    timestamp = float(data[time_col][idx])
                except Exception:
                    pass
            entry = (subscan, idx, timestamp)
            if obj == 'TSYS':
                tsys_lookup.setdefault(key, []).append(entry)
            else:
                tau_sig_lookup.setdefault(key, []).append(entry)

    logger.info(f"Found {len(tsys_lookup)} TSYS groups and {len(tau_sig_lookup)} TAU_SIG groups")

    n_tsys_linked = n_tsys_fallback = n_tau_linked = n_tau_fallback = 0

    for idx in range(len(data)):
        if object_filter not in objects[idx]:
            continue
        key        = (mission_ids[idx], data['SCAN'][idx])
        sci_sub    = data['SUBSCAN'][idx]
        sci_ts     = None
        if time_col:
            try:
                sci_ts = float(data[time_col][idx])
            except Exception:
                pass

        def _best(lookup, idx_arr, n_linked, n_fallback):
            entries = lookup.get(key)
            if not entries:
                return n_linked, n_fallback
            pattern = [e for e in entries if e[0] == sci_sub - 1]
            if pattern:
                idx_arr[idx] = pattern[0][1]
                return n_linked + 1, n_fallback
            if sci_ts is not None and any(e[2] is not None for e in entries):
                valid = [e for e in entries if e[2] is not None]
                nearest = min(valid, key=lambda e: abs(e[2] - sci_ts))
                idx_arr[idx] = nearest[1]
            else:
                idx_arr[idx] = entries[0][1]
            return n_linked, n_fallback + 1

        n_tsys_linked,  n_tsys_fallback = _best(tsys_lookup,    tsys_indices,    n_tsys_linked,  n_tsys_fallback)
        n_tau_linked,   n_tau_fallback  = _best(tau_sig_lookup,  tau_sig_indices, n_tau_linked,   n_tau_fallback)

    sci_total = int(np.sum([object_filter in o for o in objects]))
    logger.info(f"Linked {n_tsys_linked}/{sci_total} science spectra to TSYS (pattern-aware)")
    logger.info(f"Linked {n_tsys_fallback}/{sci_total} science spectra to TSYS (fallback)")
    logger.info(f"Linked {n_tau_linked}/{sci_total} science spectra to TAU_SIG (pattern-aware)")
    logger.info(f"Linked {n_tau_fallback}/{sci_total} science spectra to TAU_SIG (fallback)")

    return tsys_indices, tau_sig_indices


def _default_mission_params_file() -> Path:
    """Return the bundled mission_id_parameters.yml path."""
    return Path(__file__).parent / 'mission_id_parameters.yml'


def load_mission_parameters(mission_id: str,
                            yaml_file: Optional[Union[str, Path]] = None) -> dict:
    """
    Load mission-specific parameters from YAML file.

    Parameters
    ----------
    mission_id : str
        Mission identifier
    yaml_file : str or Path, optional
        Path to the mission parameters YAML file.  Defaults to the
        ``mission_id_parameters.yml`` bundled with the package.

    Returns
    -------
    dict
        Mission parameters including telluric_line_center and telluric_line_width
    """
    yaml_path = Path(yaml_file) if yaml_file else _default_mission_params_file()

    if not yaml_path.exists():
        logger.warning(f"Mission parameters file not found: {yaml_path}")
        return {}

    try:
        with open(yaml_path, 'r') as f:
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
                            aor_id: Optional[str] = None,
                            fill_noise: bool = False,
                            filter_missions: bool = False,
                            filter_flight: Optional[list] = None,
                            filter_below: Optional[list] = None,
                            filter_above: Optional[list] = None,
                            mission_params_file: Optional[Union[str, Path]] = None,
                            baseline_order: Optional[int] = None,
                            mission_pca_params_file: Optional[Union[str, Path]] = None) -> None:
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
            from oi_zeigt.basic_io import reconstruct_velocity_axis
            velocity_axis_kms = reconstruct_velocity_axis(hdul[1]) / 1000.0
            logger.info(f"Velocity axis reconstructed: {len(velocity_axis_kms)} channels, "
                        f"{velocity_axis_kms[0]:.1f} to {velocity_axis_kms[-1]:.1f} km/s")
        except (AttributeError, KeyError, TypeError, ValueError) as e:
            logger.warning(f"Could not reconstruct velocity axis: {e}")
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
        
        # Keep spectra where OBJECT equals pca_source (exact) or contains object_filter
        # as a substring (to catch e.g. "M82" and "M82_P290" when object_filter="M82").
        # Also always include TSYS and TAU_SIG rows.
        combined_mask = (objects == pca_source) | \
                       np.array([object_filter in o for o in objects]) | \
                       (objects == 'TSYS') | (objects == 'TAU_SIG')
        
        # Apply mission_id filter if specified (substring match: F373 matches 2017-02-10_GR_F373)
        if mission_id:
            logger.info(f"Filtering to MISSION_ID containing '{mission_id}'...")
            mission_filter = np.array([mission_id in mid for mid in mission_ids])
            before_filter = np.sum(combined_mask)
            combined_mask = combined_mask & mission_filter
            after_filter = np.sum(combined_mask)
            matched = sorted(set(mission_ids[mission_filter].tolist()))
            logger.info(f"  Matched missions: {matched}")
            logger.info(f"  After mission_id filter: {after_filter}/{before_filter} rows kept")
        
        # Apply flight filter: remove all rows whose MISSION_ID contains any of the given strings
        if filter_flight:
            logger.info(f"Filtering out flights: {filter_flight}...")
            keep_mask = np.array([not any(f in mid for f in filter_flight) for mid in mission_ids])
            removed_missions = sorted(set(mission_ids[~keep_mask].tolist()))
            before_filter = int(np.sum(combined_mask))
            combined_mask = combined_mask & keep_mask
            after_filter = int(np.sum(combined_mask))
            for mid in removed_missions:
                logger.info(f"  Removed flight: {mid}")
            logger.info(f"  After flight filter: {after_filter}/{before_filter} rows kept")

        # Apply scan filter if specified
        if scan is not None:
            logger.info(f"Filtering to SCAN = {scan}...")
            scans = np.array(data['SCAN'])
            scan_filter = scans == scan
            before_filter = np.sum(combined_mask)
            combined_mask = combined_mask & scan_filter
            after_filter = np.sum(combined_mask)
            logger.info(f"  After scan filter: {after_filter}/{before_filter} rows kept")

        # Apply AOR_ID filter if specified
        if aor_id is not None and 'AOR_ID' in data.dtype.names:
            terms = [t.strip() for t in aor_id.split(',') if t.strip()]
            logger.info(f"Filtering to AOR_ID containing any of {terms}...")
            aor_id_col = np.array([
                s.decode().strip() if isinstance(s, bytes) else str(s).strip()
                for s in data['AOR_ID']
            ])
            aor_filter = np.array([any(t in a for t in terms) for a in aor_id_col])
            before_filter = int(np.sum(combined_mask))
            combined_mask = combined_mask & aor_filter
            after_filter = int(np.sum(combined_mask))
            matched = sorted(set(aor_id_col[aor_filter].tolist()))
            logger.info(f"  Matched AOR_IDs: {matched}")
            logger.info(f"  After AOR_ID filter: {after_filter}/{before_filter} rows kept")

        # Apply mission drop rules: the drop: blocks of mission_id_parameters.yml
        # (or the bundled default) AND of the per-mission PCA YAML, plus their
        # all_missions blocks — see oi_zeigt.drop_rules for the schema.
        if filter_missions:
            from oi_zeigt.drop_rules import load_drop_rules, drop_mask
            mission_yml = (Path(mission_params_file) if mission_params_file
                           else _default_mission_params_file())
            rules = load_drop_rules([mission_yml, mission_pca_params_file], echo=logger.info)
            n_before_drop = int(np.sum(combined_mask))
            combined_mask &= ~drop_mask(data, rules, within=combined_mask, echo=logger.info)
            n_after_drop = int(np.sum(combined_mask))
            logger.info(f"Mission drop rules: removed {n_before_drop - n_after_drop} rows "
                        f"({n_before_drop} → {n_after_drop})")

        # --- Parameter-value filters on the science + PCA-reference spectra -----
        # Mirrors filter_fits' --filter-below / --filter-above (implemented by
        # filter_and_save_fits._apply_param_filters): keep rows where COLUMN <= VALUE
        # (below) or COLUMN >= VALUE (above). A row with a non-finite value in COLUMN
        # always passes.
        #
        # SCOPE: the cut applies to the science rows (OBJECT contains object_filter)
        # AND the PCA reference rows (pca_source, e.g. SKYCHOPDIFF) — a noisy
        # reference pollutes the PCA basis, so it should be droppable too. TSYS and
        # TAU_SIG calibration rows are ALWAYS exempt: their RMS/RMSRATIO live on a
        # totally different scale (TSYS RMSRATIOB ~ hundreds) and a naive cut would
        # wipe out the calibration set.
        #
        # NB on columns: SKYCHOPDIFF has RMS_THEORETICAL = 0, so its RMSRATIO/RMSRATIOB
        # is NaN and a `--filter-below RMSRATIO 1.3` is a no-op on the reference (NaN
        # passes) — it only cuts science. To actually drop noisy references, filter on
        # a column that is finite for them, e.g. `--filter-below RMS <value>`.
        # As a convenience (and unlike the raise in filter_fits) RMSRATIO<->RMSRATIOB
        # fall back to each other, because the noise-filtered nopca inputs this
        # pipeline consumes carry RMSRATIOB while other products carry RMSRATIO.
        param_value_filters = ([(c, 'below', v) for (c, v) in (filter_below or [])] +
                               [(c, 'above', v) for (c, v) in (filter_above or [])])
        if param_value_filters:
            # Everything that is not calibration is subject to the cut.
            filterable_mask = (objects != 'TSYS') & (objects != 'TAU_SIG')
            for col, direction, value in param_value_filters:
                use_col = col
                if use_col not in data.dtype.names:
                    alt = {'RMSRATIO': 'RMSRATIOB', 'RMSRATIOB': 'RMSRATIO'}.get(col)
                    if alt and alt in data.dtype.names:
                        logger.warning(f"  Filter column '{col}' not found — "
                                       f"falling back to '{alt}'")
                        use_col = alt
                    else:
                        raise ValueError(f"Filter column '{col}' not found in FITS data")
                col_vals = np.asarray(data[use_col], dtype=np.float64)
                finite = np.isfinite(col_vals)
                if direction == 'below':
                    cond = (~finite) | (col_vals <= value)
                else:
                    cond = (~finite) | (col_vals >= value)
                # Non-calibration rows are subject to the cut; calibration rows
                # always pass (cond OR not-filterable).
                subject = combined_mask & filterable_mask
                sci = subject & np.array([object_filter in o for o in objects])
                ref = subject & (objects == pca_source)
                n_before = int(np.sum(combined_mask))
                removed = int(np.sum(subject & ~cond))
                removed_sci = int(np.sum(sci & ~cond))
                removed_ref = int(np.sum(ref & ~cond))
                combined_mask = combined_mask & (cond | ~filterable_mask)
                op = '<=' if direction == 'below' else '>='
                logger.info(f"Parameter filter: keep {use_col} {op} {value} "
                            f"(science + {pca_source}; TSYS/TAU_SIG exempt) — removed "
                            f"{removed} rows [{removed_sci} science, {removed_ref} {pca_source}] "
                            f"({n_before} → {int(np.sum(combined_mask))} rows)")

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
            _yml_used = (Path(mission_params_file) if mission_params_file
                         else _default_mission_params_file())
            logger.info(f"  Mission parameters file: {_yml_used}")

            # Pre-cache mission parameters — avoids re-parsing the YAML file for every spectrum
            unique_mission_ids = list(set(filtered_mission_ids))
            mission_params_cache = {
                mid: load_mission_parameters(mid, yaml_file=mission_params_file)
                for mid in unique_mission_ids
            }
            logger.info(f"  Cached parameters for {len(mission_params_cache)} unique mission IDs")

            # Build per-spectrum velocity axes as a plain numpy array for parallel workers
            spectra_array = np.array(filtered_data['SPECTRUM'], dtype=np.float64)
            n_spectra, n_chan = spectra_array.shape

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

        # --- Optional polynomial baseline subtraction on the science spectra ---
        # The order is supplied EXCLUSIVELY on the command line (--baseline N):
        # there is no default and nothing is read from config.toml — when
        # baseline_order is None this step is skipped entirely. Only science
        # rows (OBJECT == object_filter) are baselined; TSYS/TAU_SIG and the
        # PCA reference (pca_source, e.g. SKYCHOPDIFF) are left untouched — the
        # reference is separately baselined inside pca_decompose. The emission
        # line is protected by baseline_subtract's iterative sigma-clipping
        # (no fixed line window is used here, by design — see --baseline docs).
        if baseline_order is not None:
            from oi_zeigt.reduction.core import baseline_subtract
            order = int(baseline_order)
            spectra_bl = np.array(filtered_data['SPECTRUM'], dtype=np.float64)
            is_science = np.array([obj not in ('TSYS', 'TAU_SIG')
                                   and obj != pca_source
                                   for obj in filtered_objects])
            for i in np.where(is_science)[0]:
                spectra_bl[i] = baseline_subtract(spectra_bl[i], order=order)
            filtered_data['SPECTRUM'][:] = spectra_bl.astype(filtered_data['SPECTRUM'].dtype)
            logger.info(f"Baseline-subtracted (order {order}) {int(np.sum(is_science))} science "
                        f"spectra [CLI --baseline; line protected via sigma-clipping]")
        else:
            logger.info("Skipping baseline subtraction (--baseline not set)")

        # Create index columns for linking science spectra to TSYS and TAU_SIG
        logger.info("Creating TSYS and TAU_SIG index columns...")
        tsys_indices, tau_sig_indices = build_tau_tsys_indices(filtered_data, object_filter)
        
        # Add the index columns to filtered_data by extending the table
        # Build column list from filtered data arrays. Skip any pre-existing
        # TSYS_INDEX / TAU_SIG_INDEX columns: inputs that already went through
        # filter_fits --filter-tau (e.g. the nopca post_cleaned files this
        # pipeline consumes) carry these columns already, and re-appending them
        # below would create duplicate names ("name already used as a name or
        # title"). We recompute them fresh here, so drop the old copies.
        _index_cols = {'TSYS_INDEX', 'TAU_SIG_INDEX'}
        col_list = []
        for name in filtered_data.columns.names:
            if name in _index_cols:
                continue
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
                   aor_id: Optional[str] = None,
                   fill_noise: bool = False,
                   filter_missions: bool = False,
                   filter_flight: Optional[list] = None,
                   filter_below: Optional[list] = None,
                   filter_above: Optional[list] = None,
                   mission_params_file: Optional[str] = None,
                   baseline_order: Optional[int] = None,
                   mission_pca_params_file: Optional[str] = None) -> None:
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

        # Get mission parameters file (optional).
        # Check [pca] then [input] so both config styles work.
        if not mission_params_file:
            mission_params_file = (cfg.get('pca', {}).get('mission_parameters')
                                   or cfg.get('input', {}).get('mission_parameters'))

        # Per-mission PCA YAML: --filter-missions also honours its drop: blocks.
        if not mission_pca_params_file:
            mission_pca_params_file = cfg.get('input', {}).get('mission_pca_parameters')
    
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
    if mission_params_file:
        logger.info(f"Using mission parameters file: {mission_params_file}")
    fill_telluric_with_noise(fits_file, output_fits, pca_source, object_filter,
                            mission_id=mission_id, scan=scan, aor_id=aor_id,
                            fill_noise=fill_noise, filter_missions=filter_missions,
                            filter_flight=filter_flight,
                            filter_below=filter_below, filter_above=filter_above,
                            mission_params_file=mission_params_file,
                            baseline_order=baseline_order,
                            mission_pca_params_file=mission_pca_params_file)
