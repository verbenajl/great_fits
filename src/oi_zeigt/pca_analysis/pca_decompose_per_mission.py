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


def load_spectra_by_mission(fits_file: str, mission_params: dict = None,
                            pca_source: str = 'SKYCHOPDIFF') -> Dict[str, Dict]:
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
            
            # Filter for pca_source OBJECT value if OBJECT column exists
            if 'OBJECT' in data.dtype.names:
                objects = np.array([str(x).strip() if isinstance(x, bytes) else str(x).strip()
                                   for x in mission_spectra['OBJECT']])
                sky_mask = objects == pca_source
                sky_spectra = mission_spectra[sky_mask]['SPECTRUM']
                # Get the original indices in the full data array
                mission_indices = np.where(mask)[0]
                sky_indices = mission_indices[sky_mask]
                logger.info(f"    └─ {pca_source}: {np.sum(sky_mask)} spectra")
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
                    'source': pca_source if 'OBJECT' in data.dtype.names else 'all',
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
    first_valid : int
        Index of the first valid channel in the original spectra array.
    last_valid : int
        Index of the last valid channel in the original spectra array.
    """
    # Filter out spectra that are mostly NaN (>50% NaN)
    nan_fraction_per_spectrum = np.sum(np.isnan(spectra), axis=1) / spectra.shape[1]
    good_spectra_mask = nan_fraction_per_spectrum < 0.20
    original_count = len(spectra)
    spectra = spectra[good_spectra_mask]
    removed_count = original_count - len(spectra)
    if removed_count > 0:
        logger.info(f"Filtered out {removed_count} spectra with >20% NaN values")
    
    if len(spectra) == 0:
        logger.error("No spectra remain after NaN filtering!")
        raise ValueError("All spectra have >50% NaN values")
    
    # Find valid channels (not all NaN)
    valid_channels = ~np.all(np.isnan(spectra), axis=0)
    if not np.any(valid_channels):
        logger.error("No valid channels found!")
        raise ValueError("All channels are NaN in all spectra")
    
    original_n_channels = spectra.shape[1]
    first_valid = int(np.where(valid_channels)[0][0])
    last_valid = int(np.where(valid_channels)[0][-1])
    spectra = spectra[:, first_valid:last_valid+1]
    if velocity_axis is not None and len(velocity_axis) == original_n_channels:
        velocity_axis = velocity_axis[first_valid:last_valid + 1]

    n_valid_channels = last_valid - first_valid + 1
    logger.info(f"Trimmed to {n_valid_channels} valid channels (out of {original_n_channels} original), "
                f"channel offset: [{first_valid}, {last_valid}]")
    
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

    return spectra, velocity_axis, first_valid, last_valid


def decompose_mission_spectra(mission_id: str, telescope: str, spectra: np.ndarray,
                              flight_date: str, n_components: int = 5,
                              velocity_axis: np.ndarray = None,
                              line_window_kms: tuple = None,
                              science_line_window_kms: tuple = None,
                              fit_exclude_window_kms: tuple = None,
                              spectrum_indices: np.ndarray = None,
                              smoothing_kernel_size: int = None,
                              reject_lined_rows: bool = False,
                              reject_line_sigma: float = 3.0,
                              pca_source: str = 'SKYCHOPDIFF') -> DecompositionResult:
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
        (v_min, v_max) in km/s of the atmospheric telluric line to exclude from
        the baseline fit and the PCA fit during decomposition
    science_line_window_kms : tuple, optional
        (v_min, v_max) in km/s of the [CII] science line ([reduction].line_window).
        Always held out of the baseline polynomial (a linear baseline shouldn't be
        dragged by a possible line). It is NOT, by itself, held out of the PCA fit
        — that is governed separately by fit_exclude_window_kms.
    fit_exclude_window_kms : tuple, optional
        (v_min, v_max) in km/s of the window to also hold OUT of the PCA component
        fit. None (default) = the science line is fully INCLUDED in the fit, so the
        components get a real baseline across it. Set to a NARROW core around the
        line peak (via --exclude-narrow-line-window / [pca].narrow_line_window) to
        keep the destructive peak out of the basis while keeping the broad wings as
        baseline; or to the full science window to reproduce the old leak-fix. The
        telluric window is excluded from the fit regardless.
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
    spectra, velocity_axis, channel_first, channel_last = preprocess_spectra(spectra, velocity_axis)
    logger.info(f"After preprocessing: {spectra.shape}")
    
    # Build channel-exclusion masks. The BASELINE polynomial and the PCA fit are
    # masked differently:
    #   * telluric window (line_window_kms): always excluded from BOTH — it is
    #     atmospheric junk that must never enter either fit.
    #   * science line window (science_line_window_kms, from [reduction].line_window):
    #     always excluded from the BASELINE polynomial (a linear baseline should
    #     not be dragged by a possible line).
    #   * fit_exclude_window_kms: the window (if any) additionally held out of the
    #     PCA FIT. None (default) → the science line is fully IN the fit, so the
    #     components get a real baseline across it. A NARROW core around the line
    #     peak keeps the destructive peak out of the basis (so pca_correct can't
    #     over-subtract the bright central line) while the broad wings still supply
    #     the baseline. The full science window here reproduces the old leak-fix.
    # The windows are generally disjoint, so we build explicit boolean channel
    # masks rather than a single (v_min,v_max) range.
    def _window_mask(win, label, reason):
        if win is None or velocity_axis is None:
            return np.zeros(spectra.shape[1], dtype=bool)
        w_lo, w_hi = win
        m = (velocity_axis >= w_lo) & (velocity_axis <= w_hi)
        if np.any(m):
            logger.info(f"  {reason} {label} window {w_lo:.1f}-{w_hi:.1f} km/s "
                        f"({int(np.sum(m))} channels)")
        return m

    telluric_win_mask = _window_mask(line_window_kms, 'telluric',
                                     'Excluding (baseline+fit)')
    science_win_mask = _window_mask(science_line_window_kms, 'science',
                                    'Excluding from baseline')
    fit_win_mask = _window_mask(fit_exclude_window_kms, 'science-core',
                                'Excluding from PCA fit')
    if fit_exclude_window_kms is None:
        logger.info("  Science line window INCLUDED in the PCA fit "
                    "(no --exclude-narrow-line-window)")

    baseline_exclude_mask = telluric_win_mask | science_win_mask
    baseline_exclude_mask = baseline_exclude_mask if np.any(baseline_exclude_mask) else None

    pca_exclude_mask = telluric_win_mask | fit_win_mask
    pca_exclude_mask = pca_exclude_mask if np.any(pca_exclude_mask) else None

    # Baseline-subtract each SKYCHOPDIFF spectrum before decomposition, matching
    # the original pyclass behaviour where prepare_spectrum() calls baseline()
    # prior to adding a spectrum to the PCA input set. Both line windows are held
    # out of the baseline fit so neither telluric emission nor a (possibly
    # contaminating) [CII] line biases the polynomial.
    from oi_zeigt.reduction.core import baseline_subtract
    baselined = np.empty_like(spectra)
    for s_idx in range(len(spectra)):
        baselined[s_idx] = baseline_subtract(spectra[s_idx], order=1,
                                             exclude_mask=baseline_exclude_mask)
    spectra = baselined
    logger.info(f"  ✓ Baseline-subtracted {len(spectra)} spectra (order=1, "
                f"excluded {0 if baseline_exclude_mask is None else int(np.sum(baseline_exclude_mask))} channels)")

    # ── Optional per-ROW line rejection ───────────────────────────────────────
    # Drop the individual SKYCHOPDIFF reference spectra that themselves carry a
    # [CII] line, so the PCA basis is trained ONLY on clean (line-free) rows and
    # the components never learn the leak direction. This is a per-ROW filter,
    # complementary to (and stackable with) the per-CHANNEL science-core holdout
    # (fit_exclude_window_kms / --exclude-narrow-line-window) and the per-FLIGHT
    # skip_correction pass-through. Detection is a data-driven matched filter:
    # the median line-window profile across the rows (the coherent leak shape) is
    # the template, each row is projected onto it and divided by its own line-free
    # noise, and rows with |MF/σ| > reject_line_sigma are removed before the fit.
    # NOTE: like preprocess_spectra's NaN drop above, this decouples the fitted
    # row count from spectrum_indices (metadata only), which is already the case.
    if reject_lined_rows and np.any(science_win_mask):
        line_ch = science_win_mask
        free_ch = ~(telluric_win_mask | science_win_mask)
        n_before = len(spectra)
        if np.sum(free_ch) >= 5 and np.sum(line_ch) >= 3 and n_before > 0:
            sigma = np.std(spectra[:, free_ch], axis=1)
            templ = np.median(spectra[:, line_ch], axis=0)
            tnorm = float(np.linalg.norm(templ))
            if tnorm > 0 and np.all(np.isfinite(templ)):
                templ = templ / tnorm
                mf = spectra[:, line_ch] @ templ
                with np.errstate(divide='ignore', invalid='ignore'):
                    row_snr = np.where(sigma > 0, np.abs(mf) / sigma, 0.0)
                keep = row_snr <= reject_line_sigma
                n_keep = int(np.sum(keep))
                min_keep = max(n_components + 2, 10)
                if n_keep < min_keep:
                    logger.warning(
                        f"  ⚠ Per-row line rejection would leave only {n_keep}/"
                        f"{n_before} rows (< {min_keep} needed for {n_components} "
                        f"components) — KEEPING ALL rows for {mission_id}")
                else:
                    spectra = spectra[keep]
                    frac = 100.0 * (n_before - n_keep) / n_before
                    logger.info(
                        f"  ✓ Per-row line rejection (|MF/σ| > {reject_line_sigma}): "
                        f"dropped {n_before - n_keep}/{n_before} lined SKYCHOPDIFF "
                        f"rows ({frac:.1f}%); fitting on {n_keep} clean rows")
            else:
                logger.info("  Per-row line rejection: degenerate template — skipped")
        else:
            logger.info("  Per-row line rejection: too few line/free channels — skipped")

    # Exclude the chosen channels from the PCA fit entirely, instead of injecting
    # a mean-filled placeholder. Matches legacy pyclass, which physically trims
    # the channel range before fitting (apply_pca_exclude_range) so the components
    # never encode anything about that region. Implemented as a boolean mask
    # rather than a physical trim so the rest of the pipeline (correction,
    # plotting) keeps working with full-length, fixed-shape arrays — the excluded
    # channels are interpolated back in below, after fitting/smoothing.
    telluric_mask = None
    if pca_exclude_mask is not None:
        n_excl = int(np.sum(pca_exclude_mask))
        if n_excl >= spectra.shape[1]:
            logger.warning(f"  ⚠ Exclusion windows cover the entire spectrum — ignoring exclusion")
        else:
            telluric_mask = pca_exclude_mask
            _what = 'telluric + science-core' if fit_exclude_window_kms is not None else 'telluric only'
            logger.info(f"  ✓ Excluding {n_excl} channels ({_what}) from PCA fit "
                        f"(never mean-filled or otherwise faked)")

    # Fit PCA on the non-excluded channels only
    fit_spectra = spectra[:, ~telluric_mask] if telluric_mask is not None else spectra
    logger.info(f"Fitting PCA with {n_components} components on {fit_spectra.shape[1]}/{spectra.shape[1]} channels...")
    decomposer = PCADecomposer(n_components=n_components, scale=False)
    decomposer.fit(fit_spectra)

    # Smooth PCA components with a boxcar kernel (matches legacy smooth_pca_components()
    # in pca_decompose.py, which used np.convolve with mode='same'). Smoothing runs on
    # the fitted (possibly telluric-excluded) array, so channels immediately either
    # side of an excluded window become smoothing neighbours of each other — the same
    # "stitching" legacy gets from smoothing/baselining a physically trimmed array.
    if smoothing_kernel_size:
        kernel_size = int(smoothing_kernel_size)
        kernel = np.ones(kernel_size) / kernel_size
        # np.convolve(mode='same') returns length max(len(signal), len(kernel)),
        # so when the kernel is LONGER than the fitted-channel array (which happens
        # here once wide science+telluric exclusion windows leave fewer channels
        # than smoothing_kernel_size) it silently grows the component by
        # (kernel_size - n_channels) samples and the later interp back to full
        # length fails ("fp and xp are not of the same length"). Compute the full
        # convolution and take the centred slice so the output is ALWAYS the same
        # length as the input component, for any kernel size.
        n_ch_fit = decomposer.pca_model.components_.shape[1]
        _off = (kernel_size - 1) // 2
        smoothed = []
        for component in decomposer.pca_model.components_:
            full = np.convolve(component.copy(), kernel, mode='full')
            smoothed.append(full[_off:_off + n_ch_fit])
        decomposer.pca_model.components_ = np.array(smoothed)
        logger.info(f"Smoothed {len(smoothed)} components with boxcar kernel size {kernel_size}")
        if kernel_size > n_ch_fit:
            logger.warning(f"  ⚠ smoothing_kernel_size ({kernel_size}) exceeds the "
                           f"{n_ch_fit} channels left after exclusion — smoothing is very "
                           f"aggressive (near-flat components); consider lowering it.")

    components_fit = decomposer.get_components()   # (n_comp, n_channels_fit)
    mean_fit = decomposer.pca_model.mean_           # (n_channels_fit,)

    # Pad components/mean back to the full channel length. The excluded windows
    # were never part of the fit, so PCA reports no eigenvector value there. We
    # fill them by LINEAR INTERPOLATION from the flanking fitted channels rather
    # than zero-filling. Rationale (see AskUserQuestion 2026-07-13, "off-off"
    # model): the SKYCHOPDIFF references are off-off spectra assumed free of the
    # [CII] science line, so the component's true behaviour across the window is
    # its smooth baseline continuation, NOT a line and NOT a hard-zero plateau.
    #   * Interpolation reads only from good-channel flanks, which never saw the
    #     weak leaked line (F296 ~+0.6 K), so it recovers the baseline level, not
    #     the line — the line stays mathematically untouched by correction.
    #   * Unlike zero-filling (which flatlined ~47% of the component and removed
    #     its ability to model the baseline *under* the line), the interpolated
    #     component subtracts the continuum beneath the line while leaving the
    #     narrow line itself intact.
    # np.interp clamps to the nearest endpoint for any excluded channel outside
    # the good-channel range (constant extrapolation at spectrum edges).
    #
    # CAVEAT that the science-core block below fixes: interpolation across a WIDE
    # window draws a straight CHORD between the two flanks. For the ~7-ch telluric
    # window the flanks sit at nearly the same level, so the chord is ~flat and
    # harmless. But the science core ([pca].narrow_line_window, ~25 ch) spans a
    # region where a real eigenvector genuinely differs between flanks, so the
    # chord is a SLOPE — and pca_correct subtracts coeff*component across the line,
    # injecting a fake line (one coeff sign) or a trough (the other) right where we
    # are most sensitive. So we interpolate everything (good for telluric) and then
    # force the science core to ZERO (below).
    if telluric_mask is not None:
        n_channels_full = spectra.shape[1]
        good = ~telluric_mask
        x = np.arange(n_channels_full)
        xg = x[good]
        components_full = np.empty((components_fit.shape[0], n_channels_full), dtype=components_fit.dtype)
        for _c in range(components_fit.shape[0]):
            components_full[_c] = np.interp(x, xg, components_fit[_c])
        mean_full = np.interp(x, xg, mean_fit).astype(mean_fit.dtype)
        logger.info(f"  ✓ Interpolated components/mean across {int(np.sum(telluric_mask))} "
                    f"excluded channels (smooth baseline continuation, no zero plateau)")

        # SCIENCE-CORE PROTECTION (only when --exclude-narrow-line-window is on).
        # Force the components to ZERO inside the science core, cosine-tapered from
        # the full flank value down to zero over a few channels at each edge so
        # there is no hard step. Then coeff*component == 0 across the protected
        # interior: the PCA correction is a strict NO-OP there and can neither add
        # a line nor dig a trough — the [CII] line is left exactly as the order-1
        # baseline poly leaves it. The taper keeps real correction at the core
        # edges (line wings, adjacent to fitted channels) while fully protecting
        # the centre (the line peak). Only the COMPONENTS are zeroed: the mean just
        # centres the coefficient fit and is never subtracted from the data
        # (pca_correct subtracts coeff*component only — pca_correct_fits.py:566,
        # 1503), so a zero component alone fully protects the core. Supersedes the
        # 2026-07-13 "interpolate, don't zero" choice, which assumed similar flanks.
        if fit_win_mask is not None and np.any(fit_win_mask):
            idx = np.where(fit_win_mask)[0]
            lo, hi = int(idx.min()), int(idx.max())
            width = hi - lo + 1
            taper = np.zeros(n_channels_full, dtype=components_full.dtype)
            taper[:lo] = 1.0
            taper[hi + 1:] = 1.0
            ramp = int(min(4, max(1, width // 4)))   # cosine-ramp channels per edge
            for k in range(ramp):
                w = 0.5 * (1.0 + np.cos(np.pi * k / ramp))   # 1 at boundary → 0 inward
                taper[lo + k] = w
                taper[hi - k] = w
            components_full = components_full * taper[None, :]
            logger.info(f"  ✓ Zeroed components across the {width}-channel science core "
                        f"(cosine-tapered {ramp} ch/edge) — correction is a no-op there, "
                        f"[CII] line left intact")
    else:
        components_full = components_fit
        mean_full = mean_fit

    # Create result
    if spectrum_indices is not None:
        logger.info(f"Creating DecompositionResult with {len(spectrum_indices)} spectrum indices")
    else:
        logger.info(f"Creating DecompositionResult with NO spectrum indices (spectrum_indices is None)")

    result = DecompositionResult(
        mean_spectrum=mean_full,
        components=components_full,
        explained_variance=decomposer.get_explained_variance(),
        explained_variance_ratio=decomposer.get_explained_variance_ratio(),
        config=DecompositionConfig(n_components=n_components),
        metadata={
            'n_reference_spectra': len(spectra),
            'n_channels': spectra.shape[1],
            'channel_first': channel_first,
            'channel_last': channel_last,
            'source': pca_source,
            'mission_id': mission_id,
            'velocity_axis_kms': velocity_axis,
            'telluric_mask': telluric_mask,
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


def save_results(result: DecompositionResult, mission_id: str, telescop: str,
                 output_dir: Path = None) -> Path:
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
    output_dir : Path, optional
        Directory to write output files. Defaults to output/pca_components.

    Returns
    -------
    Path
        Path to saved pickle file
    """
    if output_dir is None:
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
        help="Limit decomposition to MISSION_IDs containing this substring (e.g. F373). "
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
        "--filter-flight",
        type=str,
        nargs='+',
        default=None,
        dest='filter_flight',
        help="Exclude all entries whose MISSION_ID contains these strings "
             "(space-separated, e.g. F528 F299)."
    )
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Verbose output"
    )
    parser.add_argument(
        "--plot",
        action="store_true",
        dest="plot_components",
        help="Generate visualization plots of PCA components"
    )
    parser.add_argument(
        "--exclude-narrow-line-window",
        action="store_true",
        dest="exclude_narrow_line",
        help="Hold a NARROW core around the [CII] line peak out of the PCA "
             "component fit (window from [pca].narrow_line_window, else the "
             "central 40%% of [reduction].line_window). The broad line wings "
             "stay in the fit as baseline, but the destructive peak never enters "
             "the basis — prevents pca_correct over-subtracting the bright "
             "central line. Without this flag the whole line window is included "
             "in the fit."
    )
    parser.add_argument(
        "--reject-lined-references",
        action="store_true",
        dest="reject_lined_references",
        help="Drop individual SKYCHOPDIFF reference spectra that themselves "
             "carry a [CII] line (per-ROW matched-filter detection over the "
             "science line window) BEFORE the component fit, so the basis is "
             "trained only on clean reference rows. Complementary to "
             "--exclude-narrow-line-window (per-channel) and skip_correction "
             "(per-flight). Enable via [pca].reject_lined_references too."
    )
    parser.add_argument(
        "--reject-line-sigma",
        type=float,
        default=None,
        dest="reject_line_sigma",
        help="Threshold |matched-filter/σ| above which a SKYCHOPDIFF row is "
             "treated as lined and dropped (default 3.0, or [pca].reject_line_sigma). "
             "Per-mission override via the YAML key reject_line_sigma. Lower = "
             "more aggressive rejection."
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
        
        # Resolve components_dir: config [output][components_dir] > default
        output_config = config.get('output', {})
        components_dir = Path(output_config.get('components_dir', 'output/pca_components'))
        logger.info(f"✓ Components directory = {components_dir}")

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
        pca_source = pca_config.get('pca_source', 'SKYCHOPDIFF')
        logger.info(f"✓ PCA source = {pca_source} (from config [pca][pca_source])" if pca_config.get('pca_source') else f"✓ PCA source = {pca_source} (default)")
        mission_data = load_spectra_by_mission(fits_file, pca_source=pca_source)
        
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
        
        # The [CII] science line window ([reduction].line_window). It is always
        # held out of the baseline polynomial; whether a (narrow) core is also
        # held out of the PCA fit is governed by --exclude-narrow-line-window
        # (default: whole line window IS included in the component fit). See
        # decompose_mission_spectra for the trade-off.
        _rw = config.get('reduction', {}).get('line_window')
        science_line_window_kms = (float(_rw[0]), float(_rw[1])) if _rw and len(_rw) == 2 else None

        # Decide which window (if any) is also held out of the PCA FIT. Default:
        # None → the science line is fully IN the fit. With --exclude-narrow-line-window
        # we hold out only a NARROW core around the line peak, taken from
        # [pca].narrow_line_window; if that key is absent we fall back to the
        # central 40% of the wide science window so the flag works out-of-the-box.
        fit_exclude_window_kms = None
        if args.exclude_narrow_line:
            _nw = pca_config.get('narrow_line_window')
            if _nw and len(_nw) == 2:
                fit_exclude_window_kms = (float(_nw[0]), float(_nw[1]))
                _src = '[pca].narrow_line_window'
            elif science_line_window_kms is not None:
                _lo, _hi = science_line_window_kms
                _span = _hi - _lo
                fit_exclude_window_kms = (_lo + 0.3 * _span, _hi - 0.3 * _span)
                _src = 'central 40% of [reduction].line_window (no [pca].narrow_line_window set)'
            else:
                logger.warning("  ⚠ --exclude-narrow-line-window set but no "
                               "narrow_line_window and no line_window — ignoring")
                _src = None
            if fit_exclude_window_kms is not None:
                logger.info(f"✓ --exclude-narrow-line-window: holding "
                            f"{fit_exclude_window_kms[0]:.1f}-{fit_exclude_window_kms[1]:.1f} "
                            f"km/s out of the PCA fit (from {_src})")

        if science_line_window_kms is not None:
            _fit_state = (f"core {fit_exclude_window_kms[0]:.1f}-{fit_exclude_window_kms[1]:.1f} "
                          f"km/s excluded from PCA fit" if fit_exclude_window_kms is not None
                          else "fully INCLUDED in PCA fit")
            logger.info(f"✓ Science line window ({science_line_window_kms[0]:.1f}-"
                        f"{science_line_window_kms[1]:.1f} km/s): held out of baseline, "
                        f"{_fit_state}")

        # Per-ROW line rejection: drop lined SKYCHOPDIFF reference spectra before
        # the component fit. Enabled by CLI flag OR [pca].reject_lined_references.
        # Base sigma from CLI > [pca].reject_line_sigma > 3.0 default; a per-mission
        # reject_line_sigma in the YAML can override it below (see loop).
        reject_lined_rows = bool(args.reject_lined_references
                                 or pca_config.get('reject_lined_references', False))
        if args.reject_line_sigma is not None:
            base_reject_sigma = float(args.reject_line_sigma)
        else:
            base_reject_sigma = float(pca_config.get('reject_line_sigma', 3.0))
        if reject_lined_rows:
            if science_line_window_kms is None:
                logger.warning("  ⚠ reject-lined-references set but no "
                               "[reduction].line_window — per-row rejection will be skipped")
            else:
                logger.info(f"✓ Per-row line rejection ENABLED (base |MF/σ| > "
                            f"{base_reject_sigma}); lined SKYCHOPDIFF rows dropped before fit")

        # Load mission-specific parameters from mission_id_parameters.yml
        mission_line_windows = {}
        try:
            # Check [pca] then [input] then fall back to bundled file.
            _mission_params_path = (pca_config.get('mission_parameters')
                                    or config.get('input', {}).get('mission_parameters'))
            if _mission_params_path:
                mission_yml = Path(_mission_params_path)
            else:
                mission_yml = Path(__file__).parent / "mission_id_parameters.yml"
            if mission_yml.exists():
                logger.info(f"Loading mission parameters from {mission_yml}")
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

        # Flight filter: remove all entries whose MISSION_ID contains any of the given strings
        if args.filter_flight:
            kept = {}
            removed_missions = set()
            for k, v in mission_data.items():
                mid = v['mission_id']
                if any(f in mid for f in args.filter_flight):
                    removed_missions.add(mid)
                else:
                    kept[k] = v
            for mid in sorted(removed_missions):
                logger.info(f"  Removed flight: {mid}")
            logger.info(
                f"Flight filter {args.filter_flight}: removed {len(mission_data) - len(kept)} "
                f"combinations ({len(mission_data)} → {len(kept)})"
            )
            mission_data = kept
            if not mission_data:
                logger.error("No data remaining after --filter-flight")
                sys.exit(1)

        # Filter to a single mission_id / telescope if requested
        if args.mission_id or args.telescope:
            before = len(mission_data)
            mission_data_full = mission_data  # keep for error reporting
            mission_data = {
                k: v for k, v in mission_data.items()
                if (args.mission_id is None or args.mission_id in v['mission_id'])
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

        # Check for missions with no pca_source spectra
        missions_no_sky = [m for m, d in mission_data.items() if d['metadata']['n_spectra'] == 0]
        if missions_no_sky:
            logger.warning(f"Missions with NO {pca_source} spectra: {', '.join(str(m) for m in missions_no_sky)}")

        logger.info(f"\n{'='*80}")
        logger.info(f"PROCESSING {len(mission_data)} MISSION/TELESCOPE COMBINATIONS")
        logger.info(f"{'='*80}\n")

        # Per-mission n_components: reuse the SAME per-mission YAML + key matching
        # that pca_correct --per-mission-parameters uses (config[input][mission_pca_parameters]).
        # Precedence for n_components: CLI --n-components  >  per-mission YAML  >
        # config [pca][n_components]  >  default. `n_components` computed above is
        # the base/fallback; a CLI value overrides the YAML globally.
        from oi_zeigt.pca_analysis.pca_correct_fits import (
            load_mission_pca_parameters, _resolve_mission_param)
        mission_pca_params = {}
        _pca_yaml = config.get('input', {}).get('mission_pca_parameters')
        if args.n_components is not None:
            logger.info(f"--n-components given on CLI ({n_components}); "
                        f"per-mission n_components in the YAML will be ignored")
        elif _pca_yaml:
            mission_pca_params = load_mission_pca_parameters(_pca_yaml)
            if mission_pca_params:
                logger.info(f"✓ Loaded per-mission PCA parameters from {_pca_yaml} "
                            f"({len(mission_pca_params)} missions) — honouring per-mission n_components")

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

            # Resolve per-mission n_components (CLI override already handled above:
            # when --n-components is given, mission_pca_params is empty so this
            # returns the base value).
            eff_n_components = int(_resolve_mission_param(
                mission_id, 'n_components', mission_pca_params, n_components))
            if eff_n_components != n_components:
                logger.info(f"  {mission_id}: n_components={eff_n_components} "
                            f"(per-mission YAML override; base {n_components})")

            # Per-mission reject_line_sigma override (falls back to base when absent)
            eff_reject_sigma = float(_resolve_mission_param(
                mission_id, 'reject_line_sigma', mission_pca_params, base_reject_sigma))
            if reject_lined_rows and eff_reject_sigma != base_reject_sigma:
                logger.info(f"  {mission_id}: reject_line_sigma={eff_reject_sigma} "
                            f"(per-mission YAML override; base {base_reject_sigma})")

            # Decompose with velocity axis and line masking if available
            result = decompose_mission_spectra(
                mission_id,
                telescop,
                data['spectra'],
                data['date'],
                n_components=eff_n_components,
                velocity_axis=velocity_axis,
                line_window_kms=mission_line_windows.get(mission_id),
                science_line_window_kms=science_line_window_kms,
                fit_exclude_window_kms=fit_exclude_window_kms,
                spectrum_indices=spectrum_indices,
                smoothing_kernel_size=pca_config.get('smoothing_kernel_size'),
                reject_lined_rows=reject_lined_rows,
                reject_line_sigma=eff_reject_sigma,
                pca_source=pca_source,
            )
            
            if result is None:
                continue
            
            results[key] = result
            
            # Save
            output_file = save_results(result, mission_id, telescop, output_dir=components_dir)
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
                            components_path = components_dir / f"pca_components_{mission_id}_{safe_telescop}.png"
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
                            variance_path = components_dir / f"pca_variance_{mission_id}_{safe_telescop}.png"
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
