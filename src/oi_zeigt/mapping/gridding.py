"""
Core gridding and mapping functions for spectral data using WCS and cygrid.

Provides functionality to create proper WCS-based spatial maps from spectral 
observations using cygrid for optimal gridding with Gaussian kernel weighting.
Falls back to scipy griddata if cygrid is not available.
"""

from typing import Optional, Tuple, Dict, Any, List
import os
from pathlib import Path
import numpy as np
from astropy.io import fits
from astropy import units as u
from astropy.wcs import WCS
from astropy.table import Table
import matplotlib.pyplot as plt
import warnings

# Register GILDAS LUTs (e.g. 'rainbow3') so they are selectable via --colormap.
from ..basic_io import table_crpix1
from .gildas_luts import register_gildas_luts
register_gildas_luts()

# Try to import cygrid for optimal gridding
HAS_CYGRID = False
try:
    import cygrid
    HAS_CYGRID = True
except (ImportError, ValueError) as e:
    import logging as _logging
    _logging.getLogger(__name__).warning(f"cygrid import failed: {e}")


def get_gridding_params_from_config(config_path: Optional[str] = None,
                                    beamsize_deg: Optional[float] = None,
                                    pixsize_deg: Optional[float] = None) -> Tuple[float, float]:
    """
    Load gridding parameters from config.toml or use provided values.
    
    Reads beamsize_arcsec and pixel_size_arcsec from [gridding] section of config.
    If pixsize is not provided or is 0, calculates it as beamsize / 3.
    
    Parameters
    ----------
    config_path : str, optional
        Path to config.toml file. If None, searches for config.toml in standard locations.
    beamsize_deg : float, optional
        Beam size in degrees. Overrides config value if provided.
    pixsize_deg : float, optional
        Pixel size in degrees. Overrides config value if provided.
    
    Returns
    -------
    beamsize_deg : float
        Beam size in degrees
    pixsize_deg : float
        Pixel size in degrees (calculated as beamsize/3 if not specified)
    """
    gridding_cfg = {}
    if beamsize_deg is None or pixsize_deg is None:
        # Only load config when we actually need a value from it
        try:
            from oi_zeigt.basic_io import get_config
            cfg = get_config(config_path)
            gridding_cfg = cfg.get('gridding', {})
        except Exception:
            pass  # No config available — fall back to hard defaults below

    # Get beamsize from config (in arcseconds) or use provided value
    if beamsize_deg is None:
        beamsize_arcsec = gridding_cfg.get('beamsize_arcsec', 15.0)  # Default 15 arcsec
        beamsize_deg = beamsize_arcsec / 3600.0  # Convert arcseconds to degrees

    # Get pixel size from config (in arcseconds) or calculate it
    if pixsize_deg is None:
        pixsize_arcsec = gridding_cfg.get('pixel_size_arcsec', None)
        if pixsize_arcsec is None or pixsize_arcsec == 0:
            pixsize_deg = beamsize_deg / 3.0
        else:
            pixsize_deg = pixsize_arcsec / 3600.0

    return beamsize_deg, pixsize_deg


def _get_wcs_params_from_fits(hdul: fits.HDUList, beamsize_deg: float, 
                              pixsize: Optional[float] = None,
                              kernelsize_sigma: Optional[float] = None) -> Dict[str, Any]:
    """
    Extract and calculate WCS parameters for map creation.
    
    Parameters
    ----------
    hdul : astropy.io.fits.HDUList
        FITS file containing header information.
    beamsize_deg : float
        Beam size in degrees (for Gaussian kernel).
    pixsize : float, optional
        Pixel size in degrees. If None, uses beamsize/3.
    kernelsize_sigma : float, optional
        Kernel sigma in degrees. If None, uses beamsize/2.355.
    
    Returns
    -------
    dict
        Dictionary with WCS parameters and kernel settings.
    """
    if pixsize is None:
        pixsize = beamsize_deg / 3.0
    if kernelsize_sigma is None:
        kernelsize_sigma = beamsize_deg / 2.355

    header = hdul[0].header if len(hdul) > 0 else hdul[1].header
    
    # Get reference coordinate and velocity information
    crval2 = header.get('CRVAL2', 0.0)  # RA reference value
    crval3 = header.get('CRVAL3', 0.0)  # Dec reference value
    crpix1 = header.get('CRPIX1', 1.0)  # Spectral reference pixel
    velo_lsr = header.get('VELO-LSR', header.get('CRVAL3', 0.0))
    restfreq = header.get('RESTFREQ', header.get('RESTFRQ', 1.0))
    
    return {
        'pixsize': pixsize,
        'kernelsize_sigma': kernelsize_sigma,
        'crval_ra': crval2,
        'crval_dec': crval3,
        'crpix_spec': crpix1,
        'velo_lsr': velo_lsr,
        'restfreq': restfreq,
        'beamsize_deg': beamsize_deg,
    }


def create_wcs_header(naxis1: int, naxis2: int, naxis3: int,
                      crval_ra: float, crval_dec: float,
                      crpix1: float, crpix2: float, crpix3: float,
                      cdelt1: float, cdelt2: float, cdelt3: float,
                      beamsize_deg: float, velo_lsr: float,
                      restfreq: float, object_name: str = "",
                      telescop: str = "",
                      extra_header: Optional[Dict] = None) -> fits.Header:
    """
    Create a WCS-compliant FITS header for a 3D map.
    
    Parameters
    ----------
    naxis1, naxis2, naxis3 : int
        Map dimensions (spatial x, spatial y, spectral).
    crval_ra, crval_dec : float
        Reference RA/Dec in degrees.
    crpix1, crpix2, crpix3 : float
        Reference pixels (1-indexed).
    cdelt1, cdelt2, cdelt3 : float
        Pixel scales (RA, Dec, velocity).
    beamsize_deg : float
        Beam size in degrees.
    velo_lsr : float
        LSR velocity in m/s.
    restfreq : float
        Rest frequency in Hz.
    object_name : str, optional
        Object name for OBJECT keyword.
    extra_header : dict, optional
        Additional header keywords.
    
    Returns
    -------
    fits.Header
        WCS-compliant FITS header.
    """
    header = fits.Header()
    
    # Basic structure
    header['NAXIS'] = 3
    header['NAXIS1'] = naxis1
    header['NAXIS2'] = naxis2
    header['NAXIS3'] = naxis3
    
    # Beam info
    header['BMAJ'] = beamsize_deg
    header['BMIN'] = beamsize_deg
    header['BPA'] = 0.0
    header['BTYPE'] = 'Intensity'
    
    # Data info
    header['BUNIT'] = 'K'
    header['OBJECT'] = object_name.ljust(8)[:8]
    
    # WCS axes
    header['WCSAXES'] = 3
    header['CTYPE1'] = 'RA---SIN'
    header['CTYPE2'] = 'DEC--SIN'
    header['CTYPE3'] = 'VELO-LSR'
    header['CUNIT1'] = 'deg'
    header['CUNIT2'] = 'deg'
    header['CUNIT3'] = 'm/s'
    
    # Reference values
    header['CRVAL1'] = crval_ra
    header['CRVAL2'] = crval_dec
    header['CRVAL3'] = velo_lsr
    header['CRPIX1'] = crpix1
    header['CRPIX2'] = crpix2
    header['CRPIX3'] = crpix3
    
    # Pixel scales
    header['CDELT1'] = cdelt1
    header['CDELT2'] = cdelt2
    header['CDELT3'] = cdelt3
    
    # WCS projection
    header['PV2_1'] = 0.0
    header['PV2_2'] = 0.0
    header['LONPOLE'] = 180.0
    header['LATPOLE'] = 90.0
    
    # Reference system
    header['RADESYS'] = 'FK5'
    header['EQUINOX'] = 2000.0
    header['SPECSYS'] = 'LSRK'
    
    # Frequency reference
    header['RESTFRQ'] = restfreq
    header['TIMESYS'] = 'UTC'
    
    # Instrument info
    header['TELESCOP'] = telescop
    header['ORIGIN'] = 'OI-ZEIGT'
    
    # Add extra keywords if provided
    if extra_header:
        for key, value in extra_header.items():
            header[key] = value
    
    return header


def format_ra_hms(ra_deg: float) -> str:
    """
    Convert RA from degrees to hours:minutes:seconds format.
    
    Parameters
    ----------
    ra_deg : float
        Right ascension in degrees
    
    Returns
    -------
    str
        Formatted string like "08h30m45.2s"
    """
    from astropy.coordinates import Angle
    ra_angle = Angle(ra_deg, unit='deg')
    hms = ra_angle.hms
    return f"{int(hms.h):02d}h{int(hms.m):02d}m{hms.s:05.2f}s"


def ra_degrees_to_hms_formatter(axes_obj, ra_degrees):
    """
    Matplotlib formatter callback for RA axis to display hours:minutes:seconds.
    
    Parameters
    ----------
    axes_obj : matplotlib.axes.Axes
        The axes object
    ra_degrees : float
        RA value in degrees
    
    Returns
    -------
    str
        Formatted RA string in hms format
    """
    return format_ra_hms(ra_degrees)


def format_dec_dms(dec_deg: float) -> str:
    """
    Convert Dec from degrees to degrees:arcminutes:arcseconds format.
    
    Parameters
    ----------
    dec_deg : float
        Declination in degrees
    
    Returns
    -------
    str
        Formatted string like "+47°30'45.2\""
    """
    from astropy.coordinates import Angle
    dec_angle = Angle(dec_deg, unit='deg')
    dms = dec_angle.dms
    sign = '+' if dec_deg >= 0 else '-'
    return f"{sign}{int(abs(dms.d)):02d}°{int(abs(dms.m)):02d}'{abs(dms.s):05.2f}\""


def dec_degrees_to_dms_formatter(axes_obj, dec_degrees):
    """
    Matplotlib formatter callback for Dec axis to display degrees:arcminutes:arcseconds.
    
    Parameters
    ----------
    axes_obj : matplotlib.axes.Axes
        The axes object
    dec_degrees : float
        Dec value in degrees
    
    Returns
    -------
    str
        Formatted Dec string in dms format
    """
    return format_dec_dms(dec_degrees)


def _check_velocity_grid_uniform(
    table,
    velo_ref: float,
    deltav: float,
    velo_tol_frac: float = 0.1,
    deltav_rtol: float = 1e-3,
) -> None:
    """
    Warn if the VELOCITY/DELTAV columns are not uniform across the table.

    Gridding builds a single velocity axis from the first row
    (``velo_ref``/``deltav``) and applies it to every spectrum. That is only
    valid when all spectra were resampled onto a common LSR grid, which
    kalibrate normally does. If spectra carry different velocity references -
    e.g. flights observed with slightly different tuning frequencies that were
    never reconciled to a common systemic velocity - the single-axis assumption
    misregisters the data and blurs the co-added cube with no error raised.
    This emits a ``UserWarning`` describing the spread so the mismatch can be
    caught and fixed upstream (reconcile to a common Vsys before gridding).

    Parameters
    ----------
    table : FITS_rec
        Binary table with (optionally) VELOCITY and DELTAV columns, in m/s.
    velo_ref : float
        Reference velocity taken from row 0 (m/s) - the value gridding uses.
    deltav : float
        Channel spacing taken from row 0 (m/s) - the value gridding uses.
    velo_tol_frac : float, optional
        Allowed VELOCITY spread as a fraction of ``|deltav|`` (default 0.1,
        i.e. one tenth of a channel).
    deltav_rtol : float, optional
        Allowed relative spread in DELTAV (default 1e-3).
    """
    names = getattr(table, 'names', []) or []

    if 'DELTAV' in names and deltav != 0.0:
        dv = np.asarray(table['DELTAV'], dtype=np.float64)
        dv = dv[np.isfinite(dv)]
        if dv.size and float(np.ptp(dv)) > deltav_rtol * abs(deltav):
            warnings.warn(
                f"DELTAV is not uniform across the table (spread "
                f"{float(np.ptp(dv)):.4g} m/s, min {dv.min():.6g}, max {dv.max():.6g}); "
                f"gridding applies the row-0 value {deltav:.6g} m/s to all spectra. "
                f"Spectra do not share a common channel spacing - the cube will be "
                f"misgridded.",
                UserWarning,
            )

    if 'VELOCITY' in names:
        vel = np.asarray(table['VELOCITY'], dtype=np.float64)
        vel = vel[np.isfinite(vel)]
        tol = velo_tol_frac * abs(deltav)
        if vel.size and float(np.ptp(vel)) > tol:
            spread = float(np.ptp(vel))
            chans = f" = {spread / abs(deltav):.2f} channels" if deltav != 0.0 else ""
            warnings.warn(
                f"VELOCITY reference is not uniform across the table (spread "
                f"{spread:.4g} m/s{chans}; min {vel.min():.6g}, max {vel.max():.6g}); "
                f"gridding applies the row-0 value {velo_ref:.6g} m/s to all spectra. "
                f"Spectra are not registered to a common velocity grid (e.g. distinct "
                f"tuning frequencies) - the co-added cube will be blurred. Reconcile to "
                f"a common Vsys before gridding.",
                UserWarning,
            )


def _get_spectral_axis_params(hdul: fits.HDUList) -> Tuple[float, float, float, str, float, float]:
    """
    Extract spectral axis parameters from FITS file for velocity axis reconstruction.
    
    Extracts velocity reference, channel spacing, and related parameters from:
    - VELOCITY column (reference velocity in m/s, LSR convention)
    - DELTAV column (velocity spacing per channel in m/s)
    - VELDEF column (velocity definition, typically 'RADI-LSR')
    - RESTFREQ column (rest frequency in Hz)
    - SPECTRUM column shape (number of channels)
    - CRPIX1 from header (reference pixel for spectral axis, if available)
    
    Parameters
    ----------
    hdul : astropy.io.fits.HDUList
        FITS file containing spectral observations with SPECTRUM column.
    
    Returns
    -------
    velo_ref : float
        Reference velocity in m/s (from VELOCITY column)
    deltav : float
        Velocity spacing per channel in m/s (from DELTAV column)
    restfreq : float
        Rest frequency in Hz (from RESTFREQ column)
    veldef : str
        Velocity definition (from VELDEF column, typically 'RADI-LSR')
    nchans : int
        Number of spectral channels (from SPECTRUM column shape)
    crpix1_spec : float
        Reference pixel for spectral axis (from CRPIX1 in header, or default 1.0)
        FITS convention: 1-indexed, so channel_index = crpix1_spec - 1
    
    Notes
    -----
    All parameters are extracted as single values from the first observation,
    assuming they are consistent across all observations in the FITS file.
    
    The velocity of channel i is calculated as:
        v(i) = velo_ref + (i - (crpix1_spec - 1)) * deltav
    
    Example:
        If crpix1_spec=506, velo_ref=470000, deltav=500:
        - Channel 505: v = 470000 + (505 - 505) * 500 = 470000 m/s (reference)
        - Channel 0: v = 470000 + (0 - 505) * 500 = 217500 m/s
        - Channel 1263: v = 470000 + (1263 - 505) * 500 = 849000 m/s
    """
    # Get the binary table and its header (CLASS FITS stores spectral params here)
    table = hdul[1].data
    table_header = hdul[1].header
    primary_header = hdul[0].header if len(hdul) > 0 else fits.Header()

    def _header_get(key, fallback_keys=()):
        """Try table header, then primary header, then fallback keys in both."""
        for k in (key,) + tuple(fallback_keys):
            if k in table_header:
                return table_header[k]
        for k in (key,) + tuple(fallback_keys):
            if k in primary_header:
                return primary_header[k]
        return None

    # Extract reference velocity
    # Column takes precedence; fall back to VELO-LSR or VELOCITY header keyword
    if 'VELOCITY' in table.names:
        velo_ref = float(table['VELOCITY'][0])
    else:
        val = _header_get('VELO-LSR', ('VELOCITY',))
        if val is not None:
            velo_ref = float(val)
        else:
            warnings.warn("VELOCITY not found in columns or headers - using 0.0 m/s", UserWarning)
            velo_ref = 0.0

    # Extract velocity spacing per channel
    if 'DELTAV' in table.names:
        deltav = float(table['DELTAV'][0])
    else:
        val = _header_get('DELTAV')
        if val is not None:
            deltav = float(val)
        else:
            warnings.warn("DELTAV not found in columns or headers - using 1.0 m/s", UserWarning)
            deltav = 1.0

    # Guard: gridding applies this single velo_ref/deltav axis to every spectrum.
    # Verify the columns are actually uniform so an un-registered dataset (e.g.
    # flights with different tuning frequencies) is flagged rather than silently
    # blurring the cube.
    _check_velocity_grid_uniform(table, velo_ref, deltav)

    # Extract rest frequency
    if 'RESTFREQ' in table.names:
        restfreq = float(table['RESTFREQ'][0])
    else:
        val = _header_get('RESTFREQ', ('RESTFRQ',))
        if val is not None:
            restfreq = float(val)
        else:
            warnings.warn("RESTFREQ not found in columns or headers - using 1.0 Hz", UserWarning)
            restfreq = 1.0

    # Extract velocity definition
    if 'VELDEF' in table.names:
        veldef_val = table['VELDEF'][0]
        veldef = veldef_val.decode().strip() if isinstance(veldef_val, bytes) else str(veldef_val).strip()
    else:
        val = _header_get('VELDEF')
        if val is not None:
            veldef = str(val).strip()
        else:
            warnings.warn("VELDEF not found in columns or headers - using 'RADI-LSR'", UserWarning)
            veldef = 'RADI-LSR'

    # Get number of spectral channels from SPECTRUM
    if 'SPECTRUM' in table.names:
        nchans = table['SPECTRUM'].shape[1]
    else:
        raise ValueError("SPECTRUM column not found in FITS table")

    # Extract reference pixel for spectral axis: per-row CRPIX1 column
    # (median over rows) if present, else the headers.
    if 'CRPIX1' in table.names:
        crpix1_spec = table_crpix1(table)
    elif 'CRPIX1' in table_header:
        crpix1_spec = float(table_header['CRPIX1'])
    elif 'CRPIX1' in primary_header:
        crpix1_spec = float(primary_header['CRPIX1'])
    else:
        warnings.warn("CRPIX1 not found in headers - using 1.0 (first channel)", UserWarning)
        crpix1_spec = 1.0

    return velo_ref, deltav, restfreq, veldef, nchans, crpix1_spec


def _get_celestial_coords(hdul: fits.HDUList) -> Tuple[np.ndarray, np.ndarray]:
    """
    Extract celestial RA and Dec coordinates from FITS data.
    
    Uses CDELT2/CRVAL2 for RA and CDELT3/CRVAL3 for Dec, as these represent
    the actual celestial coordinates of the observations. The LONGITUDE/LATITUDE
    columns represent telescope position, not celestial position.
    
    CDELT2/CDELT3 are kalibrate's map offsets: offsets in the CLASS "radio"
    projection (kalibrate writes position.proj = PROJ_RADIO), i.e. projected
    sky offsets x = (RA - RA0) cos(Dec), y = Dec - Dec0. The exact inverse, as
    in GILDAS rel_to_abs (kernel/lib/gwcs/projec.f90, case p_radio), is:
        Dec = CRVAL3 + CDELT3
        RA = CRVAL2 + CDELT2 / cos(Dec)
    Adding CDELT2 to RA without the 1/cos(Dec) compresses every RA offset by
    cos(Dec) (68% at M51, Dec 47.2 deg): the LFA hexagon comes out distorted.
    See docs/ra_offset_projection_bug.md.
    
    Parameters
    ----------
    hdul : astropy.io.fits.HDUList
        FITS file with binary table containing CDELT2, CDELT3, etc.
    
    Returns
    -------
    ras : np.ndarray
        Right ascension values in degrees
    decs : np.ndarray
        Declination values in degrees
    """
    # Find ALL binary tables with CDELT2/CDELT3 columns (supports multi-HDU files)
    matrix_hdus = []
    for hdu in hdul:
        if hasattr(hdu, 'data') and hdu.data is not None and hasattr(hdu, 'header'):
            if 'CDELT2' in hdu.data.dtype.names and 'CDELT3' in hdu.data.dtype.names:
                matrix_hdus.append(hdu)
    
    if not matrix_hdus:
        raise ValueError("CDELT2 and CDELT3 columns not found in FITS table")
    
    # Concatenate data from all HDUs if multiple exist
    all_data = [hdu.data for hdu in matrix_hdus]
    data = np.concatenate(all_data) if len(all_data) > 1 else all_data[0]
    
    # CRVAL2/CRVAL3 are in the MATRIX table header (should be same across all HDUs)
    header = matrix_hdus[0].header
    
    # Extract reference values from header
    crval2 = header.get('CRVAL2', 0.0)  # RA reference in degrees
    crval3 = header.get('CRVAL3', 0.0)  # Dec reference in degrees
    
    # Extract delta values from data columns (radio projection, see docstring)
    # Dec = CRVAL3 + CDELT3 (in degrees)
    # RA = CRVAL2 + CDELT2 / cos(Dec) (in degrees)
    # float64: the columns are float32, which resolves RA ~200 deg only to ~0.05".
    decs = crval3 + np.asarray(data['CDELT3'], dtype=np.float64)
    ras = crval2 + np.asarray(data['CDELT2'], dtype=np.float64) / np.cos(np.radians(decs))

    return ras, decs


def grid_to_map(ras: np.ndarray, decs: np.ndarray, values: np.ndarray,
                wcs_header_dict: dict,
                naxis1: int, naxis2: int,
                beamsize_deg: float,
                ra_min: float, ra_max: float,
                dec_min: float, dec_max: float,
                weights: Optional[np.ndarray] = None,
                kernel_fwhm_deg: Optional[float] = None) -> np.ndarray:
    """
    Grid irregularly-sampled (ra, dec, value) triplets onto a regular 2D map.

    This is the single shared gridding primitive used by all higher-level
    mapping functions.  It handles the cygrid / scipy fallback and optional
    per-sample weight normalisation in one place.

    Parameters
    ----------
    ras, decs : ndarray
        Observation coordinates in degrees (already filtered to finite values).
    values : ndarray
        Data values to grid (raw, unscaled). When *weights* is provided the
        function internally computes ``values * weights`` as the numerator and
        ``weights`` as the denominator, yielding a proper weighted average.
    wcs_header_dict : dict
        WCS header produced by :func:`create_wcs_header` (used by cygrid).
    naxis1, naxis2 : int
        Output map dimensions (pixels).
    beamsize_deg : float
        Gaussian kernel FWHM in degrees.
    ra_min, ra_max, dec_min, dec_max : float
        Map bounding box in degrees (used by the scipy fallback).
    weights : ndarray, optional
        Per-sample weights.  When provided the output is
        ``sum(value * weight) / sum(weight)`` per pixel.

    Returns
    -------
    grid_map : ndarray, shape (naxis2, naxis1)
        Gridded 2D map as float32.
    """
    # Default kernel FWHM = full beam. GILDAS's xy_map defaults to beam/3
    # (xymap.f90: "fwhm = map%beam/3.0") to minimize resolution loss for
    # densely/Nyquist-sampled, high-S/N data — but that trades sensitivity
    # for resolution: fewer raw samples contribute to each output pixel, so
    # noise per pixel goes up. Tried beam/3 here against real M51 [OI] data
    # (faint, noise-limited, and not as densely cross-scan-sampled as GILDAS
    # assumes): peak SNR dropped from 21.4 to 13.5. For this kind of
    # noise-limited line, keeping the kernel >= beam (more smoothing, better
    # sensitivity) beats GILDAS's resolution-preserving default. Pass
    # --kernel-fwhm explicitly for the opposite tradeoff if a given dataset
    # is instead resolution-limited (high S/N, well-sampled).
    effective_fwhm = kernel_fwhm_deg if kernel_fwhm_deg is not None else beamsize_deg
    kernelsize_sigma = effective_fwhm / 2.355
    kernel_type = 'gauss1d'
    kernel_params = (kernelsize_sigma,)
    # Support (truncation radius) of 3*FWHM, matching GILDAS's xy_map
    # convention exactly (Gaussian value there is ~1.45e-11, fully
    # negligible) — wider than a naive 3*sigma (~3*FWHM/2.355), which would
    # still retain ~1% of the peak at the truncation edge.
    kernel_support = 3.0 * effective_fwhm
    hpx_maxres = kernelsize_sigma / 2.0

    if HAS_CYGRID:
        # Numerator: sum(kernel * value * weight) — pre-multiply so the kernel
        # applies to the weighted values, not just the raw values.
        gridded_values = values * weights if weights is not None else values
        gridder = cygrid.WcsGrid(wcs_header_dict)
        gridder.set_kernel(kernel_type, kernel_params, kernel_support, hpx_maxres)
        gridder.grid(ras, decs, gridded_values)
        with np.errstate(invalid='ignore'):
            grid_map = gridder.get_datacube().squeeze().astype(np.float64)

        # Always grid a coverage map to identify empty pixels (cygrid fills them with 0)
        gridder_cov = cygrid.WcsGrid(wcs_header_dict)
        gridder_cov.set_kernel(kernel_type, kernel_params, kernel_support, hpx_maxres)
        gridder_cov.grid(ras, decs, np.ones(len(ras), dtype=np.float64))
        with np.errstate(invalid='ignore'):
            coverage = gridder_cov.get_datacube().squeeze()

        if weights is not None:
            # Denominator: sum(kernel * weight)
            gridder_w = cygrid.WcsGrid(wcs_header_dict)
            gridder_w.set_kernel(kernel_type, kernel_params, kernel_support, hpx_maxres)
            gridder_w.grid(ras, decs, weights)
            with np.errstate(invalid='ignore'):
                grid_weights = gridder_w.get_datacube().squeeze().astype(np.float64)
            out = np.full_like(grid_map, np.nan)
            grid_map = np.divide(grid_map, grid_weights,
                                 where=grid_weights > 0,
                                 out=out)

        # Set pixels with no coverage to NaN
        grid_map[coverage <= 0] = np.nan
        grid_map = grid_map.astype(np.float32)
    else:
        warnings.warn(
            "cygrid not available - using scipy griddata instead "
            "(less accurate for astronomical data)", UserWarning)
        from scipy.interpolate import griddata
        ra_grid = np.linspace(ra_min, ra_max, naxis1)
        dec_grid = np.linspace(dec_min, dec_max, naxis2)
        ra_mesh, dec_mesh = np.meshgrid(ra_grid, dec_grid)
        points = np.column_stack([ras, decs])

        gridded_values = values * weights if weights is not None else values
        grid_map = griddata(points, gridded_values, (ra_mesh, dec_mesh),
                            method='linear').astype(np.float32)
        if weights is not None:
            grid_weights = griddata(points, weights, (ra_mesh, dec_mesh),
                                    method='linear').astype(np.float32)
            out = np.full_like(grid_map, np.nan)
            grid_map = np.divide(grid_map, grid_weights,
                                 where=grid_weights > 0,
                                 out=out)

    return grid_map


def grid_spectra_with_cygrid(hdul: fits.HDUList,
                             data_column: str,
                             ra_col: str, dec_col: str,
                             channel_idx: Optional[int] = None,
                             object_filter: Optional[str] = None,
                             beamsize_deg: float = 0.25,
                             pixsize: Optional[float] = None,
                             figsize: Tuple[int, int] = (10, 8),
                             show_scatter: bool = False) -> Tuple[np.ndarray, fits.Header, plt.Figure]:
    """
    Grid spectral data using cygrid with proper WCS.
    
    Parameters
    ----------
    hdul : astropy.io.fits.HDUList
        FITS HDUList containing the data.
    data_column : str
        Column name containing the data to grid.
    ra_col, dec_col : str
        Column names for RA and Dec coordinates.
    channel_idx : int, optional
        Channel index for 2D spectral data.
    object_filter : str, optional
        Filter by object name (substring match).
    beamsize_deg : float, optional
        Beam size in degrees (default: 0.25).
    pixsize : float, optional
        Map pixel size in degrees. If None, uses beamsize/3.
    figsize : tuple, optional
        Figure size (default: (10, 8)).
    show_scatter : bool, optional
        Whether to show observation points as scatter plot (default: False).
    
    Returns
    -------
    grid_map : np.ndarray
        Gridded 2D map array.
    wcs_header : fits.Header
        WCS-compliant FITS header.
    fig : matplotlib.figure.Figure
        Visualization figure.
    """
    # Find ALL binary tables with the required data column (supports multi-HDU files)
    matrix_hdus = []
    for hdu in hdul:
        if hasattr(hdu, 'data') and hdu.data is not None:
            if data_column in hdu.data.dtype.names:
                matrix_hdus.append(hdu)
    
    if not matrix_hdus:
        raise ValueError(f"Column '{data_column}' not found")
    
    # Concatenate data from all HDUs if multiple exist
    all_data = [hdu.data for hdu in matrix_hdus]
    data = np.concatenate(all_data) if len(all_data) > 1 else all_data[0]
    
    # Extract celestial coordinates (RA/Dec from CDELT2/3 and CRVAL2/3)
    ras, decs = _get_celestial_coords(hdul)
    
    # Extract data column values
    values = data[data_column]
    
    # Apply object filter if specified
    if object_filter:
        if 'OBJECT' not in data.dtype.names:
            raise ValueError("OBJECT column not found")
        
        objects = data['OBJECT']
        mask = np.array([
            object_filter.upper() in (obj.decode().strip() if isinstance(obj, bytes) else str(obj).strip()).upper()
            for obj in objects
        ])
        
        if not np.any(mask):
            raise ValueError(f"No objects matching '{object_filter}' found")
        
        ras = ras[mask]
        decs = decs[mask]
        values = values[mask]
    
    # Handle 2D data
    if values.ndim == 2:
        if channel_idx is None:
            raise ValueError("channel_idx required for 2D data")
        values = values[:, channel_idx]
    
    # Filter valid points
    valid = np.isfinite(ras) & np.isfinite(decs) & np.isfinite(values)
    ras = ras[valid]
    decs = decs[valid]
    values = values[valid]
    
    if len(values) < 3:
        raise ValueError(f"Not enough valid data points ({len(values)})")
    
    # Calculate map bounds and WCS parameters
    if pixsize is None:
        pixsize = beamsize_deg / 3.0
    
    ra_min, ra_max = np.min(ras), np.max(ras)
    dec_min, dec_max = np.min(decs), np.max(decs)
    
    naxis1 = int((ra_max - ra_min) / abs(pixsize)) + 1
    naxis2 = int((dec_max - dec_min) / abs(pixsize)) + 1
    
    # Map center
    ra_center = (ra_max + ra_min) / 2.0
    dec_center = (dec_max + dec_min) / 2.0
    
    crpix1 = (ra_max - ra_center) / pixsize + 1
    crpix2 = (dec_center - dec_min) / pixsize + 1
    
    # Get additional WCS info from header
    header = hdul[0].header if len(hdul) > 0 else {}
    crpix_spec = header.get('CRPIX1', 1.0)
    velo_lsr = header.get('VELO-LSR', 0.0)
    restfreq = header.get('RESTFREQ', 1.0)
    object_name = header.get('OBJECT', '').strip()
    
    wcs_header_dict = create_wcs_header(
        naxis1=naxis1, naxis2=naxis2, naxis3=1,
        crval_ra=ra_center, crval_dec=dec_center,
        crpix1=crpix1, crpix2=crpix2, crpix3=crpix_spec,
        cdelt1=-pixsize, cdelt2=pixsize, cdelt3=1.0,
        beamsize_deg=beamsize_deg,
        velo_lsr=velo_lsr,
        restfreq=restfreq,
        object_name=object_name,
    )

    grid_map = grid_to_map(
        ras, decs, values,
        wcs_header_dict=wcs_header_dict,
        naxis1=naxis1, naxis2=naxis2,
        beamsize_deg=beamsize_deg,
        ra_min=ra_min, ra_max=ra_max,
        dec_min=dec_min, dec_max=dec_max,
    )

    # Create WCS header
    wcs_header = create_wcs_header(
        naxis1=naxis1, naxis2=naxis2, naxis3=1,
        crval_ra=ra_center, crval_dec=dec_center,
        crpix1=crpix1, crpix2=crpix2, crpix3=crpix_spec,
        cdelt1=-pixsize, cdelt2=pixsize, cdelt3=1.0,
        beamsize_deg=beamsize_deg,
        velo_lsr=velo_lsr,
        restfreq=restfreq,
        object_name=object_name,
    )
    
    # Create figure
    fig, ax = plt.subplots(figsize=figsize)
    
    # Flip the map horizontally so it aligns with inverted RA axis (RA increases right to left)
    grid_map_display = np.fliplr(grid_map)
    
    im = ax.imshow(grid_map_display, origin='lower', extent=[ra_min, ra_max, dec_min, dec_max],
                   cmap='viridis', aspect='auto')
    if show_scatter:
        ax.scatter(ras, decs, c='red', s=5, alpha=0.3, label='Observations')
    
    # Format RA axis in hours:minutes:seconds and Dec axis in degrees:arcminutes:arcseconds
    from matplotlib.ticker import FuncFormatter, MaxNLocator
    ra_formatter = FuncFormatter(lambda x, p: format_ra_hms(x))
    dec_formatter = FuncFormatter(lambda x, p: format_dec_dms(x))
    ax.xaxis.set_major_formatter(ra_formatter)
    ax.yaxis.set_major_formatter(dec_formatter)
    # Reduce number of x-axis ticks to avoid crowding
    ax.xaxis.set_major_locator(MaxNLocator(nbins=5, integer=False))
    
    # Invert x-axis so RA increases from right to left (standard sky projection)
    ax.invert_xaxis()
    
    ax.set_xlabel('RA (h:m:s)', fontsize=12)
    ax.set_ylabel('Dec (°:′:″)', fontsize=12)
    gridding_method = 'cygrid' if HAS_CYGRID else 'scipy'
    ax.set_title(f'Spatial Map ({gridding_method}, beam={beamsize_deg*3600:.1f}″)', 
                fontsize=13, fontweight='bold')
    if show_scatter:
        ax.legend(loc='best')
    
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label('Intensity (K)', fontsize=11)
    
    fig.tight_layout()
    
    return grid_map, wcs_header, fig


def create_map_from_column(hdul: fits.HDUList,
                          column_name: str,
                          ra_col: str = 'LONGITUDE',
                          dec_col: str = 'LATITUDE',
                          object_filter: Optional[str] = None,
                          beamsize_deg: float = 0.25,
                          pixsize: Optional[float] = None,
                          figsize: Tuple[int, int] = (10, 8),
                          show_scatter: bool = False) -> Tuple[np.ndarray, fits.Header, plt.Figure]:
    """
    Create a spatial map from a scalar column (quality metric, temperature, etc.).
    Uses proper WCS-based gridding with cygrid.
    """
    return grid_spectra_with_cygrid(
        hdul,
        data_column=column_name,
        ra_col=ra_col,
        dec_col=dec_col,
        object_filter=object_filter,
        beamsize_deg=beamsize_deg,
        pixsize=pixsize,
        figsize=figsize,
        show_scatter=show_scatter
    )


def create_integrated_map(hdul: fits.HDUList,
                         spectrum_column: str = 'SPECTRUM',
                         ra_col: str = 'LONGITUDE',
                         dec_col: str = 'LATITUDE',
                         object_filter: Optional[str] = None,
                         beamsize_deg: float = 0.25,
                         pixsize: Optional[float] = None,
                         figsize: Tuple[int, int] = (10, 8),
                         show_scatter: bool = False,
                         velocity_range: Optional[Tuple[float, float]] = None,
                         weight_column: Optional[str] = None,
                         channel_weights: bool = False) -> Tuple[np.ndarray, fits.Header, plt.Figure]:
    """
    Create a spatial map from integrated spectral intensity with optional per-spectrum weighting.

    Integrates spectra across all channels (or specified velocity range), then grids with proper WCS.

    Parameters
    ----------
    velocity_range : tuple or None
        (v_min, v_max) in km/s to restrict integration to this velocity range.
        If None, integrates entire spectrum.
    weight_column : str, optional
        Column name for per-spectrum weighting (e.g., 'RMSRATIO' or 'RMSRATIOB').
        If None, all spectra are weighted equally (uniform weighting).
        For RMSRATIO/RMSRATIOB, weights are calculated as:
            weight = exp(-(value - 1.0)^2 / (2*sigma^2))  with sigma=0.5
        Final map is normalized: map = sum(intensity*weight) / sum(weights).
    channel_weights : bool, optional
        If True, compute per-channel weights from the TSYS and TAU_SIG calibration
        spectra linked via TSYS_INDEX / TAU_SIG_INDEX columns (added by prepare_for_pca).
        Weight per channel: w_{i,c} = exp(-tau_{i,c}) / tsys_{i,c}
        The integration becomes: sum_c(w_{i,c} * T_{i,c}) / sum_c(w_{i,c}).
        Falls back silently to uniform channel weights when the index columns are
        absent, the index is -1, or the calibration spectrum contains non-finite values.
    """
    # Find ALL binary tables with the required spectrum column (supports multi-HDU files)
    matrix_hdus = []
    for hdu in hdul:
        if hasattr(hdu, 'data') and hdu.data is not None:
            if spectrum_column in hdu.data.dtype.names:
                matrix_hdus.append(hdu)
    
    if not matrix_hdus:
        raise ValueError(f"Column '{spectrum_column}' not found")
    
    # Concatenate data from all HDUs if multiple exist
    all_data = [hdu.data for hdu in matrix_hdus]
    data = np.concatenate(all_data) if len(all_data) > 1 else all_data[0]
    
    # Extract celestial coordinates (RA/Dec from CDELT2/3 and CRVAL2/3)
    ras, decs = _get_celestial_coords(hdul)
    spectra = data[spectrum_column]
    
    # Extract coordinates
    if object_filter:
        if 'OBJECT' not in data.dtype.names:
            raise ValueError("OBJECT column not found")
        
        objects = data['OBJECT']
        mask = np.array([
            object_filter.upper() in (obj.decode().strip() if isinstance(obj, bytes) else str(obj).strip()).upper()
            for obj in objects
        ])
        
        if not np.any(mask):
            raise ValueError(f"No objects matching '{object_filter}' found")
        
        ras = ras[mask]
        decs = decs[mask]
        spectra = spectra[mask]
    
    # Apply velocity range filtering if specified
    if velocity_range is not None:
        v_min, v_max = velocity_range
        
        # Get header for velocity axis computation
        header = hdul[0].header if len(hdul) > 0 else {}
        
        # Extract velocity calibration parameters (prefer VELOCITY/DELTAV columns, fallback to header)
        try:
            n_channels = spectra.shape[1]
            
            # Try to get velocity from data columns first (units: m/s)
            if 'VELOCITY' in data.dtype.names and 'DELTAV' in data.dtype.names:
                velo_ref = float(data['VELOCITY'][0]) / 1000.0  # Convert m/s to km/s
                deltav = float(data['DELTAV'][0]) / 1000.0      # Convert m/s to km/s
                # Same single-axis assumption as create_spectral_datacube: flag
                # non-uniform grids (checker works in m/s, matching the columns).
                _check_velocity_grid_uniform(data, velo_ref * 1000.0, deltav * 1000.0)
            else:
                # Fallback to header (CRVAL3/CDELT3 units depend on CUNIT3)
                velo_ref = header.get('CRVAL3', 0.0)
                deltav = header.get('CDELT3', 1.0)
                
                # Check units and convert if necessary
                cunit3 = header.get('CUNIT3', 'm/s')
                if 'm/s' in str(cunit3).lower():
                    velo_ref = velo_ref / 1000.0  # Convert m/s to km/s
                    deltav = deltav / 1000.0
            
            crpix1_spec = table_crpix1(data, header)
            
            # Compute velocity axis in km/s
            channel_indices = np.arange(n_channels, dtype=np.float64)
            velocity_axis = velo_ref + (channel_indices - (crpix1_spec - 1.0)) * deltav
            
            # Find channel range for this velocity window
            ch_min = np.argmin(np.abs(velocity_axis - v_min))
            ch_max = np.argmin(np.abs(velocity_axis - v_max))
            ch_min, ch_max = min(ch_min, ch_max), max(ch_min, ch_max)
            
            # Include both boundary channels
            ch_max = min(ch_max + 1, n_channels)
            
            # Slice spectra to velocity range
            spectra = spectra[:, ch_min:ch_max]
            
        except (KeyError, ValueError, TypeError, AttributeError) as e:
            raise ValueError(f"Cannot compute velocity axis from data/header: {e}")
    
    # Extract per-spectrum weights if requested
    weights = None
    if weight_column is not None:
        # Fallback: if RMSRATIOB requested but absent, try RMSRATIO
        if weight_column not in data.dtype.names and weight_column == 'RMSRATIOB' and 'RMSRATIO' in data.dtype.names:
            warnings.warn("RMSRATIOB column not found, falling back to RMSRATIO for per-spectrum weighting.", UserWarning)
            weight_column = 'RMSRATIO'
        if weight_column not in data.dtype.names:
            raise ValueError(f"Weight column '{weight_column}' not found in FITS file")
        
        weight_values = data[weight_column]
        
        # Apply object filter to weights if applicable
        if object_filter:
            weight_values = weight_values[mask]
        
        # Calculate weights based on quality metric
        # For RMSRATIO: 1.0 is the theoretical radiometer floor. Spectra at or
        # below it (RMSRATIO <= 1) keep full weight 1.0; only excess noise
        # (RMSRATIO > 1) is penalised via a one-sided Gaussian falloff.
        # sigma=0.5 gives 0.61 at 1.5, 0.14 at 2.0. Deliberately NOT symmetric
        # about 1.0 — quieter-than-theoretical spectra are never down-weighted.
        weights = np.zeros_like(weight_values, dtype=np.float64)
        valid_weights = np.isfinite(weight_values)

        if 'RMSRATIO' in weight_column.upper():
            # One-sided Gaussian weighting for RMSRATIO (penalise > 1 only)
            sigma = 0.5
            excess = np.clip(weight_values[valid_weights] - 1.0, 0.0, None)
            weights[valid_weights] = np.exp(-(excess ** 2) / (2 * sigma ** 2))
        else:
            # Generic inverse weighting for other columns (lower is better)
            weights[valid_weights & (weight_values != 0)] = 1.0 / np.abs(weight_values[valid_weights & (weight_values != 0)])
            # Normalize to [0, 1] range
            max_weight = np.max(weights[weights > 0])
            if max_weight > 0:
                weights[weights > 0] /= max_weight
    
    # Integrated intensity — always original unmodified values (nansum over channels).
    # Per-channel and per-spectrum weights affect only how spectra are blended
    # spatially by the gridder (via the weights argument to grid_to_map), following
    # the same principle as CLASS map.template.class:
    #   grid(T * w) / grid(w)  →  weighted average preserving original T units.
    integrated = np.nansum(spectra, axis=1)

    # --- Per-channel gridding weight from TSYS / TAU_SIG calibration spectra ---
    # For each spectrum i: W_i = nansum_c( exp(-tau_{i,c}) / T_sys_{i,c} )
    # This collapses the per-channel quality to a per-spectrum scalar weight.
    # Spectra with lower opacity and lower T_sys across more channels get higher weight.
    # The weight is used only for spatial gridding; integrated intensities are untouched.
    chan_quality = None
    if channel_weights:
        has_cols = ('TSYS_INDEX' in data.dtype.names and 'TAU_SIG_INDEX' in data.dtype.names)
        if not has_cols:
            warnings.warn(
                "channel_weights=True requested but TSYS_INDEX / TAU_SIG_INDEX columns "
                "not found in file (available after prepare_for_pca). "
                "Falling back to uniform channel weights.", UserWarning)
        else:
            all_spectra = data[spectrum_column]
            tsys_idx = data['TSYS_INDEX']
            tau_idx  = data['TAU_SIG_INDEX']
            if object_filter:
                tsys_idx = tsys_idx[mask]
                tau_idx  = tau_idx[mask]
            n_spec, n_chan = spectra.shape
            # NaN sentinel: spectra with no calibration match get NaN initially
            chan_quality = np.full(n_spec, np.nan, dtype=np.float64)
            for i in range(n_spec):
                ti = int(tsys_idx[i])
                ai = int(tau_idx[i])
                if ti < 0 or ai < 0 or ti >= len(all_spectra) or ai >= len(all_spectra):
                    continue  # stale / missing index → leave as NaN
                tsys_spec = np.asarray(all_spectra[ti], dtype=np.float64)
                tau_spec  = np.asarray(all_spectra[ai], dtype=np.float64)
                if velocity_range is not None:
                    tsys_spec = tsys_spec[ch_min:ch_max]
                    tau_spec  = tau_spec[ch_min:ch_max]
                w_ch = np.where(
                    np.isfinite(tsys_spec) & np.isfinite(tau_spec) & (tsys_spec > 0),
                    (np.exp(-tau_spec) / tsys_spec) ** 2,
                    np.nan,
                )
                w_sum = np.nansum(w_ch)
                if np.isfinite(w_sum) and w_sum > 0:
                    chan_quality[i] = w_sum

            # Normalize: divide by median of calibrated values so the typical
            # spectrum gets weight=1.0. Spectra with no calibration match also
            # get weight=1.0 (median = neutral), avoiding the magnitude mismatch
            # that would otherwise suppress or dominate entire sky regions.
            calibrated_mask = np.isfinite(chan_quality)
            if np.any(calibrated_mask):
                median_cq = np.median(chan_quality[calibrated_mask])
                if median_cq > 0:
                    chan_quality /= median_cq
            chan_quality = np.where(np.isfinite(chan_quality), chan_quality, 1.0)

    # Combine per-spectrum (RMSRATIO) and per-channel (tau/T_sys) gridding weights.
    # grid_to_map computes: Σ(kernel * gridding_weight * integrated) / Σ(kernel * gridding_weight)
    # which is a weighted average of the original integrated intensities.
    if weights is not None and chan_quality is not None:
        gridding_weights = weights * chan_quality
    elif chan_quality is not None:
        gridding_weights = chan_quality
    else:
        gridding_weights = weights  # may be None → uniform

    # Filter valid points
    valid = np.isfinite(ras) & np.isfinite(decs) & np.isfinite(integrated)
    ras = ras[valid]
    decs = decs[valid]
    integrated = integrated[valid]

    if gridding_weights is not None:
        gridding_weights = gridding_weights[valid]
    
    if len(integrated) < 3:
        raise ValueError(f"Not enough valid data points ({len(integrated)})")
    
    # Calculate map bounds and WCS parameters
    if pixsize is None:
        pixsize = beamsize_deg / 3.0
    
    ra_min, ra_max = np.min(ras), np.max(ras)
    dec_min, dec_max = np.min(decs), np.max(decs)
    
    naxis1 = int((ra_max - ra_min) / abs(pixsize)) + 1
    naxis2 = int((dec_max - dec_min) / abs(pixsize)) + 1
    
    # Map center
    ra_center = (ra_max + ra_min) / 2.0
    dec_center = (dec_max + dec_min) / 2.0
    
    crpix1 = (ra_max - ra_center) / pixsize + 1
    crpix2 = (dec_center - dec_min) / pixsize + 1
    
    # Get additional WCS info from header
    header = hdul[0].header if len(hdul) > 0 else {}
    crpix_spec = header.get('CRPIX1', 1.0)
    velo_lsr = header.get('VELO-LSR', 0.0)
    restfreq = header.get('RESTFREQ', 1.0)
    object_name = header.get('OBJECT', '').strip()
    
    wcs_header_dict = create_wcs_header(
        naxis1=naxis1, naxis2=naxis2, naxis3=1,
        crval_ra=ra_center, crval_dec=dec_center,
        crpix1=crpix1, crpix2=crpix2, crpix3=crpix_spec,
        cdelt1=-pixsize, cdelt2=pixsize, cdelt3=1.0,
        beamsize_deg=beamsize_deg,
        velo_lsr=velo_lsr,
        restfreq=restfreq,
        object_name=object_filter or "",
    )

    # grid_to_map expects values pre-multiplied by weights so that:
    #   output = Σ(K * W * integrated) / Σ(K * W)  →  weighted average in original units
    values_to_grid = integrated * gridding_weights if gridding_weights is not None else integrated
    grid_map = grid_to_map(
        ras, decs, values_to_grid,
        wcs_header_dict=wcs_header_dict,
        naxis1=naxis1, naxis2=naxis2,
        beamsize_deg=beamsize_deg,
        ra_min=ra_min, ra_max=ra_max,
        dec_min=dec_min, dec_max=dec_max,
        weights=gridding_weights,
    )
    
    # Create WCS header
    wcs_header = create_wcs_header(
        naxis1=naxis1, naxis2=naxis2, naxis3=1,
        crval_ra=ra_center, crval_dec=dec_center,
        crpix1=crpix1, crpix2=crpix2, crpix3=crpix_spec,
        cdelt1=-pixsize, cdelt2=pixsize, cdelt3=1.0,
        beamsize_deg=beamsize_deg,
        velo_lsr=velo_lsr,
        restfreq=restfreq,
        object_name=object_name,
    )
    
    # Create figure
    fig, ax = plt.subplots(figsize=figsize)
    
    # Flip the map horizontally so it aligns with inverted RA axis (RA increases right to left)
    grid_map_display = np.fliplr(grid_map)
    
    im = ax.imshow(grid_map_display, origin='lower', extent=[ra_min, ra_max, dec_min, dec_max],
                   cmap='viridis', aspect='auto')
    if show_scatter:
        ax.scatter(ras, decs, c='red', s=5, alpha=0.3, label='Observations')
    
    # Format RA axis in hours:minutes:seconds and Dec axis in degrees:arcminutes:arcseconds
    from matplotlib.ticker import FuncFormatter, MaxNLocator
    ra_formatter = FuncFormatter(lambda x, p: format_ra_hms(x))
    dec_formatter = FuncFormatter(lambda x, p: format_dec_dms(x))
    ax.xaxis.set_major_formatter(ra_formatter)
    ax.yaxis.set_major_formatter(dec_formatter)
    # Reduce number of x-axis ticks to avoid crowding
    ax.xaxis.set_major_locator(MaxNLocator(nbins=5, integer=False))
    
    # Invert x-axis so RA increases from right to left (standard sky projection)
    ax.invert_xaxis()
    
    ax.set_xlabel('RA (h:m:s)', fontsize=12)
    ax.set_ylabel('Dec (°:′:″)', fontsize=12)
    gridding_method = 'cygrid' if HAS_CYGRID else 'scipy'
    weight_info = f', weighted by {weight_column}' if weight_column else ''
    ax.set_title(f'Integrated Intensity Map ({gridding_method}, beam={beamsize_deg*3600:.1f}″{weight_info})',
                fontsize=13, fontweight='bold')
    if show_scatter:
        ax.legend(loc='best')
    
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label('Integrated Intensity (K·m/s)', fontsize=11)
    
    fig.tight_layout()
    
    return grid_map, wcs_header, fig


def save_map_to_fits(grid_map: np.ndarray, wcs_header: fits.Header,
                     output_file: str, beam_maj_deg: float, beam_min_deg: float = None,
                     beam_pa_deg: float = 0.0, overwrite: bool = True,
                     spectral_params: Optional[Dict[str, Any]] = None,
                     coverage_map: Optional[np.ndarray] = None,
                     weight_map: Optional[np.ndarray] = None,
                     weight_cube: Optional[np.ndarray] = None,
                     weight_column: Optional[str] = None) -> None:
    """
    Save a gridded map to a FITS file with proper WCS and beam information.
    
    Optionally includes a binary table with spectral axis information (CHANNEL and VELOCITY)
    for easy direct access without requiring WCS libraries.
    
    Parameters
    ----------
    grid_map : np.ndarray
        2D or 3D map data to save. Can be 2D (spatial only) or 3D (spatial + spectral).
    wcs_header : astropy.io.fits.Header
        FITS header containing WCS information (from create_wcs_header).
    output_file : str
        Output FITS filename.
    beam_maj_deg : float
        Beam major axis in degrees (BMAJ).
    beam_min_deg : float, optional
        Beam minor axis in degrees (BMIN). If None, uses beam_maj_deg.
    beam_pa_deg : float, optional
        Beam position angle in degrees (BPA). Default is 0.0.
    overwrite : bool, optional
        Whether to overwrite existing file. Default is True.
    spectral_params : dict, optional
        Dictionary with spectral axis parameters:
        - 'nvel': number of velocity channels
        - 'velo_ref': reference velocity (m/s)
        - 'crpix1_spec': reference pixel (FITS 1-indexed)
        - 'deltav': velocity step per channel (m/s)
        - 'restfreq': rest frequency (Hz)
        - 'veldef': velocity definition string
        
        If provided, creates a binary table HDU with CHANNEL and VELOCITY columns.
    
    Returns
    -------
    None
        Writes FITS file to disk.
    """
    if beam_min_deg is None:
        beam_min_deg = beam_maj_deg
    
    # Create header from WCS information
    header = fits.Header({k: v if not isinstance(v, tuple) else v[0] 
                         for k, v in wcs_header.items()})
    
    # Add beam information
    header.set('BMAJ', beam_maj_deg, comment='beam major axis in degrees')
    header.set('BMIN', beam_min_deg, comment='beam minor axis in degrees')
    header.set('BPA', beam_pa_deg, comment='beam position angle in degrees')
    
    # Create primary HDU with data and header
    primary_hdu = fits.PrimaryHDU(data=grid_map.astype(np.float32), header=header)
    
    # Create HDU list
    hdul = fits.HDUList([primary_hdu])
    
    # Add spectral axis table if parameters provided
    if spectral_params is not None:
        nvel = spectral_params['nvel']
        velo_ref = spectral_params['velo_ref']
        crpix1_spec = spectral_params['crpix1_spec']
        deltav = spectral_params['deltav']
        restfreq = spectral_params['restfreq']
        veldef = spectral_params['veldef']
        
        # Calculate velocities for all channels
        channels = np.arange(nvel)
        velocities = velo_ref + (channels - (crpix1_spec - 1)) * deltav
        
        # Create table
        spec_table = Table()
        spec_table['CHANNEL'] = channels
        spec_table['VELOCITY'] = velocities * u.m / u.s
        
        # Create binary table HDU
        spec_bintable = fits.BinTableHDU(spec_table, name='SPECTRUM')
        spec_bintable.header['EXTVER'] = 1
        spec_bintable.header['EXTLEVEL'] = 1
        spec_bintable.header['CRPIX1'] = (crpix1_spec, 'Reference pixel (FITS 1-indexed)')
        spec_bintable.header['RESTFRQ'] = (restfreq, 'Rest frequency (Hz)')
        spec_bintable.header['VELDEF'] = (veldef.strip(), 'Velocity definition')
        spec_bintable.header['COMMENT'] = 'Spectral axis information: direct access to channel-velocity mapping'
        
        # Add to HDU list
        hdul.append(spec_bintable)

    # Add coverage map as ImageHDU if provided
    if coverage_map is not None:
        cov_hdu = fits.ImageHDU(data=coverage_map.astype(np.float32), name='COVERAGE')
        cov_hdu.header['COMMENT'] = 'Gridding coverage map: sum of kernel weights per pixel'
        hdul.append(cov_hdu)

    # Add weight map as ImageHDU if provided
    if weight_map is not None:
        wm_hdu = fits.ImageHDU(data=weight_map.astype(np.float32), name='WEIGHT_MAP')
        col_label = weight_column if weight_column else 'user-supplied'
        wm_hdu.header['WGTCOL'] = (col_label, 'Source column for per-spectrum weights')
        wm_hdu.header['COMMENT'] = ('Kernel-weighted sum of per-spectrum weights per pixel. '
                                    'Proportional to the denominator of the weighted average.')
        hdul.append(wm_hdu)

    # Add full 3D weight cube as ImageHDU if provided (only built when per-channel
    # weights are active, since that's the only case where the weight sum
    # actually varies channel to channel).
    if weight_cube is not None:
        wc_hdu = fits.ImageHDU(data=weight_cube.astype(np.float32), name='WEIGHT_CUBE')
        wc_hdu.header['COMMENT'] = ('Per-channel kernel-weighted sum of combined spectrum x '
                                    'channel weights. Denominator of the per-channel weighted '
                                    'average; varies channel to channel (unlike WEIGHT_MAP).')
        hdul.append(wc_hdu)

    # Write to file
    hdul.writeto(output_file, overwrite=overwrite)
    
    print(f"✓ FITS datacube saved: {output_file}")
    print(f"  Data shape: {grid_map.shape}")
    print(f"  Beam: BMAJ={beam_maj_deg:.4f}°, BMIN={beam_min_deg:.4f}°, BPA={beam_pa_deg:.1f}°")
    
    # Print info about spectral table if included
    if spectral_params is not None:
        nvel = spectral_params['nvel']
        velo_ref = spectral_params['velo_ref'] / 1000.0  # Convert to km/s for display
        print(f"  Spectral axis: {nvel} channels (SPECTRUM HDU)")
        print(f"    Reference velocity: {velo_ref:.2f} km/s at channel {spectral_params['crpix1_spec']-1:.0f}")

    if weight_cube is not None:
        print(f"  Weight cube: {weight_cube.shape} (WEIGHT_CUBE HDU)")


# ---------------------------------------------------------------------------
# Multiprocessing worker functions for parallel channel gridding.
# Must be module-level (not nested) so they can be pickled by multiprocessing.
# ---------------------------------------------------------------------------

_gridding_worker_state: dict = {}


def _init_gridding_worker(ras, decs, wcs_header_dict,
                          naxis1, naxis2, beamsize_deg,
                          ra_min, ra_max, dec_min, dec_max,
                          kernel_fwhm_deg=None, need_weight_grid=False):
    """Pool initializer — stores shared geometry in each worker process once."""
    global _gridding_worker_state
    _gridding_worker_state = dict(
        ras=ras, decs=decs,
        wcs_header_dict=wcs_header_dict,
        naxis1=naxis1, naxis2=naxis2,
        beamsize_deg=beamsize_deg,
        ra_min=ra_min, ra_max=ra_max,
        dec_min=dec_min, dec_max=dec_max,
        kernel_fwhm_deg=kernel_fwhm_deg,
        need_weight_grid=need_weight_grid,
    )


def _grid_channel_worker(args):
    """Worker: grid one spectral channel.

    Returns (ichannel, 2D map or None, 2D weight-sum map or None). The weight
    map is only computed when the pool was initialised with
    need_weight_grid=True (i.e. a full per-channel weight cube was
    requested) — it grids the same per-channel weight array used for the
    data, with no further re-weighting, mirroring weight_map_2d below but
    per channel instead of once.
    """
    ichannel, channel_data, weights = args
    s = _gridding_worker_state
    valid = (np.isfinite(s['ras']) & np.isfinite(s['decs']) &
             np.isfinite(channel_data))
    if not np.any(valid):
        return ichannel, None, None
    w = weights[valid] if weights is not None else None
    channel_map = grid_to_map(
        s['ras'][valid], s['decs'][valid], channel_data[valid],
        wcs_header_dict=s['wcs_header_dict'],
        naxis1=s['naxis1'], naxis2=s['naxis2'],
        beamsize_deg=s['beamsize_deg'],
        ra_min=s['ra_min'], ra_max=s['ra_max'],
        dec_min=s['dec_min'], dec_max=s['dec_max'],
        weights=w,
        kernel_fwhm_deg=s.get('kernel_fwhm_deg'),
    )
    channel_weight_map = None
    if s.get('need_weight_grid') and w is not None:
        channel_weight_map = grid_to_map(
            s['ras'][valid], s['decs'][valid], w,
            wcs_header_dict=s['wcs_header_dict'],
            naxis1=s['naxis1'], naxis2=s['naxis2'],
            beamsize_deg=s['beamsize_deg'],
            ra_min=s['ra_min'], ra_max=s['ra_max'],
            dec_min=s['dec_min'], dec_max=s['dec_max'],
            weights=None,
            kernel_fwhm_deg=s.get('kernel_fwhm_deg'),
        )
    return ichannel, channel_map, channel_weight_map


# ---------------------------------------------------------------------------


def create_spectral_datacube(hdul: fits.HDUList,
                            beamsize_deg: float = 0.25,
                            pixsize: Optional[float] = None,
                            object_filter: Optional[str] = None,
                            figsize: Tuple[float, float] = (12, 5),
                            output_file: Optional[str] = None,
                            telescop: str = "",
                            n_jobs: int = -1,
                            weight_column: Optional[str] = None,
                            channel_weights: bool = False,
                            kernel_fwhm_arcsec: Optional[float] = None,
                            create_weights_datacube: bool = False) -> Tuple[np.ndarray, fits.Header, plt.Figure]:
    """
    Create a full 3D spectral datacube by gridding spectra across spatial and spectral axes.
    
    Creates a proper 3D datacube with WCS headers including a correctly-scaled velocity axis.
    The velocity axis is reconstructed from FITS table columns (VELOCITY, DELTAV, VELDEF, RESTFREQ)
    rather than simple channel indices.
    
    Parameters
    ----------
    hdul : astropy.io.fits.HDUList
        FITS file containing spectral observations.
    beamsize_deg : float, optional
        Beam size in degrees (FWHM). Default is 0.25°.
    pixsize : float, optional
        Map pixel size in degrees. If None, uses beamsize/3.
    object_filter : str, optional
        Filter observations by OBJECT column value.
    figsize : tuple, optional
        Figure size for plot (width, height). Default is (12, 5).
    output_file : str, optional
        If provided, saves the datacube as a FITS file.
    weight_column : str, optional
        Column name for per-spectrum weighting (e.g. 'RMSRATIOB').
        Gaussian weight: exp(-(value - 1.0)^2 / (2*0.5^2)).
    channel_weights : bool, optional
        If True, compute per-channel gridding weights from TSYS / TAU_SIG
        calibration spectra (via TSYS_INDEX / TAU_SIG_INDEX columns).
        Weight for spectrum i at channel c: exp(-tau_{i,c}) / tsys_{i,c}.
        Unlike create_integrated_map (which collapses to a scalar), this uses
        the channel-specific weight for each channel slice individually.
    create_weights_datacube : bool, optional
        If True *and* channel_weights is also True, write a full 3D
        WEIGHT_CUBE extension (nvel, dec, ra) holding the gridded weight sum
        per channel — needed because per-channel weights vary channel to
        channel, so a single 2D plane (as GILDAS's xy_map writes for its
        channel-independent per-spectrum weights) can't represent them.
        If channel_weights is False, this has no effect: the existing 2D
        WEIGHT_MAP/COVERAGE extensions are written exactly as before,
        whenever weight_column is set.

    Returns
    -------
    tuple
        - datacube (np.ndarray): 3D array of shape (nvel, dec, ra) containing gridded spectral data
        - wcs_header (astropy.io.fits.Header): WCS-compliant FITS header with proper velocity axis
        - fig (matplotlib.figure.Figure): Diagnostic plot showing sample slices
    
    Notes
    -----
    Creates a datacube by:
    1. Extracting SPECTRUM (shape: nobs x nvel), LONGITUDE, LATITUDE columns
    2. Extracting spectral axis parameters (VELOCITY, DELTAV, VELDEF, RESTFREQ)
    3. Filtering by object if specified
    4. Gridding each velocity channel separately onto a spatial map
    5. Stacking into a 3D cube (velocity, dec, ra) with proper WCS velocity axis
    
    The velocity axis is constructed from the FITS table columns:
    - VELOCITY (CRVAL3): Reference velocity in m/s (LSR convention)
    - DELTAV (CDELT3): Velocity spacing per channel in m/s
    - VELDEF: Velocity definition (typically 'RADI-LSR')
    - RESTFREQ: Rest frequency of the observed line in Hz
    
    The velocity of channel i is: v(i) = VELOCITY[0] + i * DELTAV[0]
    """
    if pixsize is None:
        pixsize = beamsize_deg / 3.0

    kernel_fwhm_deg = kernel_fwhm_arcsec / 3600.0 if kernel_fwhm_arcsec is not None else None

    # Validate output path before doing any expensive work
    if output_file:
        output_parent = Path(output_file).parent
        if not output_parent.exists():
            raise FileNotFoundError(
                f"Output directory does not exist: {output_parent}\n"
                f"  Please create it first or correct the output path."
            )

    # Extract data
    table = hdul[1].data
    spectra = table['SPECTRUM']  # Shape: (nobs, nvel)
    
    # Extract celestial coordinates (RA/Dec from CDELT2/3 and CRVAL2/3)
    ras, decs = _get_celestial_coords(hdul)
    
    nobs, nvel = spectra.shape

    def _obj_str(v):
        return (v.decode().strip() if isinstance(v, bytes) else str(v).strip()).upper()

    has_object_col = 'OBJECT' in table.names

    # Always strip calibration rows (TSYS, TAU_SIG) — these are not sky positions
    # and must never enter the gridder regardless of the object filter.
    if has_object_col:
        NON_SCIENCE = {'TSYS', 'TAU_SIG'}
        science_mask = np.array([_obj_str(o) not in NON_SCIENCE for o in table['OBJECT']])
        n_dropped = int(np.sum(~science_mask))
        if n_dropped:
            print(f"Excluded {n_dropped} calibration rows (TSYS / TAU_SIG)")
        spectra  = spectra[science_mask]
        ras      = ras[science_mask]
        decs     = decs[science_mask]
        remaining_objects = table['OBJECT'][science_mask]
        nobs     = len(ras)
    else:
        remaining_objects = None

    # Further filter by object name if requested
    if object_filter:
        if has_object_col:
            mask = np.array([
                object_filter.upper() in _obj_str(obj)
                for obj in remaining_objects
            ])
            if not np.any(mask):
                raise ValueError(f"No observations found matching object filter '{object_filter}'")
            spectra = spectra[mask]
            ras     = ras[mask]
            decs    = decs[mask]
            nobs    = len(ras)
            print(f"Filtered to {nobs} observations matching '{object_filter}'")
        else:
            warnings.warn("OBJECT column not found in FITS table - not filtering by object", UserWarning)
    else:
        print(f"Using {nobs} science observations (no object filter)")
    
    # Check if we have data
    if nobs == 0:
        raise ValueError("No observations available after filtering")
    if nvel == 0:
        raise ValueError("No spectral channels found")
    
    # Extract spectral axis parameters from FITS columns and headers
    velo_ref, deltav, restfreq_from_table, veldef, nchans_check, crpix1_spec = _get_spectral_axis_params(hdul)
    
    # Verify consistency
    if nchans_check != nvel:
        warnings.warn(f"Channel count mismatch: SPECTRUM has {nvel} channels, DELTAV/VELOCITY suggest {nchans_check}", UserWarning)
    
    print(f"Spectral axis parameters:")
    print(f"  Reference velocity: {velo_ref:.2f} m/s ({velo_ref/1000:.2f} km/s)")
    print(f"  Velocity step: {deltav:.2f} m/s ({deltav/1000:.4f} km/s)")
    print(f"  Reference pixel (CRPIX1): {crpix1_spec:.2f} (FITS 1-indexed)")
    print(f"  Reference channel: {crpix1_spec - 1:.2f} (0-indexed)")
    print(f"  Velocity definition: {veldef}")
    print(f"  Rest frequency: {restfreq_from_table:.10e} Hz")
    print(f"  Total velocity range: {deltav * nvel:.2f} m/s ({deltav * nvel / 1000:.2f} km/s)")
    
    # Compute velocity range
    v_at_ref_channel = velo_ref
    v_at_channel_0 = velo_ref - (crpix1_spec - 1) * deltav
    v_at_last_channel = velo_ref + (nvel - crpix1_spec) * deltav
    print(f"  Velocity at channel 0: {v_at_channel_0:.2f} m/s = {v_at_channel_0/1000:.2f} km/s")
    print(f"  Velocity at ref channel {crpix1_spec - 1:.0f}: {v_at_ref_channel:.2f} m/s = {v_at_ref_channel/1000:.2f} km/s")
    print(f"  Velocity at last channel {nvel-1}: {v_at_last_channel:.2f} m/s = {v_at_last_channel/1000:.2f} km/s")
    
    # Get spatial WCS parameters
    header = hdul[0].header if len(hdul) > 0 else hdul[1].header
    crval_ra = header.get('CRVAL2', 0.0)
    crval_dec = header.get('CRVAL3', 0.0)
    # Use the extracted velocity reference and frequency as the authoritative source
    velo_lsr = velo_ref
    restfreq = restfreq_from_table
    
    # Calculate map bounds
    ra_min, ra_max = ras.min(), ras.max()
    dec_min, dec_max = decs.min(), decs.max()
    ra_range = ra_max - ra_min
    dec_range = dec_max - dec_min
    
    # Add padding
    padding = 0.1
    ra_min -= ra_range * padding
    ra_max += ra_range * padding
    dec_min -= dec_range * padding
    dec_max += dec_range * padding
    
    # Calculate grid dimensions
    naxis1 = int(np.ceil((ra_max - ra_min) / pixsize))
    naxis2 = int(np.ceil((dec_max - dec_min) / pixsize))
    
    ra_center = (ra_min + ra_max) / 2
    dec_center = (dec_min + dec_max) / 2
    crpix1 = naxis1 / 2.0 + 1
    crpix2 = naxis2 / 2.0 + 1
    
    # Initialize datacube
    datacube = np.full((nvel, naxis2, naxis1), np.nan, dtype=np.float32)

    # --- Compute gridding weights ---
    # spectrum_weights: per-spectrum scalar (shape: nobs), or None
    spectrum_weights = None
    spectrum_weight_kind = None  # provenance: which transform was applied
    spectrum_weight_col = None   # provenance: actual column used (post-fallback)
    if weight_column is not None:
        col = weight_column
        if col not in table.names and col == 'RMSRATIOB' and 'RMSRATIO' in table.names:
            warnings.warn("RMSRATIOB not found, falling back to RMSRATIO for per-spectrum weighting.", UserWarning)
            col = 'RMSRATIO'
        if col not in table.names:
            raise ValueError(f"Weight column '{col}' not found in FITS file")
        wvals = np.asarray(table[col], dtype=np.float64)
        if has_object_col:
            wvals = wvals[science_mask]
        if object_filter and has_object_col:
            wvals = wvals[mask]
        spectrum_weight_col = col
        if 'RMSRATIO' in col.upper():
            # One-sided Gaussian roll-off: spectra at or below the theoretical
            # radiometer noise (RMSRATIO <= 1) are already as good as it gets and
            # keep full weight 1.0; only excess noise (RMSRATIO > 1) is penalised,
            # with a Gaussian falloff (0.61 at 1.5, 0.14 at 2.0). Deliberately NOT
            # symmetric about 1.0 — a quieter-than-theoretical spectrum must never
            # be down-weighted relative to an ideal one.
            sigma = 0.5
            excess = np.clip(wvals - 1.0, 0.0, None)
            spectrum_weights = np.where(
                np.isfinite(wvals),
                np.exp(-(excess ** 2) / (2 * sigma ** 2)),
                0.0,
            )
            spectrum_weight_kind = 'exp(-max(w-1,0)^2/(2*0.5^2))'
        else:
            # Generic column: inverse-variance weighting (w = 1/value²).
            # Non-finite, zero, or negative values get weight 0.
            spectrum_weights = np.where(
                np.isfinite(wvals) & (wvals > 0),
                1.0 / (wvals ** 2),
                0.0,
            )
            spectrum_weight_kind = '1/value^2'

    # channel_weight_matrix: per-spectrum per-channel (shape: nobs, nvel), or None
    # Unlike create_integrated_map (which collapses to scalar), we keep the full
    # 2D matrix so each channel slice gets its own optimal per-spectrum weights.
    channel_weight_matrix = None
    if channel_weights:
        has_cols = ('TSYS_INDEX' in table.names and 'TAU_SIG_INDEX' in table.names)
        if not has_cols:
            warnings.warn(
                "channel_weights=True requested but TSYS_INDEX / TAU_SIG_INDEX columns "
                "not found (available after prepare_for_pca). "
                "Falling back to uniform channel weights.", UserWarning)
        else:
            all_spectra_raw = table['SPECTRUM']
            tsys_idx = np.asarray(table['TSYS_INDEX'])
            tau_idx  = np.asarray(table['TAU_SIG_INDEX'])
            if has_object_col:
                tsys_idx = tsys_idx[science_mask]
                tau_idx  = tau_idx[science_mask]
            if object_filter and has_object_col:
                tsys_idx = tsys_idx[mask]
                tau_idx  = tau_idx[mask]
            channel_weight_matrix = np.ones((nobs, nvel), dtype=np.float64)
            for i in range(nobs):
                ti = int(tsys_idx[i])
                ai = int(tau_idx[i])
                if ti < 0 or ai < 0 or ti >= len(all_spectra_raw) or ai >= len(all_spectra_raw):
                    continue
                tsys_spec = np.asarray(all_spectra_raw[ti], dtype=np.float64)
                tau_spec  = np.asarray(all_spectra_raw[ai], dtype=np.float64)
                w_ch = np.where(
                    np.isfinite(tsys_spec) & np.isfinite(tau_spec) & (tsys_spec > 0),
                    (np.exp(-tau_spec) / tsys_spec) ** 2,
                    np.nan,
                )
                finite_mask = np.isfinite(w_ch)
                channel_weight_matrix[i, finite_mask] = w_ch[finite_mask]
                # Spectra with no valid calibration for a channel get weight=1.0 (neutral)
            # Normalize each channel by its median across spectra so the typical
            # spectrum has weight=1.0 per channel
            for c in range(nvel):
                col_w = channel_weight_matrix[:, c]
                finite_vals = col_w[np.isfinite(col_w) & (col_w > 0)]
                if len(finite_vals) > 0:
                    med = np.median(finite_vals)
                    if med > 0:
                        channel_weight_matrix[:, c] /= med

    # Build WCS header dict once (same geometry for every channel)
    wcs_header_dict_ch = create_wcs_header(
        naxis1=naxis1, naxis2=naxis2, naxis3=1,
        crval_ra=ra_center, crval_dec=dec_center,
        crpix1=crpix1, crpix2=crpix2, crpix3=crpix1_spec,
        cdelt1=-pixsize, cdelt2=pixsize, cdelt3=1.0,
        beamsize_deg=beamsize_deg,
        velo_lsr=velo_lsr,
        restfreq=restfreq,
        object_name=object_filter or "",
    )

    # Resolve n_jobs: -1 means use all available CPUs
    n_workers = os.cpu_count() if n_jobs == -1 else max(1, n_jobs)
    n_workers = min(n_workers, nvel)  # No point having more workers than channels

    # Effective output resolution after convolving the (already beam-sized)
    # data with the gridding kernel: sqrt(beam^2 + kernel_fwhm^2), same
    # quantity GILDAS's xy_map tracks as map%reso (vs. map%beam) and writes
    # to BMAJ/BMIN — see save_map_to_fits call below.
    _grid_kernel_fwhm_deg = kernel_fwhm_deg if kernel_fwhm_deg is not None else beamsize_deg
    effective_resolution_deg = float(np.sqrt(beamsize_deg**2 + _grid_kernel_fwhm_deg**2))
    if kernel_fwhm_deg is not None:
        print(f"Gridding kernel FWHM: {kernel_fwhm_arcsec:.1f}\" "
              f"(beam: {beamsize_deg * 3600:.1f}\", effective resolution: "
              f"{effective_resolution_deg * 3600:.1f}\")")
    else:
        print(f"Gridding kernel FWHM: {_grid_kernel_fwhm_deg * 3600:.1f}\" (= beam; effective "
              f"resolution: {effective_resolution_deg * 3600:.1f}\". Pass --kernel-fwhm "
              f"<smaller value> for sharper-but-noisier resolution-limited cases.)")
    print(f"Creating datacube: {nvel} channels × {naxis2} × {naxis1} pixels "
          f"(n_jobs={n_workers})...")

    # A full 3D weight cube only makes sense when weights actually vary by
    # channel (channel_weights=True). Otherwise the per-spectrum weight is
    # identical for every channel, so the existing 2D WEIGHT_MAP below
    # already captures it — same reasoning GILDAS's xy_map uses to force its
    # own .wei file to 2D (see xymap.f90: hwei%gil%ndim = 2).
    need_weight_cube = bool(channel_weights and create_weights_datacube)
    weight_cube_3d = (np.full((nvel, naxis2, naxis1), np.nan, dtype=np.float32)
                       if need_weight_cube else None)

    def _get_channel_weights(ichannel):
        """Combine spectrum and channel weights for a single channel slice."""
        sw = spectrum_weights  # per-spectrum scalar (nobs,) or None
        cw = channel_weight_matrix[:, ichannel] if channel_weight_matrix is not None else None
        if sw is not None and cw is not None:
            return sw * cw
        elif cw is not None:
            return cw
        else:
            return sw  # may be None → uniform

    channel_args = [
        (ichannel,
         np.array(spectra[:, ichannel], dtype=np.float32),
         _get_channel_weights(ichannel))
        for ichannel in range(nvel)
    ]

    if n_workers == 1:
        # Sequential path — useful for debugging
        for ichannel, channel_data, weights in channel_args:
            if ichannel % max(1, nvel // 10) == 0:
                print(f"  Channel {ichannel+1}/{nvel}...", end='\r')
            valid = np.isfinite(ras) & np.isfinite(decs) & np.isfinite(channel_data)
            if not np.any(valid):
                continue
            w = weights[valid] if weights is not None else None
            datacube[ichannel] = grid_to_map(
                ras[valid], decs[valid], channel_data[valid],
                wcs_header_dict=wcs_header_dict_ch,
                naxis1=naxis1, naxis2=naxis2,
                beamsize_deg=beamsize_deg,
                ra_min=ra_min, ra_max=ra_max,
                dec_min=dec_min, dec_max=dec_max,
                weights=w,
                kernel_fwhm_deg=kernel_fwhm_deg,
            )
            if need_weight_cube and w is not None:
                weight_cube_3d[ichannel] = grid_to_map(
                    ras[valid], decs[valid], w,
                    wcs_header_dict=wcs_header_dict_ch,
                    naxis1=naxis1, naxis2=naxis2,
                    beamsize_deg=beamsize_deg,
                    ra_min=ra_min, ra_max=ra_max,
                    dec_min=dec_min, dec_max=dec_max,
                    weights=None,
                    kernel_fwhm_deg=kernel_fwhm_deg,
                )
    else:
        from multiprocessing import Pool
        with Pool(
            processes=n_workers,
            initializer=_init_gridding_worker,
            initargs=(
                np.array(ras, dtype=np.float32),
                np.array(decs, dtype=np.float32),
                wcs_header_dict_ch,
                naxis1, naxis2, beamsize_deg,
                ra_min, ra_max, dec_min, dec_max,
                kernel_fwhm_deg,
                need_weight_cube,
            ),
        ) as pool:
            n_done = 0
            for ichannel, channel_map, channel_weight_map in pool.imap_unordered(
                    _grid_channel_worker, channel_args):
                if channel_map is not None:
                    datacube[ichannel] = channel_map
                if need_weight_cube and channel_weight_map is not None:
                    weight_cube_3d[ichannel] = channel_weight_map
                n_done += 1
                if n_done % max(1, nvel // 10) == 0:
                    print(f"  {n_done}/{nvel} channels done...", end='\r')

    print(f"  {nvel}/{nvel} channels done.    ")

    # Compute 2D coverage map (sum of kernel weights per pixel, independent of channel)
    coverage_map_2d = None
    if HAS_CYGRID:
        # Must match grid_to_map's kernel exactly (beam default, 3*FWHM support).
        effective_fwhm = kernel_fwhm_deg if kernel_fwhm_deg is not None else beamsize_deg
        _ks = effective_fwhm / 2.355
        _cov_gridder = cygrid.WcsGrid(wcs_header_dict_ch)
        _cov_gridder.set_kernel('gauss1d', (_ks,), 3.0 * effective_fwhm, _ks / 2.0)
        _valid = np.isfinite(ras) & np.isfinite(decs)
        _cov_gridder.grid(
            np.array(ras[_valid], dtype=np.float64),
            np.array(decs[_valid], dtype=np.float64),
            np.ones(int(_valid.sum()), dtype=np.float64),
        )
        # Use the *summed* kernel weights (Σ K), NOT get_datacube(). Gridding a
        # field of ones and reading get_datacube() returns the kernel-weighted
        # AVERAGE of ones = 1.0 in every covered pixel, which carries no coverage
        # gradient (so --coverage-threshold could never fire). get_weights()
        # returns Σ K per pixel — the true sample/coverage density that tapers
        # toward the map edges.
        coverage_map_2d = _cov_gridder.get_weights().squeeze().astype(np.float32)
        print(f"  Coverage map (Σ kernel weights): peak={coverage_map_2d.max():.3f}, "
              f"pixels with coverage>0: {(coverage_map_2d > 0).sum()}")

    # Compute 2D weight map: grid spectrum_weights with the same kernel so each
    # pixel shows the kernel-weighted sum of spectrum weights (denominator of the
    # weighted average).  Only computed when per-spectrum weights are present,
    # and skipped when the full 3D weight cube was built instead (that case
    # already captures everything this 2D map would show, plus the per-channel
    # variation it can't represent).
    weight_map_2d = None
    if HAS_CYGRID and spectrum_weights is not None and not need_weight_cube:
        # Must match grid_to_map's kernel exactly (beam default, 3*FWHM support).
        _weff = kernel_fwhm_deg if kernel_fwhm_deg is not None else beamsize_deg
        _wks = _weff / 2.355
        _wm_gridder = cygrid.WcsGrid(wcs_header_dict_ch)
        _wm_gridder.set_kernel('gauss1d', (_wks,), 3.0 * _weff, _wks / 2.0)
        _valid_w = np.isfinite(ras) & np.isfinite(decs)
        _wm_gridder.grid(
            np.array(ras[_valid_w], dtype=np.float64),
            np.array(decs[_valid_w], dtype=np.float64),
            np.array(spectrum_weights[_valid_w], dtype=np.float64),
        )
        # get_unweighted_datacube() = Σ K·wᵢ (the actual denominator of the
        # science weighted-average), NOT get_datacube() which would return the
        # normalized MEAN weight Σ K·wᵢ / Σ K (≈1.0 everywhere when the input
        # weights are uniform, hiding all spatial structure).
        weight_map_2d = _wm_gridder.get_unweighted_datacube().squeeze().astype(np.float32)
        if coverage_map_2d is not None:
            weight_map_2d[coverage_map_2d <= 0] = np.nan
        print(f"  Weight map ({weight_column}, Σ kernel·weight): peak={np.nanmax(weight_map_2d):.4f}, "
              f"pixels with weight>0: {np.isfinite(weight_map_2d).sum()}")

    # Create WCS header for 3D cube (spectral, dec, ra)
    # Note: CDELT3 should be deltav (velocity spacing per channel in m/s)
    # CRVAL3 is the reference velocity in m/s
    # CRPIX3 is the reference pixel from FITS header (usually around 505-506)
    wcs_header = create_wcs_header(
        naxis1=naxis1, naxis2=naxis2, naxis3=nvel,
        crval_ra=ra_center, crval_dec=dec_center,
        crpix1=crpix1, crpix2=crpix2, crpix3=crpix1_spec,
        cdelt1=-pixsize, cdelt2=pixsize, cdelt3=deltav,  # Use actual velocity spacing!
        beamsize_deg=beamsize_deg,
        velo_lsr=velo_lsr,
        restfreq=restfreq,
        object_name=object_filter or "",
        telescop=telescop,
    )

    # --- Gridding-provenance cards: make the cube self-documenting about how it
    # was weighted, so channel-weighted and plain cubes can be told apart from
    # the header alone (not just the filename). ---
    chan_applied = bool(channel_weights and channel_weight_matrix is not None)
    grid_method = 'cygrid' if HAS_CYGRID else 'scipy'
    wcs_header.set('GRIDMETH', grid_method, 'gridding engine (cygrid or scipy)')
    wcs_header.set('CHANWEI', chan_applied,
                   'per-channel weighting actually applied')
    if channel_weights and not chan_applied:
        # requested but fell back to uniform (missing TSYS_INDEX/TAU_SIG_INDEX) —
        # exactly the silent no-op we want on the record.
        wcs_header.set('CHWEIREQ', True,
                       'per-channel weighting requested but fell back to uniform')
    if chan_applied:
        wcs_header.set('CHWEIFRM', '(exp(-tau)/Tsys)^2',
                       'per-channel weight formula (median-normalised)')
    wcs_header.set('SPECWEI', spectrum_weight_col or 'NONE',
                   'per-spectrum weight column')
    if spectrum_weight_kind is not None:
        wcs_header.set('SPECWFRM', spectrum_weight_kind,
                       'per-spectrum weight transform')
    wcs_header.set('WEICUBE', bool(need_weight_cube),
                   '3D WEIGHT_CUBE HDU present')

    # Create visualization of a sample slice
    mid_channel = nvel // 2
    fig, axes = plt.subplots(1, 3, figsize=figsize)
    
    # Format RA axis in hours:minutes:seconds and Dec axis in degrees:arcminutes:arcseconds
    from matplotlib.ticker import FuncFormatter, MaxNLocator
    ra_formatter = FuncFormatter(lambda x, p: format_ra_hms(x))
    dec_formatter = FuncFormatter(lambda x, p: format_dec_dms(x))
    
    # First channel
    datacube_display_1 = np.fliplr(datacube[0])
    im1 = axes[0].imshow(datacube_display_1, origin='lower', extent=[ra_min, ra_max, dec_min, dec_max],
                        cmap='viridis', aspect='auto')
    axes[0].set_title(f'Channel 1 (first)')
    axes[0].xaxis.set_major_formatter(ra_formatter)
    axes[0].yaxis.set_major_formatter(dec_formatter)
    axes[0].xaxis.set_major_locator(MaxNLocator(nbins=2, integer=False))
    axes[0].invert_xaxis()  # RA increases from right to left
    axes[0].set_xlabel('RA (h:m:s)')
    axes[0].set_ylabel('Dec (°:′:″)')
    plt.colorbar(im1, ax=axes[0])
    
    # Middle channel
    datacube_display_2 = np.fliplr(datacube[mid_channel])
    im2 = axes[1].imshow(datacube_display_2, origin='lower', extent=[ra_min, ra_max, dec_min, dec_max],
                        cmap='viridis', aspect='auto')
    axes[1].set_title(f'Channel {mid_channel+1} (middle)')
    axes[1].xaxis.set_major_formatter(ra_formatter)
    axes[1].yaxis.set_major_formatter(dec_formatter)
    axes[1].xaxis.set_major_locator(MaxNLocator(nbins=2, integer=False))
    axes[1].invert_xaxis()  # RA increases from right to left
    axes[1].set_xlabel('RA (h:m:s)')
    axes[1].set_ylabel('Dec (°:′:″)')
    plt.colorbar(im2, ax=axes[1])
    
    # Last channel
    datacube_display_3 = np.fliplr(datacube[-1])
    im3 = axes[2].imshow(datacube_display_3, origin='lower', extent=[ra_min, ra_max, dec_min, dec_max],
                        cmap='viridis', aspect='auto')
    axes[2].set_title(f'Channel {nvel} (last)')
    axes[2].xaxis.set_major_formatter(ra_formatter)
    axes[2].yaxis.set_major_formatter(dec_formatter)
    axes[2].xaxis.set_major_locator(MaxNLocator(nbins=2, integer=False))
    axes[2].invert_xaxis()  # RA increases from right to left
    axes[2].set_xlabel('RA (h:m:s)')
    axes[2].set_ylabel('Dec (°:′:″)')
    plt.colorbar(im3, ax=axes[2])
    
    fig.suptitle(f'3D Spectral Datacube: {nvel} channels × {naxis2} × {naxis1} pixels, beam={beamsize_deg*3600:.1f}″',
                fontsize=13, fontweight='bold')
    fig.tight_layout()
    
    # Save to FITS if requested
    if output_file:
        # Prepare spectral axis parameters for the FITS file
        spectral_params = {
            'nvel': nvel,
            'velo_ref': velo_ref,
            'crpix1_spec': crpix1_spec,
            'deltav': deltav,
            'restfreq': restfreq_from_table,
            'veldef': veldef,
        }
        save_map_to_fits(datacube, wcs_header, output_file, beam_maj_deg=effective_resolution_deg,
                        spectral_params=spectral_params,
                        coverage_map=coverage_map_2d,
                        weight_map=weight_map_2d,
                        weight_cube=weight_cube_3d,
                        weight_column=weight_column)

    return datacube, wcs_header, fig


# ---------------------------------------------------------------------------
# Peak-aligned stacking — visual aid for averaging spectra whose line centre
# genuinely differs from spaxel to spaxel (e.g. a rotation curve), without
# that velocity spread itself broadening the averaged line shape.
# ---------------------------------------------------------------------------

def _peak_aligned_mean_spectrum(spectra_2d: np.ndarray, channels_kms: np.ndarray,
                                search_range: Optional[Tuple[float, float]] = None) -> np.ndarray:
    """
    Average a set of spectra after aligning each one's own peak to a common velocity.

    For each input spectrum, the peak channel is located (within `search_range`
    if given) and refined to sub-channel precision with a 3-point parabolic fit.
    Every spectrum is then shifted, by linear interpolation, so its own peak
    lands on the median peak velocity across all of them, before averaging.

    This is purely a visual aid for inspecting line *shape* (width, asymmetry,
    wings) from a stack of spectra that legitimately have different centroid
    velocities — e.g. spaxels spanning a galaxy's rotation curve — where a
    plain mean would smear the profile out by that velocity spread alone.
    It is not appropriate for anything that needs the true per-spaxel velocity
    preserved (use the un-aligned mean for that).

    Parameters
    ----------
    spectra_2d : np.ndarray, shape (n_spectra, n_channels)
        Spectra to align and average — already restricted to whichever pixels
        are being stacked (zoom box, circular aperture, or the whole map).
    channels_kms : np.ndarray, shape (n_channels,)
        Velocity axis, assumed evenly spaced.
    search_range : (float, float), optional
        Velocity range to search for the peak in (e.g. the same window used
        for the moment-0 integration). Defaults to the full axis, which risks
        picking a noise spike in low-S/N spectra.

    Returns
    -------
    np.ndarray, shape (n_channels,)
        The peak-aligned mean spectrum.
    """
    spectra_2d = np.atleast_2d(spectra_2d)
    n_spec, n_chan = spectra_2d.shape

    if search_range is not None:
        v_lo, v_hi = sorted(search_range)
        search_idx = np.where((channels_kms >= v_lo) & (channels_kms <= v_hi))[0]
    else:
        search_idx = np.arange(n_chan)

    if len(search_idx) < 3:
        # Not enough channels to refine a peak — fall back to a plain mean.
        return np.nanmean(spectra_2d, axis=0)

    dv = float(channels_kms[1] - channels_kms[0])

    peak_velocities = np.full(n_spec, np.nan)
    for i in range(n_spec):
        row = spectra_2d[i]
        sub = row[search_idx]
        if np.all(np.isnan(sub)):
            continue
        k = search_idx[int(np.nanargmax(sub))]
        # 3-point parabolic refinement around the peak for sub-channel precision.
        k_lo, k_hi = max(k - 1, 0), min(k + 1, n_chan - 1)
        y_lo, y0, y_hi = row[k_lo], row[k], row[k_hi]
        denom = y_lo - 2 * y0 + y_hi
        if k_lo < k < k_hi and np.isfinite(denom) and denom != 0:
            delta = float(np.clip(0.5 * (y_lo - y_hi) / denom, -1.0, 1.0))
        else:
            delta = 0.0
        peak_velocities[i] = channels_kms[k] + delta * dv

    valid = np.isfinite(peak_velocities)
    if not np.any(valid):
        return np.nanmean(spectra_2d, axis=0)
    reference_v = float(np.nanmedian(peak_velocities[valid]))

    chan_idx = np.arange(n_chan)
    shifted = np.empty_like(spectra_2d, dtype=float)
    for i in range(n_spec):
        if not valid[i]:
            shifted[i] = spectra_2d[i]
            continue
        shift_chan = (reference_v - peak_velocities[i]) / dv
        shifted[i] = np.interp(chan_idx - shift_chan, chan_idx, spectra_2d[i],
                               left=np.nan, right=np.nan)

    return np.nanmean(shifted, axis=0)


# ---------------------------------------------------------------------------
# Interactive polygon selection for restricting the rendered area
# ---------------------------------------------------------------------------

def _interactive_polygon(display_map, colormap='rainbow', title=None):
    """Let the user draw a polygon on a moment-0 map and return its vertices.

    Interaction (as requested):
      * MIDDLE mouse button  — drop a polygon vertex
      * LEFT   mouse button  — close the polygon (finish; needs >= 3 vertices)
      * RIGHT  mouse button  — undo the last vertex

    Returns a list of (x, y) pixel-coordinate vertices (x = column, y = row),
    or None if fewer than three vertices were placed. The window is drawn in
    plain pixel space (origin='lower'), so the returned coordinates map directly
    onto array indices.
    """
    import matplotlib.pyplot as plt

    finite = np.isfinite(display_map)
    if finite.any():
        vmin, vmax = np.nanpercentile(display_map[finite], [1, 99])
        if not np.isfinite(vmin) or not np.isfinite(vmax) or vmin == vmax:
            vmin, vmax = None, None
    else:
        vmin, vmax = None, None

    fig, ax = plt.subplots(figsize=(9, 8))
    ax.imshow(display_map, origin='lower', cmap=colormap, vmin=vmin, vmax=vmax)
    ax.set_xlabel('X pixel')
    ax.set_ylabel('Y pixel')
    ax.set_title(title or
                 'Middle-click: add vertex   |   Left-click: close   |   '
                 'Right-click: undo')

    verts = []
    state = {'closed': False}
    line, = ax.plot([], [], '-o', color='magenta', lw=1.5, ms=5, mfc='yellow')

    def _redraw():
        if verts:
            xs = [v[0] for v in verts]
            ys = [v[1] for v in verts]
            if state['closed']:
                xs = xs + [xs[0]]
                ys = ys + [ys[0]]
            line.set_data(xs, ys)
        else:
            line.set_data([], [])
        fig.canvas.draw_idle()

    def _on_click(event):
        # Ignore clicks outside the axes or while a toolbar tool (zoom/pan) is active.
        if event.inaxes is not ax or event.xdata is None:
            return
        tb = getattr(fig.canvas, 'toolbar', None)
        if tb is not None and getattr(tb, 'mode', ''):
            return
        if event.button == 2:            # middle → add vertex
            verts.append((float(event.xdata), float(event.ydata)))
            _redraw()
        elif event.button == 3:          # right → undo
            if verts:
                verts.pop()
                _redraw()
        elif event.button == 1:          # left → close/finish
            if len(verts) >= 3:
                state['closed'] = True
                _redraw()
                plt.close(fig)

    cid = fig.canvas.mpl_connect('button_press_event', _on_click)
    plt.show()
    fig.canvas.mpl_disconnect(cid)

    if len(verts) >= 3:
        return verts
    return None


def _polygon_mask(vertices, ny, nx):
    """Boolean (ny, nx) mask that is True OUTSIDE the given polygon.

    Vertices are (x, y) pixel coordinates (x = column, y = row); pixel centres
    are tested for containment.
    """
    from matplotlib.path import Path as _Path
    path = _Path(np.asarray(vertices, dtype=float))
    yy, xx = np.mgrid[0:ny, 0:nx]
    pts = np.column_stack([xx.ravel(), yy.ravel()])
    inside = path.contains_points(pts).reshape(ny, nx)
    return ~inside


# ---------------------------------------------------------------------------
# Collapse 3D datacube to 2D integrated map
# ---------------------------------------------------------------------------

def _spectral_smooth_nanaware(cube: np.ndarray, sigma_chan: float) -> np.ndarray:
    """NaN-aware Gaussian smoothing of a cube along the velocity (axis-0) direction.

    Convolves each spaxel's spectrum with a 1D Gaussian of standard deviation
    ``sigma_chan`` channels, ignoring NaNs (each output channel is renormalised by
    the summed kernel weight of the finite input channels, so map edges and blanked
    voxels do not bleed zeros in). Vectorised over the spatial axes; loops only over
    the ~6·sigma+1 kernel taps. Uses no scipy (see the spatial-smoothing path, which
    likewise uses astropy) — the kernel is built explicitly here.

    Parameters
    ----------
    cube : (nvel, ny, nx) ndarray
    sigma_chan : float
        Gaussian standard deviation in channels.

    Returns
    -------
    (nvel, ny, nx) ndarray of float64; channels with no finite support are NaN.
    """
    if sigma_chan <= 0:
        return cube
    half = int(np.ceil(3.0 * sigma_chan))
    taps = np.arange(-half, half + 1)
    kernel = np.exp(-0.5 * (taps / sigma_chan) ** 2)
    kernel /= kernel.sum()

    finite = np.isfinite(cube)
    c0 = np.where(finite, cube, 0.0).astype(np.float64)
    wv = finite.astype(np.float64)
    num = np.zeros_like(c0)
    den = np.zeros_like(c0)
    nvel = cube.shape[0]
    for t, wgt in zip(taps, kernel):
        if t == 0:
            num += wgt * c0
            den += wgt * wv
        elif t > 0:
            # output channel i pulls from input channel i+t (edges lose support → lower den)
            num[:nvel - t] += wgt * c0[t:]
            den[:nvel - t] += wgt * wv[t:]
        else:  # t < 0
            num[-t:] += wgt * c0[:nvel + t]
            den[-t:] += wgt * wv[:nvel + t]
    with np.errstate(invalid='ignore', divide='ignore'):
        out = np.where(den > 0, num / den, np.nan)
    return out


def _velocity_field_from_cube(cube, channels_kms, line_mask, deltav_kms,
                              snr_min, smooth_pix):
    """Build a smooth per-pixel line-velocity field v0(x,y) for shuffling.

    Computes the intensity-weighted mean velocity (moment-1) within ``line_mask``,
    keeps only pixels above ``snr_min`` (so faint/noise pixels don't contribute a
    noisy centroid), then spatially smooths and gap-fills so a value is defined at
    EVERY pixel — faint arm pixels inherit v0 from their bright neighbours. That
    robust, smooth field is what lets the subsequent shuffle centre the line
    correctly at low S/N without a per-pixel self-peak (unlike peak-range-int).

    Returns (v0_raw, v0_field, snr_map): v0_raw = moment-1 on the good pixels only
    (NaN elsewhere, for inspection); v0_field = smoothed+filled field used to shuffle.
    """
    from astropy.convolution import Gaussian2DKernel, convolve as _convolve

    line = cube[line_mask, :, :]                       # (nline, ny, nx)
    nline = int(line_mask.sum())
    off = cube[~line_mask, :, :]
    with np.errstate(invalid='ignore'):
        sigma = np.nanstd(off, axis=0) if off.shape[0] >= 2 else np.full(cube.shape[1:], np.nan)

    v = channels_kms[line_mask]
    I = np.where(np.isfinite(line), line, 0.0)
    Isum = I.sum(axis=0)                                # (ny, nx)
    mom0 = Isum * deltav_kms
    with np.errstate(invalid='ignore', divide='ignore'):
        mom1 = np.where(Isum != 0, (v[:, None, None] * I).sum(axis=0) / Isum, np.nan)
        snr_map = np.where((sigma > 0) & np.isfinite(sigma),
                           mom0 / (sigma * np.sqrt(nline) * deltav_kms), np.nan)

    good = np.isfinite(mom1) & np.isfinite(snr_map) & (snr_map >= snr_min)
    v0_raw = np.where(good, mom1, np.nan)

    # Smooth + gap-fill: interpolate a smooth field into every pixel.
    kernel = Gaussian2DKernel(x_stddev=smooth_pix)
    v0_field = _convolve(v0_raw, kernel, nan_treatment='interpolate',
                         preserve_nan=False, boundary='extend')
    # Any pixel still without support falls back to the global median velocity
    # (near-systemic → ~no shift; these are empty/edge pixels anyway).
    if np.any(np.isfinite(v0_raw)):
        v0_field = np.where(np.isfinite(v0_field), v0_field, np.nanmedian(v0_raw))
    return v0_raw, v0_field.astype(np.float64), snr_map


def _shuffle_cube(cube, channels_kms, v0_field, ref_v):
    """Shift every spaxel's spectrum so its line at v0(x,y) lands at ``ref_v``.

    new[v] = old[v + (v0 - ref_v)]  →  the feature at v0 appears at ref_v.
    Vectorised fractional-channel shift via linear interpolation between the two
    bracketing integer shifts (no scipy). Channels shifted past the band edge, and
    NaN source voxels, become NaN in the aligned cube.
    """
    nvel, ny, nx = cube.shape
    dv_per_index = channels_kms[1] - channels_kms[0]    # signed km/s per channel
    shift_ch = (v0_field - ref_v) / dv_per_index         # (ny, nx), signed, fractional
    i0 = np.floor(shift_ch).astype(np.int64)             # (ny, nx)
    w = (shift_ch - i0)[None, :, :]                      # (1, ny, nx) fractional part

    idx = np.arange(nvel)[:, None, None]
    src0 = idx + i0[None, :, :]                          # (nvel, ny, nx)
    src1 = src0 + 1
    valid0 = (src0 >= 0) & (src0 < nvel)
    valid1 = (src1 >= 0) & (src1 < nvel)
    v0v = np.take_along_axis(cube, np.clip(src0, 0, nvel - 1), axis=0)
    v1v = np.take_along_axis(cube, np.clip(src1, 0, nvel - 1), axis=0)
    v0v = np.where(valid0, v0v, np.nan)
    v1v = np.where(valid1, v1v, np.nan)
    out = (1.0 - w) * v0v + w * v1v
    # If one neighbour is out of range but the other is in, keep the in-range one.
    out = np.where(valid0 & ~valid1, v0v, out)
    out = np.where(valid1 & ~valid0, v1v, out)
    return out.astype(cube.dtype)


def collapse_cube(
    cube_fits: str,
    velocity_range: Optional[Tuple[float, float]] = None,
    zoom_size_arcmin: Optional[float] = None,
    zoom_x: Optional[float] = None,
    zoom_y: Optional[float] = None,
    zoom_ra: Optional[float] = None,
    zoom_dec: Optional[float] = None,
    region_x: Optional[float] = None,
    region_y: Optional[float] = None,
    region_ra: Optional[float] = None,
    region_dec: Optional[float] = None,
    region_radius_arcmin: Optional[float] = None,
    use_wcs: bool = False,
    fits_output: Optional[str] = None,
    noise_map_output: Optional[str] = None,
    plot_output: Optional[str] = None,
    overwrite: bool = True,
    colormap: str = 'rainbow',
    coverage_threshold: float = 0.0,
    mask_ra: Optional[float] = None,
    mask_dec: Optional[float] = None,
    mask_radius_arcmin: Optional[float] = None,
    trim_edges: Optional[str] = None,
    suppress_negative: bool = False,
    suppress_high: Optional[float] = None,
    mode: str = 'moment0',
    peak_range_channels: int = 5,
    peak_smooth_channels: Optional[float] = None,
    shuffle: bool = False,
    shuffle_snr: float = 5.0,
    shuffle_field_smooth: float = 6.0,
    shuffle_window_kms: float = 20.0,
    shuffle_ref_velocity: Optional[float] = None,
    shuffle_field_output: Optional[str] = None,
    hex_plot: bool = False,
    contour: bool = False,
    stretch: str = 'linear',
    gamma: float = 1.0,
    smooth_sigma: Optional[float] = None,
    interpolation: str = 'nearest',
    percentile_clip: Optional[Tuple[float, float]] = None,
    snr_threshold: Optional[float] = None,
    drop_rms: Optional[float] = None,
    drop_window: Optional[Tuple[float, float]] = None,
    align_peaks: bool = False,
    collapse_weights: bool = False,
    select_polygon: bool = False,
    polygon: Optional[List[Tuple[float, float]]] = None,
) -> Tuple[np.ndarray, fits.Header, plt.Figure]:
    """
    Collapse a 3D spectral datacube to a 2D integrated intensity map (moment-0).

    Reads an already-gridded FITS cube (NAXIS1=RA, NAXIS2=Dec, NAXIS3=velocity),
    sums along the velocity axis within an optional velocity range, and returns
    the 2D map with a corresponding 2D WCS header.

    Weights are not applied here — they were already baked into the cube during
    gridding.  The resulting map has units of K km/s.

    With ``collapse_weights=True`` the function instead reads the 3D
    ``WEIGHT_CUBE`` extension (written by ``create_datacube --weight-channels
    --create-weights-datacube``) and integrates the per-channel gridding weights
    (``Σ K·w``) over the velocity window into a 2D weight/sensitivity map — a
    plain sum, NOT multiplied by channel width (units: ``weight``).  All the
    spatial/display options (zoom, region spectrum, circular mask, trim-edges,
    suppress-high, stretch/gamma, smooth, percentile-clip, hex, contour, WCS
    axes, polygon, colormap) apply unchanged; the line/noise-specific options
    (noise map, ``snr_threshold``, ``drop_rms``, ``shuffle``, ``align_peaks``,
    peak modes) are meaningless for a weight cube and are ignored with a warning.

    Parameters
    ----------
    cube_fits : str
        Path to the 3D FITS datacube.
    velocity_range : (float, float), optional
        Integration range in km/s (v_min, v_max).  If None, all channels are used.
    zoom_size_arcmin : float, optional
        Side length in arcmin of the zoom region shown in the plot.
        The zoom panel uses a colour scale recalculated for that region only.
    zoom_x : float, optional
        X pixel coordinate of the zoom region centre.  Defaults to map centre.
        Takes priority over zoom_ra/zoom_dec.
    zoom_y : float, optional
        Y pixel coordinate of the zoom region centre.
    zoom_ra : float, optional
        RA in degrees of the zoom region centre.
    zoom_dec : float, optional
        Dec in degrees of the zoom region centre.
    region_x : float, optional
        X pixel coordinate (as shown in the plot) of the extraction aperture
        centre.  Takes priority over region_ra/region_dec if both are given.
        If None but region_radius_arcmin is given, map centre is used.
    region_y : float, optional
        Y pixel coordinate of the extraction aperture centre.
    region_ra : float, optional
        RA in degrees of the extraction aperture centre (used when region_x/y
        are not given).
    region_dec : float, optional
        Dec in degrees of the extraction aperture centre.
    region_radius_arcmin : float, optional
        Radius in arcmin of a circular aperture for spectrum extraction.
        Triggers an extra spectrum panel and a circle overlay on the map(s).
    align_peaks : bool, optional
        Before averaging spectra into the zoom/region/full-map spectrum
        panel(s), shift each one (sub-channel, by interpolation) so its own
        peak — found within `velocity_range` if given, else the full axis —
        lands on the median peak velocity of the stack. A visual aid for
        inspecting line shape from spaxels whose true centroid velocity
        differs (e.g. a rotation curve), without that spread itself
        broadening the averaged profile. Default False (plain mean).
    use_wcs : bool
        If True, display map axes in RA/Dec instead of pixel indices.
        Default False.
    fits_output : str, optional
        If given, save the 2D map as a FITS file at this path.
    plot_output : str, optional
        If given, save the figure to this path.
    overwrite : bool
        Overwrite existing output files.  Default True.

    Returns
    -------
    collapsed : np.ndarray  shape (ny, nx)
    header2d  : fits.Header  (2D WCS)
    fig       : matplotlib.Figure
    """
    # ------------------------------------------------------------------
    # 1. Load cube
    # ------------------------------------------------------------------
    weight_map_only = False
    with fits.open(cube_fits) as hdul:
        if collapse_weights:
            ext_names = [h.name for h in hdul]
            if 'WEIGHT_CUBE' in ext_names:
                cube = hdul['WEIGHT_CUBE'].data.astype(float)   # (nvel, ny, nx)
            elif 'WEIGHT_MAP' in ext_names:
                # Per-spectrum weighting (no per-channel weights) records a single
                # 2D WEIGHT_MAP (Σ K·w per pixel) — the analog of CLASS's .wei
                # file. Load it as a 1-channel cube so the shared collapse/display
                # machinery renders it unchanged (moment0 of one plane = the plane).
                cube = hdul['WEIGHT_MAP'].data.astype(float)[np.newaxis, :, :]
                weight_map_only = True
            else:
                raise ValueError(
                    f"No WEIGHT_CUBE or WEIGHT_MAP extension in {cube_fits} "
                    f"(found: {', '.join(n for n in ext_names if n) or 'none'}). "
                    f"Re-run create_datacube with --weight-spectra (writes a 2D "
                    f"WEIGHT_MAP) or --weight-channels --create-weights-datacube "
                    f"(writes a 3D WEIGHT_CUBE) to record the weights used."
                )
        else:
            cube = hdul[0].data.astype(float)   # (nvel, ny, nx)
        header = hdul[0].header.copy()          # WCS lives in the primary header

    if cube.ndim != 3:
        raise ValueError(f"Expected a 3-D datacube, got shape {cube.shape}")

    nvel, ny, nx = cube.shape

    # Weight-collapse ignores options defined in terms of a spectral line + noise
    # (a weight cube has neither).  Disable them here — warning if the user set
    # any — so the rest of the (shared) collapse/display machinery runs unchanged.
    if collapse_weights:
        _incompat = []
        if shuffle:              _incompat.append('--shuffle');          shuffle = False
        if snr_threshold is not None: _incompat.append('--snr-threshold'); snr_threshold = None
        if drop_rms is not None: _incompat.append('--drop-rms');         drop_rms = None
        if align_peaks:          _incompat.append('--align-peaks');      align_peaks = False
        if noise_map_output:     _incompat.append('--output-noise-map'); noise_map_output = None
        if mode != 'moment0':    _incompat.append(f'--mode {mode}');     mode = 'moment0'
        # A 2D WEIGHT_MAP (loaded as a single channel) has no spectral axis to
        # window over, so a velocity range would blank it entirely.
        if weight_map_only and velocity_range is not None:
            _incompat.append('--velocity-range'); velocity_range = None
        if _incompat:
            warnings.warn(
                "collapse_weights: ignoring line/noise-specific option(s) "
                + ", ".join(_incompat)
                + " — they have no meaning for a weight cube.", UserWarning)

    # ------------------------------------------------------------------
    # 2. Reconstruct velocity axis from WCS keywords
    # ------------------------------------------------------------------
    crval3 = float(header.get('CRVAL3', 0.0))   # m/s
    crpix3 = float(header.get('CRPIX3', 1.0))   # 1-indexed reference pixel
    cdelt3 = float(header.get('CDELT3', 1.0))   # m/s per channel

    channels_kms = (crval3 + (np.arange(nvel) - (crpix3 - 1)) * cdelt3) / 1e3

    # ------------------------------------------------------------------
    # 3. Select channels within velocity range
    # ------------------------------------------------------------------
    if velocity_range is not None:
        v_min, v_max = float(velocity_range[0]), float(velocity_range[1])
        chan_mask = (channels_kms >= v_min) & (channels_kms <= v_max)
        if not np.any(chan_mask):
            raise ValueError(
                f"No channels found in velocity range [{v_min}, {v_max}] km/s. "
                f"Cube covers {channels_kms.min():.1f} – {channels_kms.max():.1f} km/s."
            )
        selected_channels = channels_kms[chan_mask]
        print(f"  Integrating {chan_mask.sum()} channels: "
              f"{selected_channels.min():.1f} – {selected_channels.max():.1f} km/s")
    else:
        chan_mask = np.ones(nvel, dtype=bool)
        print(f"  Integrating all {nvel} channels: "
              f"{channels_kms.min():.1f} – {channels_kms.max():.1f} km/s")

    # ------------------------------------------------------------------
    # 4. Collapse
    # ------------------------------------------------------------------
    deltav_kms = abs(cdelt3) / 1e3

    # --- Optional velocity-field SHUFFLE (align lines, then integrate narrow) ---
    # Build a smooth per-pixel line-velocity field from the cube itself (moment-1
    # on high-S/N pixels, smoothed+filled), shift every spectrum so its line lands
    # at a common ref velocity, then reset the integration to a NARROW window
    # around that ref velocity. This removes galaxy rotation from the *velocity
    # axis* so a tight moment-0 window catches the line everywhere → far less
    # integrated noise → faint extended/spiral structure survives. The field is
    # smooth (not a per-pixel self-peak), so faint pixels get the right shift and
    # need no snr-threshold. Meant for mode=moment0.
    if shuffle:
        if velocity_range is None:
            raise ValueError("--shuffle requires --velocity-range (the line window "
                             "used to measure the velocity field).")
        line_mask0 = chan_mask.copy()
        v0_raw, v0_field, snr_field = _velocity_field_from_cube(
            cube, channels_kms, line_mask0, deltav_kms, shuffle_snr, shuffle_field_smooth)
        n_good = int(np.sum(np.isfinite(v0_raw)))
        if n_good == 0:
            raise ValueError(f"--shuffle: no pixels reached S/N ≥ {shuffle_snr} in "
                             f"[{v_min:.0f}, {v_max:.0f}] km/s — cannot build a velocity "
                             f"field. Lower --shuffle-snr or widen --velocity-range.")
        ref_v = (float(shuffle_ref_velocity) if shuffle_ref_velocity is not None
                 else float(np.nanmedian(v0_raw)))
        print(f"  Shuffle: velocity field from {n_good} pixels (S/N ≥ {shuffle_snr}), "
              f"v0 range {np.nanmin(v0_raw):.1f}–{np.nanmax(v0_raw):.1f} km/s, "
              f"field-smooth {shuffle_field_smooth:.1f} px")
        print(f"  Shuffle: aligning lines to ref velocity {ref_v:.1f} km/s")

        # Optionally write the velocity field for inspection (bad field = main risk).
        _vf_out = shuffle_field_output
        if _vf_out is None:
            _base = (plot_output or fits_output or cube_fits)
            _vf_out = str(Path(_base).with_suffix('')) + '_vfield.fits'
        try:
            _vf_hdr = WCS(header).dropaxis(2).to_header()
            _vf_hdr['BUNIT'] = 'km/s'
            _vf_hdr['OBJECT'] = header.get('OBJECT', '').strip()
            _vf_hdr['REFV'] = (ref_v, '[km/s] shuffle reference velocity')
            fits.writeto(_vf_out, v0_field.astype(np.float32), _vf_hdr, overwrite=True)
            print(f"  Shuffle: velocity field written → {_vf_out}")
        except Exception as _e:
            print(f"  Shuffle: could not write velocity field ({_e})")

        cube = _shuffle_cube(cube, channels_kms, v0_field, ref_v)

        # Reset integration to a narrow window around the common ref velocity.
        v_min = ref_v - shuffle_window_kms
        v_max = ref_v + shuffle_window_kms
        velocity_range = (v_min, v_max)
        chan_mask = (channels_kms >= v_min) & (channels_kms <= v_max)
        print(f"  Shuffle: integrating aligned window {v_min:.1f} – {v_max:.1f} km/s "
              f"({int(chan_mask.sum())} channels)")

    # --- Optional spectral (velocity-axis) smoothing used only to LOCATE the
    #     per-pixel peak, never to set the reported intensity. A narrow residual
    #     spike (e.g. an un-refilled telluric wing) can win a raw np.nanmax, but a
    #     Gaussian-smoothed spectrum discriminates by width: the broad real line
    #     survives while a 1–2 channel spike is crushed. We take argmax on the
    #     smoothed cube, then read the RAW cube at that channel so the map keeps
    #     true Kelvin amplitudes. astropy (not scipy) convolution to match the
    #     spatial-smoothing path and the codebase's cygrid/astropy convention.
    cube_peakfind = cube
    if peak_smooth_channels is not None and peak_smooth_channels > 0 \
            and mode in ('peak-intensity', 'peak-range-int'):
        cube_peakfind = _spectral_smooth_nanaware(cube, float(peak_smooth_channels))
        print(f"  Peak location via spectral smoothing: sigma = "
              f"{peak_smooth_channels:.1f} channels ({peak_smooth_channels*deltav_kms:.1f} km/s); "
              f"reported intensities are raw (unsmoothed)")

    if mode == 'peak-intensity':
        # Locate the peak channel on the (optionally smoothed) cube, report raw value there.
        peak_idx = np.argmax(np.where(np.isnan(cube_peakfind), -np.inf, cube_peakfind), axis=0)
        _ry = np.arange(ny)[:, None]
        _rx = np.arange(nx)[None, :]
        collapsed = cube[peak_idx, _ry, _rx].astype(np.float32)
        _map_units = 'K'
        _map_label = 'Peak intensity (K)'
        print(f"  Peak-intensity: full spectrum ({nvel} channels, "
              f"{channels_kms.min():.1f} – {channels_kms.max():.1f} km/s)")
    elif mode == 'peak-range-int':
        sub = cube[chan_mask, :, :]          # (nchan_sel, ny, nx)
        sub_pf = cube_peakfind[chan_mask, :, :]
        nchan_sel = sub.shape[0]
        # Per-pixel peak channel: replace NaN with -inf so argmax works on all-NaN pixels
        peak_idx = np.argmax(np.where(np.isnan(sub_pf), -np.inf, sub_pf), axis=0)  # (ny, nx)
        n = peak_range_channels
        print(f"  Peak-range: ±{n} channels around per-pixel peak "
              f"(window up to {2*n+1} ch)")
        acc = np.zeros((ny, nx), dtype=np.float64)
        row_idx = np.arange(ny)[:, None]
        col_idx = np.arange(nx)[None, :]
        for k in range(-n, n + 1):
            chan_k = peak_idx + k
            valid = (chan_k >= 0) & (chan_k < nchan_sel)
            vals = sub[np.clip(chan_k, 0, nchan_sel - 1), row_idx, col_idx]
            acc += np.where(valid & ~np.isnan(vals), vals, 0.0)
        collapsed = (acc * deltav_kms).astype(np.float32)
        _map_units = 'K km/s'
        _map_label = f'Peak ±{n}ch (K km/s)'
    elif collapse_weights:  # weight cube: plain sum of Σ K·w, no channel-width factor
        collapsed = np.nansum(cube[chan_mask, :, :], axis=0)
        _map_units = 'weight'
        _map_label = 'Integrated weight (Σ K·w)'
    else:  # moment0
        collapsed = np.nansum(cube[chan_mask, :, :], axis=0) * deltav_kms
        _map_units = 'K km/s'
        _map_label = 'K km/s'

    # ------------------------------------------------------------------
    # 4b. Build display mask (NaN = not shown in plot; FITS always unmasked)
    # ------------------------------------------------------------------
    display_mask = np.zeros((ny, nx), dtype=bool)  # True = masked out
    # Pixels rejected by --snr-threshold / --drop-rms are shown as the zero
    # colour rather than blanked (NaN); tracked here so the interactive
    # velocity-recompute path (which re-NaNs display_mask) can re-zero them.
    zero_display_mask = np.zeros((ny, nx), dtype=bool)

    # --- Coverage-based edge masking ---
    # Try to load the COVERAGE extension saved by create_spectral_datacube.
    coverage_map = None
    if coverage_threshold > 0:
        try:
            with fits.open(cube_fits) as _hdul:
                if 'COVERAGE' in [h.name for h in _hdul]:
                    coverage_map = _hdul['COVERAGE'].data.astype(float)
        except Exception:
            pass

        if coverage_map is not None:
            peak_cov = coverage_map.max()
            if peak_cov > 0:
                # Simply ignore (blank to no-data) pixels below the coverage
                # threshold — no zero-colour fill.
                edge_mask = coverage_map < coverage_threshold * peak_cov
                display_mask |= edge_mask
                print(f"  Coverage masking: {int(edge_mask.sum())} pixels below "
                      f"{coverage_threshold*100:.0f}% of peak coverage blanked")
        else:
            print("  Note: no COVERAGE extension found in cube — "
                  "re-run create_datacube to enable coverage masking")

    # --- Circular aperture mask (blank everything outside the circle) ---
    if mask_radius_arcmin is not None:
        pixel_scale_deg = abs(float(header.get('CDELT2', 1.0)))
        radius_pix = (mask_radius_arcmin / 60.0) / pixel_scale_deg
        if mask_ra is not None and mask_dec is not None:
            from astropy.wcs import WCS as _WCS
            from astropy.coordinates import SkyCoord
            import astropy.units as _u
            _wcs2d = WCS(header).dropaxis(2)
            _sky = SkyCoord(ra=mask_ra * _u.deg, dec=mask_dec * _u.deg)
            cx_m, cy_m = _wcs2d.world_to_pixel(_sky)
            cx_m, cy_m = float(cx_m), float(cy_m)
        else:
            cx_m, cy_m = nx / 2.0, ny / 2.0
        yy, xx = np.mgrid[0:ny, 0:nx]
        outside_circle = (xx - cx_m) ** 2 + (yy - cy_m) ** 2 > radius_pix ** 2
        display_mask |= outside_circle
        print(f"  Circular mask: blanking outside r={mask_radius_arcmin:.1f}′ "
              f"(centre pixel {cx_m:.1f}, {cy_m:.1f})")

    # --- Edge trim: erode the irregular data footprint inward ---
    # Unlike a rectangular crop, this follows the (irregular) telescope-coverage
    # boundary and peels a border of uniform thickness off it, so the map keeps
    # its general shape. trim_edges is a string: a bare number or "N%" peels a
    # border of N% of the smaller map axis; "Narcsec"/"Narcmin" (or ", ') peels
    # a fixed angular border.
    if trim_edges is not None:
        s = str(trim_edges).strip().lower()
        pixel_scale_arcsec = abs(float(header.get('CDELT2', 1.0))) * 3600.0
        try:
            if s.endswith('arcmin') or s.endswith('am') or s.endswith("'"):
                val = float(s.rstrip("'").replace('arcmin', '').replace('am', ''))
                border_px = int(round(val * 60.0 / pixel_scale_arcsec))
                desc = f"{val:g}′ ({border_px}px)"
            elif s.endswith('arcsec') or s.endswith('as') or s.endswith('"'):
                val = float(s.rstrip('"').replace('arcsec', '').replace('as', ''))
                border_px = int(round(val / pixel_scale_arcsec))
                desc = f"{val:g}″ ({border_px}px)"
            else:  # percent (bare number or trailing %) of the smaller map axis
                pct = float(s.rstrip('%'))
                border_px = int(round(min(nx, ny) * pct / 100.0))
                desc = f"{pct:g}% of map ({border_px}px)"
        except ValueError:
            border_px = 0
            desc = None
            print(f"  Warning: could not parse --trim-edges '{trim_edges}'; ignoring")

        if desc is not None and border_px > 0:
            # Base footprint: pixels that actually hold integrated data. Coverage,
            # if present, refines the true telescope footprint.
            footprint = ~np.all(np.isnan(cube[chan_mask, :, :]), axis=0)
            if coverage_map is not None:
                footprint &= (coverage_map > 0)

            # Fill interior holes so we trim ONLY the outer boundary, not rings
            # around interior masked pixels. "Exterior" = background reachable
            # from the array border via a 4-connected flood fill.
            from collections import deque
            bg = ~footprint
            exterior = np.zeros_like(footprint)
            dq = deque()
            for x in range(nx):
                for yb in (0, ny - 1):
                    if bg[yb, x] and not exterior[yb, x]:
                        exterior[yb, x] = True; dq.append((yb, x))
            for y in range(ny):
                for xb in (0, nx - 1):
                    if bg[y, xb] and not exterior[y, xb]:
                        exterior[y, xb] = True; dq.append((y, xb))
            while dq:
                y, x = dq.popleft()
                for dy, dx in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                    yy2, xx2 = y + dy, x + dx
                    if 0 <= yy2 < ny and 0 <= xx2 < nx and bg[yy2, xx2] \
                            and not exterior[yy2, xx2]:
                        exterior[yy2, xx2] = True; dq.append((yy2, xx2))
            filled = ~exterior

            # 8-connected erosion by border_px: a pixel survives only if all of
            # its neighbours are inside the filled footprint. Pixels outside the
            # array count as background, so a footprint touching the array edge
            # is trimmed there too.
            eroded = filled
            for _ in range(border_px):
                e = eroded.copy()
                e[:-1, :] &= eroded[1:, :]
                e[1:, :] &= eroded[:-1, :]
                e[:, :-1] &= eroded[:, 1:]
                e[:, 1:] &= eroded[:, :-1]
                e[:-1, :-1] &= eroded[1:, 1:]
                e[:-1, 1:] &= eroded[1:, :-1]
                e[1:, :-1] &= eroded[:-1, 1:]
                e[1:, 1:] &= eroded[:-1, :-1]
                e[0, :] = False; e[-1, :] = False
                e[:, 0] = False; e[:, -1] = False
                eroded = e

            trim_ring = footprint & ~eroded
            if eroded.sum() == 0:
                print(f"  Warning: --trim-edges {trim_edges} would blank the whole "
                      f"footprint ({border_px}px erosion); ignoring")
            else:
                display_mask |= trim_ring
                print(f"  Edge trim: peeled {desc} off the coverage footprint "
                      f"following its shape ({int(trim_ring.sum())} pixels)")

    # Pixels where every channel in the integration range is NaN have no data.
    # np.nansum returns 0 for those, so we mask them explicitly here so they
    # show as white in the plot rather than as the colormap colour for zero.
    all_nan_pixels = np.all(np.isnan(cube[chan_mask, :, :]), axis=0)
    display_mask |= all_nan_pixels

    collapsed_display = collapsed.copy()
    collapsed_display[display_mask] = np.nan
    n_masked = int(display_mask.sum())
    if n_masked > 0:
        print(f"  Total masked pixels in display: {n_masked}")

    print(f"  Collapsed map ({mode}): {ny} × {nx} pixels, "
          f"range [{np.nanmin(collapsed_display):.2f}, {np.nanmax(collapsed_display):.2f}] {_map_units}")

    # ------------------------------------------------------------------
    # 4c. Noise map: per-pixel RMS from line-free channels
    # ------------------------------------------------------------------
    noise_map = None
    noise_chan_mask = ~chan_mask
    if collapse_weights:
        pass  # a weight cube has no line-free "noise" to measure
    elif noise_chan_mask.sum() >= 2:
        noise_map = np.nanstd(cube[noise_chan_mask, :, :], axis=0).astype(np.float32)
        noise_map[display_mask] = np.nan
        avg_noise = float(np.nanmean(noise_map))
        print(f"  Noise map: {noise_chan_mask.sum()} line-free channels, "
              f"mean noise = {avg_noise:.4f} K")
    elif velocity_range is None:
        print("  Noise map: no --velocity-range given — cannot identify line-free channels")
    else:
        print("  Noise map: insufficient line-free channels for noise estimation")

    # ------------------------------------------------------------------
    # 4c2. SNR threshold: mask pixels below signal-to-noise cutoff
    # ------------------------------------------------------------------
    # SNR = collapsed [K km/s] / (noise [K] × sqrt(n_signal_chan) × deltav [km/s])
    if snr_threshold is not None and noise_map is not None:
        n_signal_chan = int(chan_mask.sum())
        if mode == 'peak-intensity':
            # peak value [K] / per-channel noise [K]
            snr_map = collapsed / noise_map
        else:
            # moment-0 [K km/s] / (noise [K] × sqrt(n) × deltav [km/s])
            snr_map = collapsed / (noise_map * np.sqrt(n_signal_chan) * deltav_kms)
        snr_mask = np.isfinite(snr_map) & (snr_map < snr_threshold)
        display_mask |= snr_mask
        # Rejected pixels are shown as the zero colour (like suppress_negative),
        # NOT blanked to no-data (NaN), while still kept in display_mask so they
        # stay out of the mean-spectrum panel.
        zero_display_mask |= snr_mask
        collapsed_display[snr_mask] = 0.0
        print(f"  SNR threshold: {snr_threshold:.1f} — masked {snr_mask.sum()} "
              f"additional pixels (shown as 0)")
    elif snr_threshold is not None and noise_map is None:
        print("  Warning: --snr-threshold ignored — noise map unavailable "
              "(need --velocity-range to identify line-free channels)")

    # ------------------------------------------------------------------
    # 4c2b. Drop-RMS: blank pixels whose per-pixel spectral RMS in a chosen
    #       velocity window exceeds a limit. RMS is measured directly on the
    #       gridded/convolved cube (sqrt of the mean square over the window's
    #       channels), so it flags noisy spaxels an SNR cut can miss.
    # ------------------------------------------------------------------
    if drop_rms is not None:
        if drop_window is not None:
            dw_min, dw_max = float(drop_window[0]), float(drop_window[1])
            drop_chan_mask = (channels_kms >= dw_min) & (channels_kms <= dw_max)
            drop_win_desc = f"{dw_min:.1f}–{dw_max:.1f} km/s"
        else:
            # Fall back to the line-free channels (outside --velocity-range).
            drop_chan_mask = noise_chan_mask
            drop_win_desc = "line-free channels (outside --velocity-range)"
        if drop_chan_mask.sum() < 2:
            print("  Warning: --drop-rms ignored — need >=2 channels in the drop "
                  "window (give --drop-window, or --velocity-range for the default)")
        else:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", category=RuntimeWarning)
                rms_map = np.sqrt(np.nanmean(
                    cube[drop_chan_mask, :, :] ** 2, axis=0))
            # Only drop pixels not already masked, with a finite RMS above LIMIT.
            drop_mask = np.isfinite(rms_map) & (rms_map > drop_rms) & ~display_mask
            display_mask |= drop_mask
            # Shown as the zero colour (like --snr-threshold / suppress_negative),
            # NOT blanked to no-data (NaN); still kept in display_mask so they
            # stay out of the mean-spectrum panel.
            zero_display_mask |= drop_mask
            collapsed_display[drop_mask] = 0.0
            print(f"  Drop-RMS: RMS over {drop_win_desc} "
                  f"({int(drop_chan_mask.sum())} channels) > {drop_rms:g} K "
                  f"-> masked {int(drop_mask.sum())} additional pixels (shown as 0)")

    # ------------------------------------------------------------------
    # 4c3. Suppress negative values: set to 0 rather than NaN so they
    #      show as the zero colour rather than white (no-data).
    # ------------------------------------------------------------------
    if suppress_negative:
        neg_mask = np.isfinite(collapsed_display) & (collapsed_display < 0)
        collapsed_display[neg_mask] = 0.0
        if neg_mask.sum() > 0:
            print(f"  Suppress negative: {neg_mask.sum()} pixels clipped to 0")

    # ------------------------------------------------------------------
    # 4c4. Suppress high values: clip pixels above threshold to threshold.
    # ------------------------------------------------------------------
    if suppress_high is not None:
        high_mask = np.isfinite(collapsed_display) & (collapsed_display > suppress_high)
        collapsed_display[high_mask] = suppress_high
        if high_mask.sum() > 0:
            print(f"  Suppress high: {high_mask.sum()} pixels clipped to {suppress_high}")

    # ------------------------------------------------------------------
    # 4d. Optional Gaussian smoothing of the display map (NaN-aware)
    # ------------------------------------------------------------------
    if smooth_sigma is not None and smooth_sigma > 0:
        from astropy.convolution import Gaussian2DKernel, convolve as _convolve
        _kernel = Gaussian2DKernel(x_stddev=smooth_sigma)
        collapsed_display = _convolve(collapsed_display, _kernel,
                                      boundary='fill', fill_value=np.nan,
                                      nan_treatment='interpolate',
                                      preserve_nan=True).astype(np.float32)
        print(f"  Gaussian smoothing applied: sigma = {smooth_sigma:.1f} pixels")

    # ------------------------------------------------------------------
    # 4e. Restrict the rendered area to a polygon (interactive or explicit)
    # ------------------------------------------------------------------
    # Interactive selection draws the current display map and lets the user
    # enclose the region of interest; explicit `polygon` vertices reuse a
    # previous selection headlessly. Everything outside the polygon is blanked
    # in the display, the noise map, and any spectrum extraction (same as the
    # circular mask). The FITS output keeps the full raw map.
    if select_polygon and polygon is None:
        polygon = _interactive_polygon(
            collapsed_display, colormap=colormap,
            title=(f"{header.get('OBJECT', '').strip()}  "
                   "— Middle-click: add vertex | Left-click: close | "
                   "Right-click: undo"))
        if polygon is None:
            print("  Polygon selection: fewer than 3 vertices placed — no polygon applied")
        else:
            print("  Polygon selection: "
                  + " ".join(f"{x:.1f},{y:.1f}" for x, y in polygon))

    if polygon is not None and len(polygon) >= 3:
        outside_poly = _polygon_mask(polygon, ny, nx)
        display_mask |= outside_poly
        collapsed_display[outside_poly] = np.nan
        if noise_map is not None:
            noise_map[outside_poly] = np.nan
        print(f"  Polygon mask: {int(outside_poly.sum())} pixels blanked "
              f"outside the {len(polygon)}-vertex polygon")

    # ------------------------------------------------------------------
    # 5. Build 2D WCS header (drop spectral axis)
    # ------------------------------------------------------------------
    wcs3d = WCS(header)
    wcs2d = wcs3d.dropaxis(2)   # drop velocity (3rd WCS axis, 0-indexed)
    header2d = wcs2d.to_header()
    header2d['NAXIS']  = 2
    header2d['NAXIS1'] = nx
    header2d['NAXIS2'] = ny
    header2d['BUNIT']  = _map_units
    header2d['OBJECT'] = header.get('OBJECT', '').strip()
    for kw in ('BMAJ', 'BMIN', 'BPA', 'TELESCOP', 'RESTFRQ'):
        if kw in header:
            header2d[kw] = header[kw]
    if velocity_range is not None:
        header2d['VMIN'] = (v_min, '[km/s] integration start')
        header2d['VMAX'] = (v_max, '[km/s] integration end')

    # ------------------------------------------------------------------
    # 6. Save FITS if requested
    # ------------------------------------------------------------------
    if fits_output:
        fits.writeto(fits_output,
                     collapsed.astype(np.float32),
                     header2d,
                     overwrite=overwrite)
        print(f"  FITS saved: {fits_output}")

    if noise_map_output:
        if noise_map is not None:
            noise_header = header2d.copy()
            noise_header['BUNIT'] = 'K'
            noise_header['COMMENT'] = 'Per-pixel RMS noise from line-free channels'
            fits.writeto(noise_map_output,
                         noise_map.astype(np.float32),
                         noise_header,
                         overwrite=overwrite)
            print(f"  Noise map FITS saved: {noise_map_output}")
        else:
            print("  Warning: noise map not computed — --output-noise-map skipped")

    # ------------------------------------------------------------------
    # 7. Plot
    # ------------------------------------------------------------------
    pixel_scale_deg = abs(float(header2d.get('CDELT2', 1.0)))
    obj_name = header.get('OBJECT', '').strip()
    vrange_str = (f"  [{v_min:.0f} – {v_max:.0f} km/s]"
                  if velocity_range is not None else "")
    title_base = f"{obj_name}  integrated intensity{vrange_str}"

    from matplotlib.patches import Rectangle, Circle

    # ------------------------------------------------------------------
    # RA/Dec → pixel conversion helper (uses the 2D WCS)
    # ------------------------------------------------------------------
    def _radec_to_pix(ra_deg, dec_deg):
        """Return (x_pix, y_pix) for a sky position using the 2D WCS."""
        from astropy.wcs import WCS as _WCS
        from astropy.coordinates import SkyCoord
        import astropy.units as u
        sky = SkyCoord(ra=ra_deg * u.deg, dec=dec_deg * u.deg)
        x_pix, y_pix = wcs2d.world_to_pixel(sky)
        return float(x_pix), float(y_pix)

    # ------------------------------------------------------------------
    # Zoom region (pixel bounds + mean spectrum)
    # ------------------------------------------------------------------
    if zoom_size_arcmin is not None:
        zoom_pix = int(round((zoom_size_arcmin / 60.0) / pixel_scale_deg))
        if zoom_x is not None and zoom_y is not None:
            cx, cy = int(round(zoom_x)), int(round(zoom_y))
        elif zoom_ra is not None and zoom_dec is not None:
            _zx, _zy = _radec_to_pix(zoom_ra, zoom_dec)
            cx, cy = int(round(_zx)), int(round(_zy))
            print(f"  Zoom centre RA={zoom_ra:.5f}° Dec={zoom_dec:.5f}° "
                  f"→ pixel ({cx}, {cy})")
        else:
            cy, cx = ny // 2, nx // 2
        y0 = max(0, cy - zoom_pix // 2)
        y1 = min(ny,  cy + zoom_pix // 2)
        x0 = max(0, cx - zoom_pix // 2)
        x1 = min(nx,  cx + zoom_pix // 2)
        zoomed    = collapsed_display[y0:y1, x0:x1]
        zoom_cube_slice = cube[:, y0:y1, x0:x1]
        if align_peaks:
            zoom_spec = _peak_aligned_mean_spectrum(
                zoom_cube_slice.reshape(zoom_cube_slice.shape[0], -1).T,
                channels_kms, search_range=velocity_range)
        else:
            zoom_spec = np.nanmean(zoom_cube_slice, axis=(1, 2))
    else:
        y0 = y1 = x0 = x1 = None
        zoomed = zoom_spec = None

    # ------------------------------------------------------------------
    # Circular extraction region (pixel centre + mask + mean spectrum)
    # ------------------------------------------------------------------
    if region_radius_arcmin is not None:
        radius_pix = region_radius_arcmin / 60.0 / pixel_scale_deg
        if region_x is not None or region_y is not None:
            cx_reg = float(region_x) if region_x is not None else nx / 2.0
            cy_reg = float(region_y) if region_y is not None else ny / 2.0
        elif region_ra is not None and region_dec is not None:
            cx_reg, cy_reg = _radec_to_pix(region_ra, region_dec)
            print(f"  Region centre RA={region_ra:.5f}° Dec={region_dec:.5f}° "
                  f"→ pixel ({cx_reg:.1f}, {cy_reg:.1f})")
        else:
            cx_reg, cy_reg = nx / 2.0, ny / 2.0
        reg_label = (f"r = {region_radius_arcmin:.1f}′  "
                     f"(x={cx_reg:.1f}, y={cy_reg:.1f})")
        yy, xx = np.mgrid[0:ny, 0:nx]
        region_mask = (xx - cx_reg) ** 2 + (yy - cy_reg) ** 2 <= radius_pix ** 2
        n_reg_pix = region_mask.sum()
        if n_reg_pix == 0:
            raise ValueError("Region aperture contains no map pixels — "
                             "check pixel coordinates and radius.")
        region_cube_slice = cube[:, region_mask]  # shape (nvel, n_reg_pix)
        if align_peaks:
            region_spec = _peak_aligned_mean_spectrum(
                region_cube_slice.T, channels_kms, search_range=velocity_range)
        else:
            region_spec = np.nanmean(region_cube_slice, axis=1)
        print(f"  Region aperture: {n_reg_pix} pixels, "
              f"centre ({cx_reg:.1f}, {cy_reg:.1f}), r = {radius_pix:.1f} px")
    else:
        region_mask = region_spec = cx_reg = cy_reg = radius_pix = None
        reg_label = None

    # ------------------------------------------------------------------
    # If neither zoom nor region: use full-map mean spectrum
    # ------------------------------------------------------------------
    # The spectrum is averaged only over the pixels actually rendered
    # (~display_mask), so a polygon (--select-polygon/--polygon) or circular
    # aperture (--mask-radius), and any coverage/trim/SNR cut, restrict the
    # spectrum to the same area shown in the map.
    full_spec_label = "Mean spectrum — full map"
    if zoom_size_arcmin is None and region_radius_arcmin is None:
        keep_pix = ~display_mask                         # (ny, nx) rendered pixels
        n_keep = int(keep_pix.sum())
        restricted = bool(
            polygon is not None or mask_radius_arcmin is not None
            or trim_edges is not None or snr_threshold is not None)
        if n_keep == 0:
            full_spec = np.full(cube.shape[0], np.nan)
        else:
            sel = cube.reshape(cube.shape[0], -1)[:, keep_pix.ravel()]
            if align_peaks:
                full_spec = _peak_aligned_mean_spectrum(
                    sel.T, channels_kms, search_range=velocity_range)
            else:
                full_spec = np.nanmean(sel, axis=1)
        if restricted:
            full_spec_label = f"Mean spectrum — selected region ({n_keep} px)"
    else:
        full_spec = None

    # ------------------------------------------------------------------
    # WCS objects for map axes (built once, reused per panel)
    # ------------------------------------------------------------------
    def _make_wcs_zoom(zx0, zy0):
        # Build a WCS for the zoom panel whose reference pixel is shifted to
        # account for the zoom box's lower-left corner (zx0, zy0) in the
        # full-map pixel grid. Must be recomputed whenever the zoom box is
        # redefined interactively, otherwise the zoom panel's coordinate
        # ticks/labels are computed for the wrong pixel offset and vanish.
        h_zoom = header2d.copy()
        h_zoom['CRPIX1'] = float(header2d.get('CRPIX1', nx / 2 + 1)) - zx0
        h_zoom['CRPIX2'] = float(header2d.get('CRPIX2', ny / 2 + 1)) - zy0
        return WCS(h_zoom)

    if use_wcs:
        wcs2d_obj = WCS(header2d)
        if zoom_size_arcmin is not None:
            wcs2d_zoom = _make_wcs_zoom(x0, y0)
        else:
            wcs2d_zoom = None
    else:
        wcs2d_obj = wcs2d_zoom = None

    # ------------------------------------------------------------------
    # 2×2 figure layout
    #   [0,0] full map          (always)
    #   [0,1] zoom map          (if zoom,   else hidden)
    #   [1,0] zoom/full spectrum(always)
    #   [1,1] region spectrum   (if region, else hidden)
    # ------------------------------------------------------------------
    has_zoom   = zoom_size_arcmin is not None
    has_region = region_radius_arcmin is not None

    fig = plt.figure(figsize=(16, 11))

    # Give the map row more height than the spectra row. The maps are square
    # (equal sky aspect) and share their cell with a colorbar, so they need a
    # wide, tall top row to grow; the spectra are naturally wide-and-short, so
    # the smaller bottom row costs them little. The 16"-wide figure + 1.75:1
    # height split renders maps ≈5.7" (vs ≈4.5" at the old 14"/1:1 layout) and
    # keeps the spectra ≈3.4" tall. tight_layout() (below) preserves the ratio.
    gs = fig.add_gridspec(2, 2, height_ratios=[1.75, 1])

    def _add_map_subplot(cell, wcs_proj=None):
        if wcs_proj is not None:
            return fig.add_subplot(cell, projection=wcs_proj)
        return fig.add_subplot(cell)

    ax_full      = _add_map_subplot(gs[0, 0], wcs2d_obj)
    ax_zoom      = _add_map_subplot(gs[0, 1], wcs2d_zoom) if has_zoom else fig.add_subplot(gs[0, 1])
    ax_spec_zoom = fig.add_subplot(gs[1, 0])
    ax_spec_reg  = fig.add_subplot(gs[1, 1])

    # fig.colorbar() shrinks the host axes to make room for the colorbar.
    # When a map panel is redrawn (cla() + _plot_map), repeating that on an
    # already-shrunk axes would shrink it further each time, eventually
    # pushing the colorbar (and its label) off the figure. Restore each
    # axes to its original position before redrawing.
    _ax_full_pos = ax_full.get_position().frozen()
    _ax_zoom_pos = ax_zoom.get_position().frozen()

    if not has_zoom:
        ax_zoom.set_visible(False)
    if not has_region:
        ax_spec_reg.set_visible(False)

    # Beam size in pixels for hex grid spacing (beam/2)
    _bmaj_deg = float(header.get('BMAJ', beamsize_deg_for_hex := pixel_scale_deg * 3.0))
    _hex_spacing_pix = (_bmaj_deg / 2.0) / pixel_scale_deg  # beam/2 in pixels

    def _ax_labels(ax, wcs_proj):
        if wcs_proj is not None:
            ax.coords[0].set_axislabel('RA')
            ax.coords[1].set_axislabel('Dec')
            ax.coords[0].set_major_formatter('hh:mm:ss')
            ax.coords[1].set_major_formatter('dd:mm:ss')
        else:
            ax.set_xlabel('RA pixel')
            ax.set_ylabel('Dec pixel')

    def _vminmax(data):
        lo, hi = percentile_clip if percentile_clip is not None else (0, 100)
        # suppress_negative forces vmin=0 only for moment0: peak-intensity values
        # are always positive, so vmin=0 would compress the noise floor to black.
        if suppress_negative and mode == 'moment0':
            vmin_p = 0.0
        else:
            vmin_p = np.nanpercentile(data, lo)
        vmax_p = np.nanpercentile(data, hi)
        return vmin_p, vmax_p

    def _make_norm(vmin_p, vmax_p, data=None):
        from matplotlib import colors as mcolors
        if stretch == 'sqrt':
            # shift so vmin maps to 0, then apply power 0.5
            return mcolors.PowerNorm(gamma=0.5, vmin=vmin_p, vmax=vmax_p)
        elif stretch == 'power':
            # user-tunable power law: gamma<1 lifts faint detail, gamma>1
            # emphasises bright peaks; gamma=1 is linear, gamma=0.5 is sqrt.
            return mcolors.PowerNorm(gamma=gamma, vmin=vmin_p, vmax=vmax_p)
        elif stretch == 'log':
            safe_vmin = max(vmin_p, 1e-6 * vmax_p) if vmax_p > 0 else 1e-6
            return mcolors.LogNorm(vmin=safe_vmin, vmax=vmax_p)
        elif stretch == 'symlog':
            # linear within +-linthresh of zero, logarithmic beyond — same
            # colourbar shape as asinh but with an explicit linear threshold,
            # and it compresses negative noise dips symmetrically. linthresh is
            # set to ~10% of the colour scale (the noise-level regime).
            linthresh = max(abs(vmax_p), abs(vmin_p)) * 0.1
            if linthresh <= 0:
                linthresh = 1e-6
            return mcolors.SymLogNorm(linthresh=linthresh, vmin=vmin_p, vmax=vmax_p)
        elif stretch == 'asinh':
            # AsinhNorm available from matplotlib 3.2+; fall back to sqrt if missing
            if hasattr(mcolors, 'AsinhNorm'):
                linear_width = (vmax_p - vmin_p) * 0.1
                return mcolors.AsinhNorm(linear_width=linear_width,
                                         vmin=vmin_p, vmax=vmax_p)
            else:
                import warnings
                warnings.warn("AsinhNorm requires matplotlib >= 3.2, falling back to sqrt")
                return mcolors.PowerNorm(gamma=0.5, vmin=vmin_p, vmax=vmax_p)
        elif stretch == 'histeq':
            # Histogram equalisation: maps the value distribution to a flat
            # histogram, maximising displayed contrast at the expense of a
            # non-uniform (non-quantitative) colourbar. Built from the finite
            # data within the clip range so it respects --percentile-clip.
            try:
                from astropy.visualization import ImageNormalize, HistEqStretch
                if data is None:
                    raise ValueError("histeq needs the data array")
                fin = np.asarray(data)[np.isfinite(data)]
                sel = fin[(fin >= vmin_p) & (fin <= vmax_p)] if vmax_p > vmin_p else fin
                if sel.size < 2:
                    sel = fin
                if sel.size < 2:
                    return mcolors.Normalize(vmin=vmin_p, vmax=vmax_p)
                return ImageNormalize(vmin=vmin_p, vmax=vmax_p,
                                      stretch=HistEqStretch(sel), clip=False)
            except Exception as _e:
                import warnings
                warnings.warn(f"histeq stretch unavailable ({_e}); falling back to linear")
                return mcolors.Normalize(vmin=vmin_p, vmax=vmax_p)
        else:  # linear
            return mcolors.Normalize(vmin=vmin_p, vmax=vmax_p)

    def _imshow_map(ax, data, title, wcs_proj=None, cmap=colormap):
        vmin_p, vmax_p = _vminmax(data)
        norm = _make_norm(vmin_p, vmax_p, data)
        cmap_obj = plt.colormaps[cmap].copy() if isinstance(cmap, str) else cmap.copy()
        cmap_obj.set_bad('white')
        im = ax.imshow(data, origin='lower', cmap=cmap_obj,
                       norm=norm, interpolation=interpolation)
        ax.set_title(title, fontsize=10)
        _ax_labels(ax, wcs_proj)
        return fig.colorbar(im, ax=ax, label=_map_label, fraction=0.046, pad=0.04)

    _hex_scatters = []  # (scatter, spacing_pix, ax) — sizes updated after tight_layout

    def _hexplot_map(ax, data, title, wcs_proj=None, cmap=colormap):
        from scipy.ndimage import map_coordinates
        _ny, _nx = data.shape
        spacing = _hex_spacing_pix
        row_spacing = spacing * np.sqrt(3) / 2.0
        xs, ys = [], []
        row = 0
        y = 0.0
        while y <= _ny - 1:
            offset = (spacing / 2.0) if row % 2 else 0.0
            x = offset
            while x <= _nx - 1:
                xs.append(x)
                ys.append(y)
                x += spacing
            y += row_spacing
            row += 1
        xs = np.array(xs, dtype=float)
        ys = np.array(ys, dtype=float)
        vals = map_coordinates(data, [ys, xs], order=1, prefilter=False, cval=np.nan)
        finite = np.isfinite(vals)
        xs, ys, vals = xs[finite], ys[finite], vals[finite]
        vmin_p, vmax_p = _vminmax(vals)
        norm = _make_norm(vmin_p, vmax_p, vals)
        ax.set_xlim(-0.5, _nx - 0.5)
        ax.set_ylim(-0.5, _ny - 0.5)
        tf = ax.get_transform('pixel') if wcs_proj is not None else ax.transData
        sc = ax.scatter(xs, ys, c=vals, marker='h', s=1, cmap=cmap,
                        norm=norm, linewidths=0, transform=tf)
        _hex_scatters.append((sc, spacing, ax))
        ax.set_title(title, fontsize=10)
        _ax_labels(ax, wcs_proj)
        return fig.colorbar(sc, ax=ax, label=_map_label, fraction=0.046, pad=0.04)

    def _contour_map(ax, data, title, wcs_proj=None, cmap=colormap):
        vmin_p, vmax_p = _vminmax(data)
        levels = np.linspace(vmin_p, vmax_p, 10)
        data_filled = np.where(np.isfinite(data), data, vmin_p)
        tf = ax.get_transform('pixel') if wcs_proj is not None else ax.transData
        cs = ax.contour(data_filled, levels=levels, cmap=cmap,
                        origin='lower', transform=tf)
        ax.set_title(title, fontsize=10)
        _ax_labels(ax, wcs_proj)
        return fig.colorbar(cs, ax=ax, label=_map_label, fraction=0.046, pad=0.04)

    def _plot_map(ax, data, title, wcs_proj=None, cmap=colormap):
        if hex_plot:
            return _hexplot_map(ax, data, title, wcs_proj, cmap)
        elif contour:
            return _contour_map(ax, data, title, wcs_proj, cmap)
        else:
            return _imshow_map(ax, data, title, wcs_proj, cmap)

    # Fixed y-axis range per spectrum panel, so the scale stays put when the
    # user changes the zoom/region/velocity range interactively. Each panel's
    # range is locked to whatever spectrum is first drawn in it.
    _spec_ylims = {}

    def _fixed_ylim(ax, spec):
        if ax not in _spec_ylims:
            ymin = float(np.nanmin(spec))
            ymax = float(np.nanmax(spec))
            margin = 0.05 * (ymax - ymin)
            _spec_ylims[ax] = (ymin - margin, ymax + margin)
        return _spec_ylims[ax]

    def _plot_spec(ax, vel, spec, title, color='steelblue'):
        ax.plot(vel, spec, color=color, linewidth=1.0, drawstyle='steps-mid')
        ax.set_xlabel('Velocity (km/s)')
        ax.set_ylabel('Mean T$_A^*$ (K)')
        ax.set_title(title, fontsize=10)
        ax.axhline(0, color='gray', linewidth=0.5, linestyle=':')
        if velocity_range is not None:
            ax.axvspan(v_min, v_max, alpha=0.15, color=color,
                       label=f'{v_min:.0f}–{v_max:.0f} km/s')
            ax.legend(fontsize=8)
        ax.set_ylim(_fixed_ylim(ax, spec))

    # When something restricts the rendered area — a polygon or circular aperture
    # selection, or --trim-edges peeling the outer border — shrink the map view to
    # the surviving data's bounding box (plus a small margin) so it fills the frame
    # instead of sitting as a small blob ringed by the blanked border. Cropping to
    # collapsed_display's finite pixels naturally excludes the trimmed ring, so the
    # frame follows the map rather than leaving the trimmed-away space blank.
    crop_to_selection = ((polygon is not None and len(polygon) >= 3)
                         or mask_radius_arcmin is not None
                         or trim_edges is not None)

    def _crop_axes_to_data(ax):
        finite = np.isfinite(collapsed_display)
        if not finite.any():
            return
        ys, xs = np.where(finite)
        x0, x1 = int(xs.min()), int(xs.max())
        y0, y1 = int(ys.min()), int(ys.max())
        pad_x = max(2.0, 0.05 * (x1 - x0))
        pad_y = max(2.0, 0.05 * (y1 - y0))
        ax.set_xlim(x0 - pad_x, x1 + pad_x)
        ax.set_ylim(y0 - pad_y, y1 + pad_y)

    # Full map
    _cb_full = [_plot_map(ax_full, collapsed_display, title_base, wcs_proj=wcs2d_obj)]
    if crop_to_selection:
        _crop_axes_to_data(ax_full)
    _cb_zoom = [None]

    # Zoom box on full map + zoom panel
    _zoom_rect = [None]   # mutable ref so the interactive handler can replace it
    if has_zoom:
        pix_tf_full = ax_full.get_transform('pixel') if use_wcs else ax_full.transData
        _zoom_rect[0] = ax_full.add_patch(Rectangle(
            (x0 - 0.5, y0 - 0.5), x1 - x0, y1 - y0,
            linewidth=1.5, edgecolor='white', facecolor='none',
            linestyle='--', transform=pix_tf_full))
        _cb_zoom[0] = _plot_map(ax_zoom, zoomed,
                                f"Zoom centre  {zoom_size_arcmin:.1f}′ × {zoom_size_arcmin:.1f}′",
                                wcs_proj=wcs2d_zoom)

    # Circle overlay on full map (and zoom map if active)
    if has_region:
        pix_tf_full = ax_full.get_transform('pixel') if use_wcs else ax_full.transData
        ax_full.add_patch(Circle(
            (cx_reg, cy_reg), radius_pix,
            linewidth=1.5, edgecolor='lime', facecolor='none',
            transform=pix_tf_full))
        if has_zoom and ax_zoom is not None:
            pix_tf_zoom = ax_zoom.get_transform('pixel') if use_wcs else ax_zoom.transData
            ax_zoom.add_patch(Circle(
                (cx_reg - x0, cy_reg - y0), radius_pix,
                linewidth=1.5, edgecolor='lime', facecolor='none',
                transform=pix_tf_zoom))

    # Spectrum panels
    if has_zoom:
        _plot_spec(ax_spec_zoom, channels_kms, zoom_spec,
                   f"Mean spectrum — zoom {zoom_size_arcmin:.1f}′")
    else:
        _plot_spec(ax_spec_zoom, channels_kms, full_spec,
                   full_spec_label)

    if has_region:
        _plot_spec(ax_spec_reg, channels_kms, region_spec,
                   f"Mean spectrum — {reg_label}", color='tomato')

    fig.tight_layout()

    # Update hex marker sizes now that axes have their final dimensions
    if _hex_scatters:
        fig.canvas.draw()
        for sc, spacing, ax in _hex_scatters:
            bbox = ax.get_window_extent()
            xlim = ax.get_xlim()
            ylim = ax.get_ylim()
            px_per_data = min(
                bbox.width  / max(abs(xlim[1] - xlim[0]), 1),
                bbox.height / max(abs(ylim[1] - ylim[0]), 1),
            )
            pts_per_data = px_per_data * 72 / fig.dpi
            hex_r_pt = spacing * pts_per_data / np.sqrt(3)
            s = max(1.0, (3 * np.sqrt(3) / 2) * hex_r_pt ** 2)
            sc.set_sizes([s] * len(sc.get_offsets()))

    # ------------------------------------------------------------------
    # Interactive: click on the full map's colorbar to adjust the colour
    # scale (left click = new lower clip, right click = new upper clip).
    # The same scale is applied to the zoom map's colorbar, if present.
    # ------------------------------------------------------------------
    _cb_history = []

    def _on_colorbar_click(event):
        cb_full = _cb_full[0]
        if cb_full is None or event.inaxes is not cb_full.ax or event.ydata is None:
            return
        if event.button not in (1, 3):
            return

        norm_full = cb_full.mappable.norm
        new_vmin, new_vmax = norm_full.vmin, norm_full.vmax
        if event.button == 1:
            new_vmin = float(event.ydata)
        else:
            new_vmax = float(event.ydata)
        if new_vmin >= new_vmax:
            return

        _cb_history.append((norm_full.vmin, norm_full.vmax))

        for cb in (_cb_full[0], _cb_zoom[0]):
            if cb is None:
                continue
            norm = cb.mappable.norm
            norm.vmin, norm.vmax = new_vmin, new_vmax
            cb.mappable.set_norm(norm)
            cb.update_normal(cb.mappable)

        fig.canvas.draw_idle()

    def _on_key_press(event):
        if event.key != 'b' or not _cb_history:
            return

        prev_vmin, prev_vmax = _cb_history.pop()

        for cb in (_cb_full[0], _cb_zoom[0]):
            if cb is None:
                continue
            norm = cb.mappable.norm
            norm.vmin, norm.vmax = prev_vmin, prev_vmax
            cb.mappable.set_norm(norm)
            cb.update_normal(cb.mappable)

        fig.canvas.draw_idle()

    fig.canvas.mpl_connect('button_press_event', _on_colorbar_click)
    fig.canvas.mpl_connect('key_press_event', _on_key_press)

    def _on_spec_erase(event):
        # 'b' with the cursor over a spectrum panel erases that panel and drops
        # its cached y-limits, so the next spectrum drawn there rescales from
        # scratch (the y-limits are otherwise fixed to the first spectrum shown).
        if event.key != 'b':
            return
        ax = event.inaxes
        if ax not in (ax_spec_zoom, ax_spec_reg):
            return
        _spec_ylims.pop(ax, None)
        ax.cla()
        ax.set_axis_off()
        fig.canvas.draw_idle()

    fig.canvas.mpl_connect('key_press_event', _on_spec_erase)

    # ------------------------------------------------------------------
    # Interactive handlers (always enabled): left-click any map → pixel
    # spectrum; right-click the full map twice → define a zoom region; click
    # a spectrum panel twice → set the velocity range. When no --zoom was
    # given the right-click *creates* the zoom panel on the fly (the same
    # code path that redefines an existing one), so `_zb` starts empty and
    # `ax_zoom` starts hidden until the first right-click pair.
    # ------------------------------------------------------------------
    if True:
        _click_circle = [None]   # mutable slot for the current marker patch
        # Shared zoom bounds — updated when the user (re)defines the zoom region.
        # Empty (x0==x1) until a zoom exists; _zoom_active() gates zoom-only redraws.
        if has_zoom:
            _zb = {'x0': x0, 'y0': y0, 'x1': x1, 'y1': y1}
        else:
            _zb = {'x0': 0, 'y0': 0, 'x1': 0, 'y1': 0}

        def _zoom_active():
            return _zb['x1'] > _zb['x0'] and _zb['y1'] > _zb['y0']
        # Last pixel selected by clicking on a map panel.
        _selected_pixel = {'px': None, 'py': None}
        # Current velocity range — updated by clicking the spectrum panels.
        _vr = {
            'v_min': float(v_min) if velocity_range is not None else float(channels_kms[0]),
            'v_max': float(v_max) if velocity_range is not None else float(channels_kms[-1]),
        }

        def _on_zoom_click(event):
            if event.inaxes is not ax_zoom or event.xdata is None:
                return
            px_zoom = int(round(event.xdata))
            py_zoom = int(round(event.ydata))
            px_cube = px_zoom + _zb['x0']
            py_cube = py_zoom + _zb['y0']
            if not (0 <= px_cube < nx and 0 <= py_cube < ny):
                return

            # _show_pixel_spectrum is defined later but called at event time — ok.
            _show_pixel_spectrum(px_cube, py_cube)

            # Circle on zoom map
            if _click_circle[0] is not None:
                _click_circle[0].remove()
            pix_tf = ax_zoom.get_transform('pixel') if use_wcs else ax_zoom.transData
            c = Circle((px_zoom, py_zoom), radius=0.5,
                       linewidth=2.0, edgecolor='white', facecolor='none',
                       transform=pix_tf)
            ax_zoom.add_patch(c)
            _click_circle[0] = c

            # Mirror circle on the full map
            if _full_circle[0] is not None:
                try:
                    _full_circle[0].remove()
                except Exception:
                    pass
            pix_tf_full = ax_full.get_transform('pixel') if use_wcs else ax_full.transData
            cf = Circle((px_cube, py_cube), radius=0.5,
                        linewidth=2.0, edgecolor='white', facecolor='none',
                        transform=pix_tf_full)
            ax_full.add_patch(cf)
            _full_circle[0] = cf

            fig.canvas.draw_idle()

        fig.canvas.mpl_connect('button_press_event', _on_zoom_click)

        # Interactive: full-map clicks.
        #   Left click  (button 1) — show pixel spectrum + circle (same as zoom map)
        #   Right click (button 3) — define new zoom region by two opposite corners
        _corner_clicks = []
        _corner_marker = [None]
        _full_circle   = [None]   # circle marker on the full map

        def _show_pixel_spectrum(px_cube, py_cube):
            """Shared helper: update ax_spec_reg with the spectrum at (px_cube, py_cube)."""
            _selected_pixel['px'] = px_cube
            _selected_pixel['py'] = py_cube
            spec = cube[:, py_cube, px_cube]
            ax_spec_reg.set_visible(True)
            ax_spec_reg.cla()
            ax_spec_reg.set_axis_on()
            ax_spec_reg.plot(channels_kms, spec, color='tomato', linewidth=1.0, drawstyle='steps-mid')
            ax_spec_reg.set_xlabel('Velocity (km/s)')
            ax_spec_reg.set_ylabel('T$_A^*$ (K)')
            ax_spec_reg.axhline(0, color='gray', linewidth=0.5, linestyle=':')
            ax_spec_reg.axvspan(_vr['v_min'], _vr['v_max'], alpha=0.15, color='tomato',
                                label=f"{_vr['v_min']:.0f}–{_vr['v_max']:.0f} km/s")
            ax_spec_reg.legend(fontsize=8)
            ax_spec_reg.set_ylim(_fixed_ylim(ax_spec_reg, spec))
            title_str = f'Pixel ({px_cube}, {py_cube})'
            try:
                sky = wcs2d.pixel_to_world(px_cube, py_cube)
                ra_str  = sky.ra.to_string(unit='hourangle', sep=':', precision=1)
                dec_str = sky.dec.to_string(sep=':', precision=0, alwayssign=True)
                title_str += f'\n{ra_str}  {dec_str}'
            except Exception:
                pass
            ax_spec_reg.set_title(title_str, fontsize=9)

        def _on_full_click(event):
            nonlocal ax_zoom, wcs2d_zoom
            if event.inaxes is not ax_full or event.xdata is None:
                return
            px = int(round(event.xdata))
            py = int(round(event.ydata))
            if not (0 <= px < nx and 0 <= py < ny):
                return

            if event.button == 1:
                # Left click: show pixel spectrum and mark with a circle
                _show_pixel_spectrum(px, py)

                # Circle on the full map
                if _full_circle[0] is not None:
                    _full_circle[0].remove()
                pix_tf = ax_full.get_transform('pixel') if use_wcs else ax_full.transData
                c = Circle((px, py), radius=0.5,
                           linewidth=2.0, edgecolor='white', facecolor='none',
                           transform=pix_tf)
                ax_full.add_patch(c)
                _full_circle[0] = c

                # Also mark on zoom map if the pixel falls within the zoom region
                if _click_circle[0] is not None:
                    try:
                        _click_circle[0].remove()
                    except Exception:
                        pass
                    _click_circle[0] = None
                if (_zb['x0'] <= px < _zb['x1'] and _zb['y0'] <= py < _zb['y1']):
                    pix_tf_z = ax_zoom.get_transform('pixel') if use_wcs else ax_zoom.transData
                    cz = Circle((px - _zb['x0'], py - _zb['y0']), radius=0.5,
                                linewidth=2.0, edgecolor='white', facecolor='none',
                                transform=pix_tf_z)
                    ax_zoom.add_patch(cz)
                    _click_circle[0] = cz

                fig.canvas.draw_idle()

            elif event.button == 3:
                # Right click: corner selection for zoom region redefinition
                _corner_clicks.append((px, py))

                if len(_corner_clicks) == 1:
                    if _corner_marker[0] is not None:
                        _corner_marker[0].remove()
                    pix_tf = ax_full.get_transform('pixel') if use_wcs else ax_full.transData
                    m, = ax_full.plot(px, py, '+', color='yellow',
                                      markersize=12, markeredgewidth=2,
                                      transform=pix_tf)
                    _corner_marker[0] = m
                    fig.canvas.draw_idle()

                elif len(_corner_clicks) == 2:
                    (px0, py0), (px1, py1) = _corner_clicks
                    new_x0 = max(0, min(px0, px1))
                    new_x1 = min(nx, max(px0, px1) + 1)
                    new_y0 = max(0, min(py0, py1))
                    new_y1 = min(ny, max(py0, py1) + 1)

                    _corner_clicks.clear()
                    if _corner_marker[0] is not None:
                        _corner_marker[0].remove()
                        _corner_marker[0] = None

                    if new_x1 <= new_x0 or new_y1 <= new_y0:
                        fig.canvas.draw_idle()
                        return

                    _zb['x0'] = new_x0
                    _zb['y0'] = new_y0
                    _zb['x1'] = new_x1
                    _zb['y1'] = new_y1

                    new_zoomed    = collapsed_display[new_y0:new_y1, new_x0:new_x1]
                    _new_zoom_slice = cube[:, new_y0:new_y1, new_x0:new_x1]
                    if align_peaks:
                        new_zoom_spec = _peak_aligned_mean_spectrum(
                            _new_zoom_slice.reshape(_new_zoom_slice.shape[0], -1).T,
                            channels_kms, search_range=velocity_range)
                    else:
                        new_zoom_spec = np.nanmean(_new_zoom_slice, axis=(1, 2))

                    if _zoom_rect[0] is not None:
                        _zoom_rect[0].remove()
                    pix_tf_full = ax_full.get_transform('pixel') if use_wcs else ax_full.transData
                    _zoom_rect[0] = ax_full.add_patch(Rectangle(
                        (new_x0 - 0.5, new_y0 - 0.5),
                        new_x1 - new_x0, new_y1 - new_y0,
                        linewidth=1.5, edgecolor='white', facecolor='none',
                        linestyle='--', transform=pix_tf_full))

                    if _cb_zoom[0] is not None:
                        _cb_zoom[0].remove()
                    if use_wcs:
                        # The zoom box's pixel offset changed, so the WCS
                        # reference pixel must change too. WCSAxes binds its
                        # projection at creation time, so cla() alone can't
                        # pick up the new WCS — recreate the axes instead.
                        wcs2d_zoom = _make_wcs_zoom(new_x0, new_y0)
                        fig.delaxes(ax_zoom)
                        ax_zoom = fig.add_subplot(2, 2, 2, projection=wcs2d_zoom)
                    else:
                        ax_zoom.cla()
                    ax_zoom.set_position(_ax_zoom_pos)
                    ax_zoom.set_visible(True)   # may have started hidden (no --zoom)
                    _cb_zoom[0] = _plot_map(ax_zoom, new_zoomed,
                                            f"Zoom ({new_x0}:{new_x1}, {new_y0}:{new_y1})",
                                            wcs_proj=wcs2d_zoom)

                    ax_spec_zoom.cla()
                    _plot_spec(ax_spec_zoom, channels_kms, new_zoom_spec,
                               'Mean spectrum — zoom')

                    if _click_circle[0] is not None:
                        try:
                            _click_circle[0].remove()
                        except Exception:
                            pass
                        _click_circle[0] = None

                    fig.canvas.draw_idle()

        fig.canvas.mpl_connect('button_press_event', _on_full_click)

        # Interactive: click on a spectrum panel to set velocity-range edges.
        # First click draws an orange dashed line; second click recollapses the cube.
        _vel_clicks  = []
        _vel_markers = []   # list of axvline handles (one per visible spectrum panel)

        def _recompute_velocity_range(new_v_min, new_v_max):
            _vr['v_min'] = new_v_min
            _vr['v_max'] = new_v_max
            new_chan_mask = (channels_kms >= new_v_min) & (channels_kms <= new_v_max)
            if not np.any(new_chan_mask):
                return
            new_deltav = abs(cdelt3) / 1e3
            if mode == 'peak-intensity':
                new_col = np.nanmax(cube[new_chan_mask], axis=0).astype(np.float32)
            else:
                new_col = (np.nansum(cube[new_chan_mask], axis=0) * new_deltav).astype(np.float32)
            new_col_disp = new_col.copy()
            new_col_disp[display_mask] = np.nan
            new_col_disp[zero_display_mask] = 0.0
            if suppress_negative:
                new_col_disp[np.isfinite(new_col_disp) & (new_col_disp < 0)] = 0.0
            if suppress_high is not None:
                new_col_disp[np.isfinite(new_col_disp) & (new_col_disp > suppress_high)] = suppress_high
            if smooth_sigma is not None and smooth_sigma > 0:
                from astropy.convolution import Gaussian2DKernel, convolve as _conv
                new_col_disp = _conv(new_col_disp, Gaussian2DKernel(x_stddev=smooth_sigma),
                                     boundary='fill', fill_value=np.nan,
                                     nan_treatment='interpolate',
                                     preserve_nan=True).astype(np.float32)

            new_title = f"{title_base}  [{new_v_min:.0f}–{new_v_max:.0f} km/s]"

            # Redraw full map
            if _cb_full[0] is not None:
                _cb_full[0].remove()
            ax_full.cla()
            ax_full.set_position(_ax_full_pos)
            _cb_full[0] = _plot_map(ax_full, new_col_disp, new_title, wcs_proj=wcs2d_obj)
            if crop_to_selection:
                _crop_axes_to_data(ax_full)
            pix_tf_full = ax_full.get_transform('pixel') if use_wcs else ax_full.transData
            if _zoom_active():
                _zoom_rect[0] = ax_full.add_patch(Rectangle(
                    (_zb['x0'] - 0.5, _zb['y0'] - 0.5),
                    _zb['x1'] - _zb['x0'], _zb['y1'] - _zb['y0'],
                    linewidth=1.5, edgecolor='white', facecolor='none',
                    linestyle='--', transform=pix_tf_full))
            _full_circle[0] = None
            if _selected_pixel['px'] is not None:
                c = Circle((_selected_pixel['px'], _selected_pixel['py']), radius=0.5,
                           linewidth=2.0, edgecolor='white', facecolor='none',
                           transform=pix_tf_full)
                ax_full.add_patch(c)
                _full_circle[0] = c

            # Redraw zoom map (only when a zoom region is currently active)
            if _zoom_active():
                if _cb_zoom[0] is not None:
                    _cb_zoom[0].remove()
                new_zoomed = new_col_disp[_zb['y0']:_zb['y1'], _zb['x0']:_zb['x1']]
                ax_zoom.cla()
                ax_zoom.set_position(_ax_zoom_pos)
                _cb_zoom[0] = _plot_map(ax_zoom, new_zoomed,
                                        f"Zoom ({_zb['x0']}:{_zb['x1']}, {_zb['y0']}:{_zb['y1']})",
                                        wcs_proj=wcs2d_zoom)
                _click_circle[0] = None
                if (_selected_pixel['px'] is not None and
                        _zb['x0'] <= _selected_pixel['px'] < _zb['x1'] and
                        _zb['y0'] <= _selected_pixel['py'] < _zb['y1']):
                    px_z = _selected_pixel['px'] - _zb['x0']
                    py_z = _selected_pixel['py'] - _zb['y0']
                    pix_tf_z = ax_zoom.get_transform('pixel') if use_wcs else ax_zoom.transData
                    cz = Circle((px_z, py_z), radius=0.5,
                                linewidth=2.0, edgecolor='white', facecolor='none',
                                transform=pix_tf_z)
                    ax_zoom.add_patch(cz)
                    _click_circle[0] = cz

            # Redraw mean spectrum panel with updated velocity range.
            # With a zoom active this is the zoom-region mean; otherwise it falls
            # back to the whole-map mean (matching the no-zoom initial layout).
            if _zoom_active():
                _zb_slice = cube[:, _zb['y0']:_zb['y1'], _zb['x0']:_zb['x1']]
                _zb_spec_title = 'Mean spectrum — zoom'
            else:
                _zb_slice = cube
                _zb_spec_title = 'Mean spectrum — full map'
            if align_peaks:
                new_zoom_spec = _peak_aligned_mean_spectrum(
                    _zb_slice.reshape(_zb_slice.shape[0], -1).T,
                    channels_kms, search_range=velocity_range)
            else:
                new_zoom_spec = np.nanmean(_zb_slice, axis=(1, 2))
            ax_spec_zoom.cla()
            ax_spec_zoom.plot(channels_kms, new_zoom_spec, color='steelblue', linewidth=1.0, drawstyle='steps-mid')
            ax_spec_zoom.set_xlabel('Velocity (km/s)')
            ax_spec_zoom.set_ylabel('Mean T$_A^*$ (K)')
            ax_spec_zoom.set_title(_zb_spec_title, fontsize=10)
            ax_spec_zoom.axhline(0, color='gray', linewidth=0.5, linestyle=':')
            ax_spec_zoom.axvspan(new_v_min, new_v_max, alpha=0.15, color='steelblue',
                                 label=f'{new_v_min:.0f}–{new_v_max:.0f} km/s')
            ax_spec_zoom.legend(fontsize=8)
            ax_spec_zoom.set_ylim(_fixed_ylim(ax_spec_zoom, new_zoom_spec))

            # Redraw pixel spectrum if one is selected
            if _selected_pixel['px'] is not None:
                _show_pixel_spectrum(_selected_pixel['px'], _selected_pixel['py'])

        def _on_spec_click(event):
            if event.inaxes not in (ax_spec_zoom, ax_spec_reg) or event.xdata is None:
                return
            v = float(event.xdata)
            _vel_clicks.append(v)

            if len(_vel_clicks) == 1:
                # First edge: orange dashed line on every visible spectrum panel
                for line in _vel_markers:
                    try:
                        line.remove()
                    except Exception:
                        pass
                _vel_markers.clear()
                for ax in (ax_spec_zoom, ax_spec_reg):
                    if ax.get_visible():
                        _vel_markers.append(
                            ax.axvline(v, color='orange', lw=1.5, linestyle='--', alpha=0.8))
                fig.canvas.draw_idle()

            elif len(_vel_clicks) == 2:
                new_v_min = min(_vel_clicks)
                new_v_max = max(_vel_clicks)
                _vel_clicks.clear()
                for line in _vel_markers:
                    try:
                        line.remove()
                    except Exception:
                        pass
                _vel_markers.clear()
                if new_v_min < new_v_max:
                    _recompute_velocity_range(new_v_min, new_v_max)
                fig.canvas.draw_idle()

        fig.canvas.mpl_connect('button_press_event', _on_spec_click)

        # Show a hint in the bottom-right panel when no region was pre-defined
        if not has_region:
            ax_spec_reg.set_visible(True)
            ax_spec_reg.set_axis_off()
            ax_spec_reg.text(0.5, 0.5,
                             'Left-click any map → pixel spectrum\n'
                             'Middle-click map or zoom → region vertex '
                             '(left-click closes ≥3)\n'
                             + ('Right-click full map ×2 → redefine zoom\n\n'
                                if has_zoom else
                                'Right-click full map ×2 → create zoom\n\n')
                             + 'Click spectrum ×2 → set velocity range\n'
                             "'b' → reset region",
                             ha='center', va='center',
                             transform=ax_spec_reg.transAxes,
                             fontsize=9, color='gray', style='italic')

    # ------------------------------------------------------------------
    # Interactive: draw a region polygon with the MIDDLE mouse button
    # (like --polygon / _interactive_polygon). Enabled only when no static
    # --polygon / --select-polygon region is already applied. Each middle
    # click on the full map drops a vertex; a LEFT click closes the polygon
    # (needs >=3 vertices), and only then is the mean spectrum over the
    # enclosed (rendered) pixels shown in the bottom-right panel. Press 'b'
    # to reset. Middle-clicking after a close starts a fresh polygon.
    # ------------------------------------------------------------------
    if polygon is None:
        _poly_verts = []            # accumulated (x, y) full-cube pixel vertices
        _poly_line = [None]         # Line2D currently drawn on ax_full
        _poly_line_zoom = [None]    # mirror Line2D drawn on ax_zoom (when active)
        _poly_closed = [False]      # True once the ring has been closed

        def _show_poly_hint():
            ax_spec_reg.set_visible(True)
            ax_spec_reg.cla()
            ax_spec_reg.set_axis_off()
            ax_spec_reg.text(0.5, 0.5,
                             'Middle-click the map to drop\n'
                             'region-polygon vertices;\n'
                             'left-click to close it (≥3).\n'
                             'The enclosed mean spectrum\n'
                             'is shown here.\n\n'
                             "press 'b' to reset the region",
                             ha='center', va='center',
                             transform=ax_spec_reg.transAxes,
                             fontsize=9, color='gray', style='italic')

        def _draw_poly_region_spec():
            outside = _polygon_mask(_poly_verts, ny, nx)   # True OUTSIDE polygon
            inside = (~outside) & (~display_mask)          # rendered pixels only
            n_in = int(inside.sum())
            if n_in == 0:                                  # fall back to raw polygon
                inside = ~outside
                n_in = int(inside.sum())
            if n_in == 0:
                return
            sel = cube.reshape(cube.shape[0], -1)[:, inside.ravel()]
            if align_peaks:
                spec = _peak_aligned_mean_spectrum(sel.T, channels_kms,
                                                   search_range=velocity_range)
            else:
                spec = np.nanmean(sel, axis=1)
            ax_spec_reg.set_visible(True)
            ax_spec_reg.cla()
            ax_spec_reg.set_axis_on()
            ax_spec_reg.plot(channels_kms, spec, color='magenta', linewidth=1.0,
                             drawstyle='steps-mid')
            ax_spec_reg.set_xlabel('Velocity (km/s)')
            ax_spec_reg.set_ylabel('Mean T$_A^*$ (K)')
            ax_spec_reg.axhline(0, color='gray', linewidth=0.5, linestyle=':')
            if velocity_range is not None:
                ax_spec_reg.axvspan(v_min, v_max, alpha=0.15, color='magenta')
            ax_spec_reg.set_title(f'Mean spectrum — polygon region ({n_in} px)',
                                  fontsize=9)
            finite = spec[np.isfinite(spec)]
            if finite.size:
                lo, hi = float(finite.min()), float(finite.max())
                m = 0.05 * (hi - lo) if hi > lo else 1.0
                ax_spec_reg.set_ylim(lo - m, hi + m)

        def _redraw_poly_line():
            for _ln in (_poly_line, _poly_line_zoom):
                if _ln[0] is not None:
                    try:
                        _ln[0].remove()
                    except Exception:
                        pass
                    _ln[0] = None
            if not _poly_verts:
                return
            xs = [v[0] for v in _poly_verts]
            ys = [v[1] for v in _poly_verts]
            if _poly_closed[0] and len(_poly_verts) >= 3:   # close the ring visually
                xs = xs + [xs[0]]
                ys = ys + [ys[0]]
            pix_tf = ax_full.get_transform('pixel') if use_wcs else ax_full.transData
            _poly_line[0], = ax_full.plot(xs, ys, '-o', color='magenta', lw=1.5,
                                          ms=5, mfc='yellow', transform=pix_tf)
            # Mirror the polygon onto the zoom panel (coords shifted by the zoom
            # origin) so it's visible whichever panel you're drawing on.
            if _zoom_active():
                xs_z = [x - _zb['x0'] for x in xs]
                ys_z = [y - _zb['y0'] for y in ys]
                pix_tf_z = ax_zoom.get_transform('pixel') if use_wcs else ax_zoom.transData
                _poly_line_zoom[0], = ax_zoom.plot(xs_z, ys_z, '-o', color='magenta',
                                                   lw=1.5, ms=5, mfc='yellow',
                                                   transform=pix_tf_z)

        def _on_poly_click(event):
            if event.xdata is None:
                return
            # Vertices can be dropped on the full map OR, when a zoom is active,
            # on the zoom panel — zoom-panel pixel coords are shifted back to
            # full-cube coords by the zoom origin so both share one polygon.
            if event.inaxes is ax_full:
                px, py = float(event.xdata), float(event.ydata)
            elif event.inaxes is ax_zoom and _zoom_active():
                px = float(event.xdata) + _zb['x0']
                py = float(event.ydata) + _zb['y0']
            else:
                return
            # MIDDLE click: drop a vertex (open, growing polygon).
            if event.button == 2:
                if not (0 <= px < nx and 0 <= py < ny):
                    return
                if _poly_closed[0]:                 # start a fresh polygon
                    _poly_verts.clear()
                    _poly_closed[0] = False
                _poly_verts.append((px, py))
                _redraw_poly_line()
                fig.canvas.draw_idle()
            # LEFT click: close the polygon (needs >=3 vertices) and show spectrum.
            elif event.button == 1 and not _poly_closed[0] and len(_poly_verts) >= 3:
                _poly_closed[0] = True
                _redraw_poly_line()
                _draw_poly_region_spec()
                fig.canvas.draw_idle()

        def _on_poly_reset(event):
            if event.key != 'b' or not _poly_verts:
                return
            _poly_verts.clear()
            _poly_closed[0] = False
            _redraw_poly_line()
            _show_poly_hint()
            fig.canvas.draw_idle()

        fig.canvas.mpl_connect('button_press_event', _on_poly_click)
        fig.canvas.mpl_connect('key_press_event', _on_poly_reset)

        # NB: the initial bottom-right hint is now drawn once by the always-on
        # interactive block above (the comprehensive hint, when no --region was
        # given). We only fall back to the polygon-only hint after a 'b' reset
        # (handled in _on_poly_reset), so nothing is drawn here to avoid
        # overwriting that comprehensive hint.

    if plot_output:
        fig.savefig(plot_output, dpi=150, bbox_inches='tight')
        print(f"  Plot saved: {plot_output}")

    return collapsed, header2d, fig


# ══════════════════════════════════════════════════════════════════════════════
# compare_maps — collapse N cubes to maps and render them in one comparison grid
# ══════════════════════════════════════════════════════════════════════════════
# These reproduce the *compute* half of collapse_cube (load → shuffle → collapse →
# coverage/trim/all-nan masks → suppress → smooth), deliberately WITHOUT the large
# interactive/plotting half. Kept separate rather than refactored out of
# collapse_cube so that function (1400+ lines, heavily interactive) is left
# untouched. If the collapse maths there changes, mirror it here too.

def _trim_border_px(trim_edges, pixel_scale_arcsec, nx, ny):
    """Parse a --trim-edges string to (border_px, human_desc). Mirrors collapse_cube.

    Bare number / "N%" → N% of the smaller map axis; a unit suffix
    ("arcmin"/"am"/"'" or "arcsec"/"as"/'"') → a fixed angular border.
    Returns (0, None) if unparseable.
    """
    s = str(trim_edges).strip().lower()
    try:
        if s.endswith('arcmin') or s.endswith('am') or s.endswith("'"):
            val = float(s.rstrip("'").replace('arcmin', '').replace('am', ''))
            border_px = int(round(val * 60.0 / pixel_scale_arcsec))
            desc = f"{val:g}' ({border_px}px)"
        elif s.endswith('arcsec') or s.endswith('as') or s.endswith('"'):
            val = float(s.rstrip('"').replace('arcsec', '').replace('as', ''))
            border_px = int(round(val / pixel_scale_arcsec))
            desc = f'{val:g}" ({border_px}px)'
        else:
            pct = float(s.rstrip('%'))
            border_px = int(round(min(nx, ny) * pct / 100.0))
            desc = f"{pct:g}% ({border_px}px)"
    except ValueError:
        return 0, None
    return border_px, desc


def _footprint_trim_ring(footprint, border_px):
    """Ring of footprint pixels to blank when eroding the footprint inward by
    border_px (8-connected), following the irregular coverage boundary. Interior
    holes are filled first (flood-fill from the array border) so only the OUTER
    edge is peeled. Mirrors collapse_cube's edge-trim block.
    """
    from collections import deque
    ny, nx = footprint.shape
    bg = ~footprint
    exterior = np.zeros_like(footprint)
    dq = deque()
    for x in range(nx):
        for yb in (0, ny - 1):
            if bg[yb, x] and not exterior[yb, x]:
                exterior[yb, x] = True; dq.append((yb, x))
    for y in range(ny):
        for xb in (0, nx - 1):
            if bg[y, xb] and not exterior[y, xb]:
                exterior[y, xb] = True; dq.append((y, xb))
    while dq:
        y, x = dq.popleft()
        for dy, dx in ((1, 0), (-1, 0), (0, 1), (0, -1)):
            yy2, xx2 = y + dy, x + dx
            if 0 <= yy2 < ny and 0 <= xx2 < nx and bg[yy2, xx2] and not exterior[yy2, xx2]:
                exterior[yy2, xx2] = True; dq.append((yy2, xx2))
    filled = ~exterior
    eroded = filled
    for _ in range(border_px):
        e = eroded.copy()
        e[:-1, :] &= eroded[1:, :]; e[1:, :] &= eroded[:-1, :]
        e[:, :-1] &= eroded[:, 1:]; e[:, 1:] &= eroded[:, :-1]
        e[:-1, :-1] &= eroded[1:, 1:]; e[:-1, 1:] &= eroded[1:, :-1]
        e[1:, :-1] &= eroded[:-1, 1:]; e[1:, 1:] &= eroded[:-1, :-1]
        e[0, :] = False; e[-1, :] = False; e[:, 0] = False; e[:, -1] = False
        eroded = e
    return footprint & ~eroded


def _finish_weight_map_2d(tag, weight_map, coverage_map, header,
                          coverage_threshold=0.0, trim_edges=None,
                          suppress_negative=False, suppress_high=None,
                          smooth_sigma=None, verbose=True):
    """Apply the display masks/scaling to a 2D WEIGHT_MAP and return like
    collapse_cube_to_map. Used when a cube carries only a 2D per-spectrum
    WEIGHT_MAP (--weight-spectra) rather than a 3D WEIGHT_CUBE."""
    ny, nx = weight_map.shape
    disp = weight_map.astype(np.float32).copy()
    display_mask = ~np.isfinite(disp)

    if coverage_threshold > 0 and coverage_map is not None and np.nanmax(coverage_map) > 0:
        display_mask |= coverage_map < coverage_threshold * np.nanmax(coverage_map)

    if trim_edges is not None:
        pix_as = abs(float(header.get('CDELT2', 1.0))) * 3600.0
        border_px, desc = _trim_border_px(trim_edges, pix_as, nx, ny)
        if desc is not None and border_px > 0:
            footprint = np.isfinite(weight_map)
            if coverage_map is not None:
                footprint &= (coverage_map > 0)
            ring = _footprint_trim_ring(footprint, border_px)
            if int(ring.sum()) < int(footprint.sum()):
                display_mask |= ring
            elif verbose:
                print(f"  [{tag}] --trim-edges {trim_edges} would blank everything; ignored")

    disp[display_mask] = np.nan
    if suppress_negative:
        disp[np.isfinite(disp) & (disp < 0)] = 0.0
    if suppress_high is not None:
        disp[np.isfinite(disp) & (disp > suppress_high)] = suppress_high
    if smooth_sigma is not None and smooth_sigma > 0:
        from astropy.convolution import Gaussian2DKernel, convolve as _convolve
        disp = _convolve(disp, Gaussian2DKernel(x_stddev=smooth_sigma),
                         boundary='fill', fill_value=np.nan,
                         nan_treatment='interpolate', preserve_nan=True).astype(np.float32)

    header2d = WCS(header).dropaxis(2).to_header()
    header2d['BUNIT'] = 'weight'
    if verbose:
        finite = np.isfinite(disp)
        rng = (f"[{np.nanmin(disp):.2f}, {np.nanmax(disp):.2f}]" if finite.any() else "[empty]")
        print(f"  [{tag}] weights (2D WEIGHT_MAP): {ny}×{nx}, "
              f"{int(finite.sum())} valid px, range {rng} weight")
    return disp, header2d, 'Per-spectrum weight (Σ K·w)', 'weight', 'weights'


def collapse_cube_to_map(
    cube_fits: str,
    velocity_range: Optional[Tuple[float, float]] = None,
    mode: str = 'moment0',
    peak_range_channels: int = 5,
    peak_smooth_channels: Optional[float] = None,
    shuffle: bool = False,
    shuffle_snr: float = 5.0,
    shuffle_field_smooth: float = 6.0,
    shuffle_window_kms: float = 20.0,
    shuffle_ref_velocity: Optional[float] = None,
    coverage_threshold: float = 0.0,
    trim_edges: Optional[str] = None,
    suppress_negative: bool = False,
    suppress_high: Optional[float] = None,
    smooth_sigma: Optional[float] = None,
    snr_threshold: Optional[float] = None,
    collapse_weights: bool = False,
    verbose: bool = True,
):
    """Collapse ONE gridded cube to a 2D display map — compute only, no plotting.

    Applies the same steps as collapse_cube: velocity-range select, optional
    velocity-field shuffle, moment0 / peak-range-int / peak-intensity collapse,
    coverage + edge-trim + all-NaN masking (masked → NaN), optional per-cube
    SNR threshold (sub-threshold pixels shown as 0), suppress-negative,
    suppress-high, and Gaussian display smoothing.

    With ``collapse_weights=True`` the weights are shown instead of the primary
    science cube: the 3D ``WEIGHT_CUBE`` extension is used if present (channel sum
    of the gridded weights Σ K·w, no channel-width factor), otherwise it falls
    back to the 2D ``WEIGHT_MAP`` (per-spectrum weights from --weight-spectra, the
    analog of CLASS's .wei file) shown directly. Line/peak-specific options
    (shuffle, peak modes/smoothing) are ignored. The velocity axis is read from
    the primary header (WEIGHT_CUBE shares the same geometry).

    Returns
    -------
    (display_map, header2d, map_label, map_units, mode)
    """
    tag = Path(cube_fits).name
    with fits.open(cube_fits) as hdul:
        header = hdul[0].header.copy()
        if collapse_weights:
            ext_names = [e.name for e in hdul]
            cov_ext = (hdul['COVERAGE'].data.astype(float)
                       if 'COVERAGE' in ext_names else None)
            if 'WEIGHT_CUBE' in ext_names:
                cube = hdul['WEIGHT_CUBE'].data.astype(float)
            elif 'WEIGHT_MAP' in ext_names:
                # Per-spectrum weighting (no per-channel weights) records a single
                # 2D WEIGHT_MAP (Σ K·w per pixel) — the analog of CLASS's .wei
                # file. There is no spectral axis to collapse, so display it
                # directly, still honouring the coverage/edge/suppress/smooth
                # display controls.
                return _finish_weight_map_2d(
                    tag, hdul['WEIGHT_MAP'].data.astype(float), cov_ext, header,
                    coverage_threshold=coverage_threshold, trim_edges=trim_edges,
                    suppress_negative=suppress_negative, suppress_high=suppress_high,
                    smooth_sigma=smooth_sigma, verbose=verbose)
            else:
                raise ValueError(
                    f"{tag}: no WEIGHT_CUBE or WEIGHT_MAP extension. Rebuild the cube "
                    f"with --weight-spectra (writes a 2D WEIGHT_MAP) or with "
                    f"--weight-channels --create-weights-datacube (writes a 3D "
                    f"WEIGHT_CUBE) to record the weights used.")
        else:
            cube = hdul[0].data.astype(float)
    if cube.ndim != 3:
        raise ValueError(f"{tag}: expected a 3-D datacube, got shape {cube.shape}")
    nvel, ny, nx = cube.shape

    crval3 = float(header.get('CRVAL3', 0.0))
    crpix3 = float(header.get('CRPIX3', 1.0))
    cdelt3 = float(header.get('CDELT3', 1.0))
    channels_kms = (crval3 + (np.arange(nvel) - (crpix3 - 1)) * cdelt3) / 1e3
    deltav_kms = abs(cdelt3) / 1e3

    if velocity_range is not None:
        v_min, v_max = float(velocity_range[0]), float(velocity_range[1])
        chan_mask = (channels_kms >= v_min) & (channels_kms <= v_max)
        if not np.any(chan_mask):
            raise ValueError(f"{tag}: no channels in [{v_min}, {v_max}] km/s "
                             f"(cube covers {channels_kms.min():.1f}–{channels_kms.max():.1f}).")
    else:
        chan_mask = np.ones(nvel, dtype=bool)

    # --- Optional velocity-field shuffle (same as collapse_cube) ---
    if shuffle and not collapse_weights:
        if velocity_range is None:
            raise ValueError(f"{tag}: --shuffle requires --velocity-range.")
        v0_raw, v0_field, _snr = _velocity_field_from_cube(
            cube, channels_kms, chan_mask.copy(), deltav_kms, shuffle_snr, shuffle_field_smooth)
        if int(np.sum(np.isfinite(v0_raw))) == 0:
            raise ValueError(f"{tag}: --shuffle found no pixels at S/N ≥ {shuffle_snr}; "
                             f"lower --shuffle-snr or widen --velocity-range.")
        ref_v = (float(shuffle_ref_velocity) if shuffle_ref_velocity is not None
                 else float(np.nanmedian(v0_raw)))
        cube = _shuffle_cube(cube, channels_kms, v0_field, ref_v)
        v_min, v_max = ref_v - shuffle_window_kms, ref_v + shuffle_window_kms
        chan_mask = (channels_kms >= v_min) & (channels_kms <= v_max)
        if verbose:
            print(f"  [{tag}] shuffle → ref {ref_v:.1f} km/s, window "
                  f"{v_min:.1f}–{v_max:.1f} ({int(chan_mask.sum())} ch)")

    # --- Peak-location smoothing (locate on smoothed, read raw) ---
    cube_pf = cube
    if (not collapse_weights and peak_smooth_channels and peak_smooth_channels > 0
            and mode in ('peak-intensity', 'peak-range-int')):
        cube_pf = _spectral_smooth_nanaware(cube, float(peak_smooth_channels))

    # --- Collapse ---
    if collapse_weights:  # plain Σ K·w over channels, no channel-width factor
        collapsed = np.nansum(cube[chan_mask, :, :], axis=0).astype(np.float32)
        units, label = 'weight', 'Integrated weight (Σ K·w)'
    elif mode == 'peak-intensity':
        peak_idx = np.argmax(np.where(np.isnan(cube_pf), -np.inf, cube_pf), axis=0)
        _ry = np.arange(ny)[:, None]; _rx = np.arange(nx)[None, :]
        collapsed = cube[peak_idx, _ry, _rx].astype(np.float32)
        units, label = 'K', 'Peak intensity (K)'
    elif mode == 'peak-range-int':
        sub = cube[chan_mask, :, :]; sub_pf = cube_pf[chan_mask, :, :]
        nsel = sub.shape[0]
        peak_idx = np.argmax(np.where(np.isnan(sub_pf), -np.inf, sub_pf), axis=0)
        n = peak_range_channels
        acc = np.zeros((ny, nx), dtype=np.float64)
        ri = np.arange(ny)[:, None]; ci = np.arange(nx)[None, :]
        for k in range(-n, n + 1):
            ck = peak_idx + k
            valid = (ck >= 0) & (ck < nsel)
            vals = sub[np.clip(ck, 0, nsel - 1), ri, ci]
            acc += np.where(valid & ~np.isnan(vals), vals, 0.0)
        collapsed = (acc * deltav_kms).astype(np.float32)
        units, label = 'K km/s', f'Peak ±{n}ch (K km/s)'
    else:  # moment0
        collapsed = (np.nansum(cube[chan_mask, :, :], axis=0) * deltav_kms).astype(np.float32)
        units, label = 'K km/s', 'K km/s'

    # --- Masks (masked → NaN in the display) ---
    display_mask = np.zeros((ny, nx), dtype=bool)
    coverage_map = None
    if coverage_threshold > 0:
        try:
            with fits.open(cube_fits) as _h:
                if 'COVERAGE' in [e.name for e in _h]:
                    coverage_map = _h['COVERAGE'].data.astype(float)
        except Exception:
            pass
        if coverage_map is not None and coverage_map.max() > 0:
            display_mask |= coverage_map < coverage_threshold * coverage_map.max()

    if trim_edges is not None:
        pix_as = abs(float(header.get('CDELT2', 1.0))) * 3600.0
        border_px, desc = _trim_border_px(trim_edges, pix_as, nx, ny)
        if desc is not None and border_px > 0:
            footprint = ~np.all(np.isnan(cube[chan_mask, :, :]), axis=0)
            if coverage_map is not None:
                footprint &= (coverage_map > 0)
            ring = _footprint_trim_ring(footprint, border_px)
            if int(ring.sum()) < int(footprint.sum()):
                display_mask |= ring
            elif verbose:
                print(f"  [{tag}] --trim-edges {trim_edges} would blank everything; ignored")

    display_mask |= np.all(np.isnan(cube[chan_mask, :, :]), axis=0)

    # --- SNR threshold (per cube, from its own line-free channels) ---
    # Noise map = std over channels OUTSIDE --velocity-range; SNR = signal / noise
    # (peak: peak/noise; integrated: moment0 / (noise·sqrt(n_chan)·deltav)).
    # Sub-threshold pixels are shown as 0 (the zero colour), matching collapse_cube,
    # not blanked to NaN. Requires --velocity-range and does not apply to weights.
    snr_mask = None
    if snr_threshold is not None and not collapse_weights:
        noise_chan_mask = ~chan_mask
        if int(noise_chan_mask.sum()) >= 2:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", category=RuntimeWarning)
                noise_map = np.nanstd(cube[noise_chan_mask, :, :], axis=0)
                n_signal_chan = int(chan_mask.sum())
                if mode == 'peak-intensity':
                    snr_map = collapsed / noise_map
                else:
                    snr_map = collapsed / (noise_map * np.sqrt(n_signal_chan) * deltav_kms)
            snr_mask = np.isfinite(snr_map) & (snr_map < snr_threshold)
            if verbose:
                print(f"  [{tag}] SNR threshold {snr_threshold:.1f}: "
                      f"{int(np.sum(snr_mask & ~display_mask))} pixels below cutoff "
                      f"shown as 0 ({int(noise_chan_mask.sum())} line-free channels)")
        elif verbose:
            print(f"  [{tag}] --snr-threshold ignored: need --velocity-range to "
                  f"identify line-free channels for the noise map")

    disp = collapsed.copy()
    disp[display_mask] = np.nan
    if snr_mask is not None:
        disp[snr_mask & ~display_mask] = 0.0
    if suppress_negative:
        disp[np.isfinite(disp) & (disp < 0)] = 0.0
    if suppress_high is not None:
        disp[np.isfinite(disp) & (disp > suppress_high)] = suppress_high
    if smooth_sigma is not None and smooth_sigma > 0:
        from astropy.convolution import Gaussian2DKernel, convolve as _convolve
        disp = _convolve(disp, Gaussian2DKernel(x_stddev=smooth_sigma),
                         boundary='fill', fill_value=np.nan,
                         nan_treatment='interpolate', preserve_nan=True).astype(np.float32)

    effective_mode = 'weights' if collapse_weights else mode
    header2d = WCS(header).dropaxis(2).to_header()
    header2d['BUNIT'] = units
    if verbose:
        finite = np.isfinite(disp)
        rng = (f"[{np.nanmin(disp):.2f}, {np.nanmax(disp):.2f}]" if finite.any() else "[empty]")
        print(f"  [{tag}] {effective_mode}: {ny}×{nx}, {int(finite.sum())} valid px, range {rng} {units}")
    return disp, header2d, label, units, effective_mode


def _compare_build_norm(stretch, vmin, vmax, gamma, data=None):
    """Same stretch→matplotlib-norm mapping as collapse_cube._make_norm."""
    from matplotlib import colors as mcolors
    if stretch == 'sqrt':
        return mcolors.PowerNorm(gamma=0.5, vmin=vmin, vmax=vmax)
    if stretch == 'power':
        return mcolors.PowerNorm(gamma=gamma, vmin=vmin, vmax=vmax)
    if stretch == 'log':
        safe_vmin = max(vmin, 1e-6 * vmax) if vmax > 0 else 1e-6
        return mcolors.LogNorm(vmin=safe_vmin, vmax=vmax)
    if stretch == 'symlog':
        linthresh = max(abs(vmax), abs(vmin)) * 0.1 or 1e-6
        return mcolors.SymLogNorm(linthresh=linthresh, vmin=vmin, vmax=vmax)
    if stretch == 'asinh':
        if hasattr(mcolors, 'AsinhNorm'):
            return mcolors.AsinhNorm(linear_width=(vmax - vmin) * 0.1 or 1e-6, vmin=vmin, vmax=vmax)
        return mcolors.PowerNorm(gamma=0.5, vmin=vmin, vmax=vmax)
    return mcolors.Normalize(vmin=vmin, vmax=vmax)


def _crop_box_max_empty(mask, tol):
    """Largest-ish axis-aligned box with empty fraction <= tol.

    Greedily trims whichever current border (top / bottom / left / right) is the
    emptiest, one line at a time, until the box's empty (no-data) fraction is at
    or below ``tol``. RA and Dec are therefore trimmed by independent amounts.
    ``tol=0`` shrinks to a fully-covered rectangle (no empty pixels), biting into
    the rounded data edges as needed.

    Parameters
    ----------
    mask : 2D bool ndarray
        Coverage mask (True = has data).
    tol : float
        Maximum allowed empty fraction of the final box (0..1).

    Returns
    -------
    (y0, y1, x0, x1) inclusive index bounds.
    """
    ny, nx = mask.shape
    y0, y1, x0, x1 = 0, ny - 1, 0, nx - 1
    while y1 > y0 and x1 > x0:
        sub = mask[y0:y1 + 1, x0:x1 + 1]
        if (1.0 - sub.mean()) <= tol:
            break
        # Empty fraction of each border line; trim the worst offender.
        top = 1.0 - mask[y0, x0:x1 + 1].mean()
        bot = 1.0 - mask[y1, x0:x1 + 1].mean()
        lft = 1.0 - mask[y0:y1 + 1, x0].mean()
        rgt = 1.0 - mask[y0:y1 + 1, x1].mean()
        worst = max(top, bot, lft, rgt)
        if worst == top:
            y0 += 1
        elif worst == bot:
            y1 -= 1
        elif worst == lft:
            x0 += 1
        else:
            x1 -= 1
    return y0, y1, x0, x1


def compare_maps(
    cube_files: List[str],
    plot_output: Optional[str] = None,
    titles: Optional[List[str]] = None,
    use_wcs: bool = False,
    colormap: str = 'rainbow',
    stretch: str = 'linear',
    gamma: float = 1.0,
    percentile_clip: Optional[Tuple[float, float]] = None,
    per_panel_scale: bool = False,
    interpolation: str = 'nearest',
    tighten: bool = False,
    crop_empty: bool = False,
    show: bool = True,
    **collapse_kwargs,
):
    """Collapse N cubes and draw their maps in one auto-sized comparison grid.

    All the collapse controls (velocity_range, mode, peak_range_channels,
    peak_smooth_channels, shuffle*, coverage_threshold, trim_edges,
    suppress_negative, suppress_high, smooth_sigma) are forwarded to
    collapse_cube_to_map via **collapse_kwargs. By default every panel shares one
    colour scale + colorbar (comparable brightness); per_panel_scale=True scales
    each to its own range instead.

    tighten=True lays all panels in a single row, pushed together horizontally
    (wspace=0) with a shared y-axis (only the leftmost keeps y ticks/label), and
    drops both the per-panel titles and the figure suptitle for a compact strip.

    crop_empty controls trimming of empty (no-coverage) regions: True drops only
    the fully-empty outer rows/columns (margins only — keeps every data pixel and
    the rounded corners); a number crops harder until at most that percent of the
    frame is empty, trimming RA and Dec by independent amounts and biting into the
    rounded data edges. None/False = no cropping.
    """
    if not cube_files:
        raise ValueError("compare_maps: no cube files given")

    maps = []
    for f in cube_files:
        disp, hdr2d, label, units, mode = collapse_cube_to_map(f, **collapse_kwargs)
        maps.append((f, disp, hdr2d, label, units, mode))

    # --- Optional crop of empty (no-coverage) regions from the frame.
    #   crop_empty is True  -> margins only: drop rows/cols that are ENTIRELY
    #                          empty, trimming the outer white border while
    #                          keeping every data pixel (and the rounded corners).
    #   crop_empty is a number -> crop harder until at most that PERCENT of the
    #                          frame is empty; the emptiest border is trimmed
    #                          first, so RA and Dec shrink by independent amounts
    #                          (this bites into the rounded data edges).
    # Coverage is the union across panels, so they stay aligned. Skipped on shape
    # mismatch. ---
    if crop_empty is not None and crop_empty is not False:
        shapes = {m[1].shape for m in maps}
        if len(shapes) == 1:
            finite_any = np.zeros(maps[0][1].shape, dtype=bool)
            for m in maps:
                finite_any |= np.isfinite(m[1])
            if crop_empty is True:
                cols = np.where(finite_any.any(axis=0))[0]
                rows = np.where(finite_any.any(axis=1))[0]
                y0, y1 = (int(rows[0]), int(rows[-1])) if rows.size else (0, -1)
                x0, x1 = (int(cols[0]), int(cols[-1])) if cols.size else (0, -1)
                _desc = "margins only"
            else:
                _pct = max(0.0, float(crop_empty))
                y0, y1, x0, x1 = _crop_box_max_empty(finite_any, _pct / 100.0)
                _desc = f"≤{_pct:.0f}% empty"
            if x1 > x0 and y1 > y0:
                cropped = []
                for (f, disp, hdr2d, label, units, m_) in maps:
                    h = hdr2d.copy()
                    if 'CRPIX1' in h:
                        h['CRPIX1'] = float(h['CRPIX1']) - x0
                    if 'CRPIX2' in h:
                        h['CRPIX2'] = float(h['CRPIX2']) - y0
                    cropped.append((f, disp[y0:y1 + 1, x0:x1 + 1], h, label, units, m_))
                maps = cropped
                _emp = 100.0 * (1.0 - np.isfinite(maps[0][1]).mean())
                print(f"  Crop ({_desc}): {maps[0][1].shape[1]}×{maps[0][1].shape[0]} "
                      f"({_emp:.0f}% empty remaining)")
            else:
                warnings.warn("compare_maps: crop collapsed the frame; skipped.",
                              UserWarning)
        else:
            warnings.warn("compare_maps: crop skipped — panels have different "
                          "shapes, cannot crop to a common footprint.", UserWarning)

    N = len(maps)
    mode = maps[0][5]
    map_label = maps[0][3]
    suppress_negative = bool(collapse_kwargs.get('suppress_negative', False))

    # --- Auto grid: near-square for N>3, single row for N<=3 (looks nice).
    # tighten forces a single row so panels can be pushed together horizontally.
    import math
    if tighten or N <= 3:
        ncols, nrows = N, 1
    else:
        ncols = int(math.ceil(math.sqrt(N)))
        nrows = int(math.ceil(N / ncols))
    # constrained_layout keeps the (possibly wrapped, multi-line) panel titles,
    # the suptitle and the colorbar from overlapping — without it a two-line
    # title grows up into the panel/frame above. It is disabled for tighten so
    # subplots_adjust(wspace=0) can butt the panels together.
    _tighten_layout = None
    if tighten:
        # Inch-based layout so the equal-aspect panels butt together (wspace=0)
        # AND the colorbar gets its own right-hand strip (never over a panel).
        # Each panel's box is sized to the data aspect ratio; the figure width is
        # the sum of panel widths plus fixed margins for the left y-labels and the
        # colorbar strip.
        _ny0, _nx0 = maps[0][1].shape
        _aspect = (_nx0 / _ny0) if _ny0 else 1.0
        _ph = 4.3                       # panel (axes) height, inches
        _pw = _ph * _aspect             # panel width per aspect ratio
        _left, _bot, _top = 0.6, 0.55, 0.2   # margins for labels
        _cgap, _cw, _rpad = 0.12, 0.16, 0.45  # colorbar gap / width / right pad
        _figw = _left + ncols * _pw + _cgap + _cw + _rpad
        _figh = _ph + _bot + _top
        _tighten_layout = dict(ph=_ph, pw=_pw, left=_left, bot=_bot, top=_top,
                               cgap=_cgap, cw=_cw, figw=_figw, figh=_figh, ncols=ncols)
        fig = plt.figure(figsize=(_figw, _figh), constrained_layout=False)
    else:
        fig = plt.figure(figsize=(5.2 * ncols, 5.0 * nrows),
                         constrained_layout=True)

    # --- Shared colour limits from pooled data (unless per-panel) ---
    def _vlims(arrs):
        pooled = np.concatenate([a[np.isfinite(a)].ravel() for a in arrs
                                 if np.isfinite(a).any()]) if arrs else np.array([np.nan])
        lo, hi = percentile_clip if percentile_clip is not None else (0, 100)
        vmin = 0.0 if (suppress_negative and mode == 'moment0') else float(np.nanpercentile(pooled, lo))
        vmax = float(np.nanpercentile(pooled, hi))
        if not np.isfinite(vmax) or vmax <= vmin:
            vmax = vmin + 1.0
        return vmin, vmax

    shared_vmin, shared_vmax = _vlims([m[1] for m in maps])
    shared_norm = _compare_build_norm(stretch, shared_vmin, shared_vmax, gamma,
                                      data=np.concatenate([m[1][np.isfinite(m[1])].ravel() for m in maps]))
    cmap_obj = plt.colormaps[colormap].copy() if isinstance(colormap, str) else colormap.copy()
    cmap_obj.set_bad('white')

    # Wrap long panel titles onto multiple lines so they stay within their
    # panel and are not clipped/hidden by the colorbar (the rightmost column's
    # centred title would otherwise overflow rightward under the shared
    # colorbar). Breaks at spaces or underscores — filenames have no spaces —
    # keeping each delimiter on the preceding line.
    def _wrap_title(s, width=24):
        if len(s) <= width:
            return s
        lines, cur, seg = [], '', ''
        for ch in s:
            seg += ch
            if ch in ' _':
                if cur and len(cur) + len(seg) > width:
                    lines.append(cur); cur = seg
                else:
                    cur += seg
                seg = ''
        if seg:
            if cur and len(cur) + len(seg) > width:
                lines.append(cur); cur = seg
            else:
                cur += seg
        if cur:
            lines.append(cur)
        return '\n'.join(lines)

    # Each colorbar is registered with the image(s) it controls, so an
    # interactive click on it can re-clip those panels (see below).
    _cb_registry = []
    axes, images, last_im = [], [], None
    for i, (f, disp, hdr2d, label, units, _m) in enumerate(maps):
        # tighten: non-WCS panels share the leftmost panel's y-axis so they align
        # and only the first needs y ticks. (WCSAxes can't take a sharey kwarg, so
        # there we just hide the repeated Dec labels below.)
        share_kw = {}
        if tighten and not use_wcs and axes:
            share_kw['sharey'] = axes[0]
        if use_wcs:
            ax = fig.add_subplot(nrows, ncols, i + 1, projection=WCS(hdr2d))
        else:
            ax = fig.add_subplot(nrows, ncols, i + 1, **share_kw)
        if per_panel_scale:
            vmn, vmx = _vlims([disp])
            norm = _compare_build_norm(stretch, vmn, vmx, gamma, data=disp[np.isfinite(disp)])
        else:
            norm = shared_norm
        im = ax.imshow(disp, origin='lower', cmap=cmap_obj, norm=norm, interpolation=interpolation)
        last_im = im
        if not tighten:
            ttl = titles[i] if (titles and i < len(titles)) else Path(f).stem
            ax.set_title(_wrap_title(ttl), fontsize=9)
        hide_y = tighten and i > 0  # keep y ticks/label only on the leftmost panel
        if use_wcs:
            ax.set_xlabel('RA')
            if hide_y:
                try:
                    ax.coords[1].set_ticklabel_visible(False)
                    ax.coords[1].set_axislabel('')
                except Exception:
                    pass
            else:
                ax.set_ylabel('Dec')
            try:
                ax.coords.grid(color='gray', alpha=0.25, linewidth=0.5)
            except Exception:
                pass
            # Drop tick labels that would overlap — re-evaluated on every draw,
            # so RA/Dec labels re-space (don't stack) when the view is zoomed.
            try:
                ax.coords[0].set_ticklabel(exclude_overlapping=True)
                ax.coords[1].set_ticklabel(exclude_overlapping=True)
            except Exception:
                pass
        else:
            ax.set_xlabel('pixel')
            if hide_y:
                ax.tick_params(labelleft=False)
                ax.set_ylabel('')
            else:
                ax.set_ylabel('pixel')
        if per_panel_scale:
            cb = fig.colorbar(im, ax=ax, label=map_label, fraction=0.046, pad=0.04)
            _cb_registry.append({'cb': cb, 'images': [im]})
        axes.append(ax)
        images.append(im)

    from matplotlib.cm import ScalarMappable
    if tighten:
        # Exact fractional margins (from the inch layout) so panels touch at
        # wspace=0 and the shared colorbar lives in its own reserved right strip.
        L = _tighten_layout
        left_f = L['left'] / L['figw']
        right_f = (L['left'] + L['ncols'] * L['pw']) / L['figw']
        bot_f = L['bot'] / L['figh']
        top_f = 1.0 - L['top'] / L['figh']
        fig.subplots_adjust(left=left_f, right=right_f, bottom=bot_f, top=top_f, wspace=0)
        if not per_panel_scale and last_im is not None:
            sm = ScalarMappable(norm=shared_norm, cmap=cmap_obj); sm.set_array([])
            cax = fig.add_axes([right_f + L['cgap'] / L['figw'], bot_f,
                                L['cw'] / L['figw'], top_f - bot_f])
            cb = fig.colorbar(sm, cax=cax, label=map_label)
            _cb_registry.append({'cb': cb, 'images': list(images)})
    else:
        if not per_panel_scale and last_im is not None:
            sm = ScalarMappable(norm=shared_norm, cmap=cmap_obj); sm.set_array([])
            cb = fig.colorbar(sm, ax=axes, label=map_label, fraction=0.046, pad=0.04)
            _cb_registry.append({'cb': cb, 'images': list(images)})
        scale_note = 'per-panel scale' if per_panel_scale else 'shared scale'
        fig.suptitle(f"{mode} comparison — {N} cubes ({nrows}×{ncols}, {scale_note})",
                     fontsize=12, fontweight='bold')
    # Layout is handled by constrained_layout (non-tighten) or the explicit
    # subplots_adjust above (tighten); tight_layout would conflict with both.

    # ------------------------------------------------------------------
    # Interactive: click on a colorbar to adjust its colour scale
    # (left click = new lower clip, right click = new upper clip). In the
    # default shared-scale mode the single colorbar re-clips every panel at
    # once; with per_panel_scale each colorbar re-clips only its own panel.
    # Press 'b' to undo the last change.
    # ------------------------------------------------------------------
    _cb_history = []

    def _apply(entry, vmin, vmax):
        cb = entry['cb']
        cb.mappable.norm.vmin, cb.mappable.norm.vmax = vmin, vmax
        cb.mappable.set_norm(cb.mappable.norm)
        cb.update_normal(cb.mappable)
        for im in entry['images']:
            im.norm.vmin, im.norm.vmax = vmin, vmax
            im.set_norm(im.norm)
        fig.canvas.draw_idle()

    def _on_colorbar_click(event):
        if event.button not in (1, 3) or event.ydata is None:
            return
        for entry in _cb_registry:
            if event.inaxes is not entry['cb'].ax:
                continue
            norm = entry['cb'].mappable.norm
            new_vmin, new_vmax = norm.vmin, norm.vmax
            if event.button == 1:
                new_vmin = float(event.ydata)
            else:
                new_vmax = float(event.ydata)
            if new_vmin >= new_vmax:
                return
            _cb_history.append((entry, norm.vmin, norm.vmax))
            _apply(entry, new_vmin, new_vmax)
            return

    # --- Interactive region select: RIGHT-click corners on any panel to define a
    # rectangle; all panels then zoom to the SAME SKY region. Clicks are converted
    # to world coords (RA/Dec) via the clicked panel's WCS and back to each other
    # panel's pixels via its own WCS, so the panels correspond even when their
    # pixel grids differ. Right-click ≥2 corners (e.g. two opposite corners, or
    # all four) — the zoom is the bounding box of the clicked corners. Press 'r'
    # to reset. Right-clicks on a colorbar still adjust the colour scale. ---
    panel_wcs = []
    for (_f, _d, _h, _l, _u, _m) in maps:
        try:
            panel_wcs.append(WCS(_h).celestial)
        except Exception:
            panel_wcs.append(None)
    _orig_lims = [(ax.get_xlim(), ax.get_ylim()) for ax in axes]
    _corners = []   # world (ra, dec) when use_wcs else pixel (x, y)
    _markers = []   # marker artists showing the clicked corners

    def _reset_zoom():
        for ax, (xl, yl) in zip(axes, _orig_lims):
            ax.set_xlim(*xl); ax.set_ylim(*yl)
        for mk in _markers:
            try:
                mk.remove()
            except Exception:
                pass
        _markers.clear(); _corners.clear()
        fig.canvas.draw_idle()

    def _apply_region():
        if len(_corners) < 2:
            return
        for j, ax in enumerate(axes):
            if use_wcs and panel_wcs[j] is not None:
                pts = [panel_wcs[j].world_to_pixel_values(ra, dec)
                       for (ra, dec) in _corners]
                xs = [float(p[0]) for p in pts]; ys = [float(p[1]) for p in pts]
            else:
                xs = [c[0] for c in _corners]; ys = [c[1] for c in _corners]
            ax.set_xlim(min(xs), max(xs)); ax.set_ylim(min(ys), max(ys))
        fig.canvas.draw_idle()

    def _clear_markers():
        for mk in _markers:
            try:
                mk.remove()
            except Exception:
                pass
        _markers.clear()

    def _on_region_click(event):
        # Right-click on a panel drops a corner. (Right-click on a colorbar is
        # handled separately by _on_colorbar_click.) A rectangle is two opposite
        # corners: the first click shows a marker, the second zooms all panels
        # and the markers are cleared. The next right-click starts a fresh
        # selection in the (now zoomed) view, so you can keep zooming further.
        if event.button != 3 or event.inaxes not in axes or event.xdata is None:
            return
        j = axes.index(event.inaxes)
        _markers.append(event.inaxes.plot(event.xdata, event.ydata, 'x',
                                          color='magenta', ms=9, mew=2)[0])
        if use_wcs and panel_wcs[j] is not None:
            ra, dec = panel_wcs[j].pixel_to_world_values(event.xdata, event.ydata)
            _corners.append((float(ra), float(dec)))
        else:
            _corners.append((float(event.xdata), float(event.ydata)))
        if len(_corners) >= 2:
            _apply_region()          # zoom all panels
            _clear_markers()         # remove the corner ✕ marks
            _corners.clear()         # ready for the next selection
        fig.canvas.draw_idle()

    def _on_key_press(event):
        if event.key in ('r', 'R'):
            _reset_zoom()
            return
        if event.key != 'b' or not _cb_history:
            return
        entry, prev_vmin, prev_vmax = _cb_history.pop()
        _apply(entry, prev_vmin, prev_vmax)

    fig.canvas.mpl_connect('button_press_event', _on_colorbar_click)
    fig.canvas.mpl_connect('button_press_event', _on_region_click)
    fig.canvas.mpl_connect('key_press_event', _on_key_press)
    if show:
        print("  Interactive: right-click corners on any panel to zoom all panels "
              "to that region; press 'r' to reset.")

    if plot_output:
        fig.savefig(plot_output, dpi=150, bbox_inches='tight')
        print(f"  Comparison figure saved: {plot_output}")
    if show:
        plt.show()
    return fig


def _celestial_header_dict(cube_fits: str) -> dict:
    """Plain 2-D celestial WCS keys for a gridded cube, as cygrid wants them.

    Read straight off the file rather than round-tripped through
    ``WCS.to_header()``: that emits a PC matrix and drops NAXISn, whereas
    cygrid.WcsGrid wants the same CDELT-style keys create_wcs_header writes.
    """
    hdr = fits.getheader(cube_fits)
    out = {'NAXIS': 2, 'NAXIS1': int(hdr['NAXIS1']), 'NAXIS2': int(hdr['NAXIS2'])}
    for key in ('CTYPE1', 'CTYPE2', 'CRVAL1', 'CRVAL2',
                'CRPIX1', 'CRPIX2', 'CDELT1', 'CDELT2'):
        out[key] = hdr[key]
    for key in ('EQUINOX', 'RADESYS'):
        if key in hdr:
            out[key] = hdr[key]
    return out


def regrid_map_with_cygrid(src_map: np.ndarray,
                           src_cube_fits: str,
                           dst_cube_fits: str,
                           kernel_fwhm_arcsec: Optional[float] = None,
                           verbose: bool = True) -> np.ndarray:
    """Resample a 2-D map onto another cube's celestial grid using cygrid.

    The source map's pixel centres are treated as irregular samples and pushed
    through :func:`grid_to_map`, the same cygrid primitive create_datacube
    uses — no second interpolation engine, and NaN pixels drop out of the
    sample list rather than bleeding into their neighbours.

    The default kernel FWHM is ONE TARGET PIXEL, not the beam. These maps are
    already beam-convolved products; re-gridding with the beam again would
    smooth the contour map relative to the image it is drawn over, so the
    contours would sit systematically wider than the colour scale beneath
    them. A ~pixel kernel moves the samples onto the new grid without
    materially changing the resolution. Pass an explicit value to override.
    """
    src_hdr = fits.getheader(src_cube_fits)
    src_wcs = WCS(src_hdr).celestial
    ny, nx = src_map.shape

    yy, xx = np.mgrid[0:ny, 0:nx]
    sky = src_wcs.pixel_to_world(xx.ravel(), yy.ravel())
    ras = np.asarray(sky.ra.deg, dtype=np.float64)
    decs = np.asarray(sky.dec.deg, dtype=np.float64)
    vals = np.asarray(src_map, dtype=np.float64).ravel()

    good = np.isfinite(ras) & np.isfinite(decs) & np.isfinite(vals)
    if not good.any():
        raise ValueError(f"{Path(src_cube_fits).name}: contour map is entirely blank; "
                         "nothing to regrid.")
    ras, decs, vals = ras[good], decs[good], vals[good]

    dst = _celestial_header_dict(dst_cube_fits)
    pix_deg = abs(float(dst['CDELT2']))
    fwhm_deg = (kernel_fwhm_arcsec / 3600.0) if kernel_fwhm_arcsec else pix_deg

    grid_kw = dict(
        wcs_header_dict=dst,
        naxis1=dst['NAXIS1'], naxis2=dst['NAXIS2'],
        beamsize_deg=fwhm_deg,
        ra_min=float(ras.min()), ra_max=float(ras.max()),
        dec_min=float(decs.min()), dec_max=float(decs.max()),
        kernel_fwhm_deg=fwhm_deg,
    )
    out = grid_to_map(ras=ras, decs=decs, values=vals, **grid_kw)

    # The kernel has 3*FWHM support, so gridding spreads signal a few pixels
    # PAST the source footprint, and grid_to_map only blanks pixels with
    # literally zero coverage. Left as is, contours get drawn over sky where
    # the contour cube has no data at all.
    #
    # The footprint has to be measured geometrically, not by re-gridding a
    # mask: cygrid returns a weight-normalised mean, so gridding an array of
    # ones yields ~1.0 wherever anything contributed and carries no coverage
    # information to threshold on. Map each target pixel centre back to the
    # nearest source pixel instead and keep it only if that pixel held data.
    dst_wcs = WCS(fits.getheader(dst_cube_fits)).celestial
    tny, tnx = int(dst['NAXIS2']), int(dst['NAXIS1'])
    tyy, txx = np.mgrid[0:tny, 0:tnx]
    back = src_wcs.world_to_pixel(dst_wcs.pixel_to_world(txx, tyy))
    sj = np.rint(np.asarray(back[0])).astype(int)
    si = np.rint(np.asarray(back[1])).astype(int)
    inside = (si >= 0) & (si < ny) & (sj >= 0) & (sj < nx)
    footprint = np.zeros((tny, tnx), dtype=bool)
    footprint[inside] = np.isfinite(src_map[si[inside], sj[inside]])
    out = np.where(footprint, out, np.nan).astype(np.float32)

    if verbose:
        engine = 'cygrid' if HAS_CYGRID else 'scipy fallback'
        print(f"  regrid ({engine}): {ny}×{nx} → {out.shape[0]}×{out.shape[1]}, "
              f"kernel FWHM {fwhm_deg * 3600:.2f}\", "
              f"{int(np.isfinite(out).sum())} valid px "
              f"(from {int(good.sum())} source px)")
    return out


def _overlay_levels(data: np.ndarray,
                    levels: Optional[List[float]],
                    level_fractions: Optional[List[float]],
                    n_levels: int,
                    reference_percentile: float = 99.0,
                    verbose: bool = True) -> List[float]:
    """Resolve contour levels: explicit values, else fractions of a reference.

    The reference defaults to a high PERCENTILE rather than the map maximum.
    Gridded [CII] maps routinely carry a handful of extreme edge/artefact
    pixels — in M82 the data-release moment-0 peaks 3.6x its own 99th
    percentile, with only ~1% of pixels above 20% of the max — so fractions of
    the true maximum bunch every contour onto those few pixels and the galaxy
    itself is left untraced. Set reference_percentile=100 for the literal peak.
    """
    finite = data[np.isfinite(data)]
    if finite.size == 0:
        raise ValueError("regridded contour map has no finite pixels; "
                         "check --velocity-range and the two cubes' overlap.")
    if levels:
        out = sorted(float(v) for v in levels)
    else:
        pct = float(np.clip(reference_percentile, 0.0, 100.0))
        ref = float(np.percentile(finite, pct))
        peak = float(np.nanmax(finite))
        if ref <= 0:
            ref = peak
        fracs = (sorted(float(f) for f in level_fractions) if level_fractions
                 else list(np.linspace(0.2, 0.9, n_levels)))
        out = [f * ref for f in fracs]
        if verbose:
            how = "max" if pct >= 100 else f"p{pct:g}"
            print(f"  contour reference {ref:.3g} ({how}; map max {peak:.3g}); "
                  f"levels at {', '.join(f'{f * 100:.0f}%' for f in fracs)}")
    uniq = sorted(set(out))
    if len(uniq) < len(out) and verbose:
        print(f"  dropped {len(out) - len(uniq)} duplicate level(s)")
    return uniq


def overlay_maps(
    image_cube: str,
    contour_cube: str,
    velocity_range: Optional[Tuple[float, float]] = None,
    mode: str = 'moment0',
    peak_range_channels: int = 5,
    peak_smooth_channels: Optional[float] = None,
    shuffle: bool = False,
    shuffle_snr: float = 5.0,
    shuffle_field_smooth: float = 6.0,
    shuffle_window_kms: float = 20.0,
    shuffle_ref_velocity: Optional[float] = None,
    coverage_threshold: float = 0.0,
    trim_edges: Optional[str] = None,
    suppress_negative: bool = False,
    suppress_high: Optional[float] = None,
    smooth_sigma: Optional[float] = None,
    levels: Optional[List[float]] = None,
    level_fractions: Optional[List[float]] = None,
    n_levels: int = 8,
    level_reference_percentile: float = 99.0,
    regrid_kernel_arcsec: Optional[float] = None,
    contour_color: str = 'white',
    contour_linewidth: float = 0.9,
    label_contours: bool = False,
    colormap: str = 'rainbow',
    stretch: str = 'linear',
    gamma: float = 1.0,
    percentile_clip: Optional[Tuple[float, float]] = None,
    titles: Optional[List[str]] = None,
    plot_output: Optional[str] = None,
    show: bool = True,
    verbose: bool = True,
):
    """Overlay two gridded cubes on ONE axes: first as colour, second as contours.

    Both cubes are collapsed to 2-D with the *same* settings via
    collapse_cube_to_map, so the two layers always represent the same
    quantity over the same velocity window. The contour cube is then
    resampled onto the image cube's celestial grid with cygrid, so the two
    layers genuinely share a grid — they are not merely drawn on top of each
    other and assumed to line up.

    That regrid is not optional bookkeeping. Two cubes from this pipeline can
    share a pixel scale and projection and still be offset by a fraction of a
    pixel (differing CRPIX/CRVAL) and differ in NAXIS2, in which case a raw
    ``ax.contour(other_map)`` silently draws the contours shifted against the
    colour scale.

    Note the velocity window is applied in km/s on each cube's own axis, never
    by channel index: cubes resampled by filter_fits --velocity-resample are
    anchored to multiples of the channel width relative to 0 km/s, so an
    externally-produced cube can share the channel width yet sit a fractional
    channel off. Selecting by velocity is immune to that; selecting by channel
    index is not.

    Returns
    -------
    matplotlib.figure.Figure
    """
    from matplotlib import colors as mcolors  # noqa: F401  (used via _compare_build_norm)

    for tag, path in (('image', image_cube), ('contour', contour_cube)):
        if not os.path.exists(path):
            raise FileNotFoundError(f"{tag} cube not found: {path}")

    collapse_kw = dict(
        velocity_range=velocity_range, mode=mode,
        peak_range_channels=peak_range_channels,
        peak_smooth_channels=peak_smooth_channels,
        shuffle=shuffle, shuffle_snr=shuffle_snr,
        shuffle_field_smooth=shuffle_field_smooth,
        shuffle_window_kms=shuffle_window_kms,
        shuffle_ref_velocity=shuffle_ref_velocity,
        coverage_threshold=coverage_threshold, trim_edges=trim_edges,
        suppress_negative=suppress_negative, suppress_high=suppress_high,
        smooth_sigma=smooth_sigma, verbose=verbose,
    )

    if verbose:
        print(f"Image   (colour)  : {Path(image_cube).name}")
    img, img_hdr2d, img_label, img_units, _ = collapse_cube_to_map(image_cube, **collapse_kw)
    if verbose:
        print(f"Contour (overlay) : {Path(contour_cube).name}")
    con, _con_hdr2d, con_label, con_units, _ = collapse_cube_to_map(contour_cube, **collapse_kw)

    if con_units != img_units and verbose:
        print(f"  ! unit mismatch: image is {img_units}, contours are {con_units}")

    con_r = regrid_map_with_cygrid(con, contour_cube, image_cube,
                                   kernel_fwhm_arcsec=regrid_kernel_arcsec,
                                   verbose=verbose)
    if con_r.shape != img.shape:
        raise ValueError(f"regrid produced {con_r.shape}, expected {img.shape}")

    lvls = _overlay_levels(con_r, levels, level_fractions, n_levels,
                           reference_percentile=level_reference_percentile,
                           verbose=verbose)

    finite = img[np.isfinite(img)]
    if finite.size == 0:
        raise ValueError(f"{Path(image_cube).name}: collapsed image map is entirely blank.")
    if percentile_clip:
        p_lo, p_hi = float(percentile_clip[0]), float(percentile_clip[1])
        lo = float(np.nanpercentile(finite, p_lo))
        hi = float(np.nanpercentile(finite, p_hi))
    else:
        lo, hi = float(np.nanmin(finite)), float(np.nanmax(finite))
    if suppress_negative:
        lo = 0.0
    norm = _compare_build_norm(stretch, lo, hi, gamma, data=finite)

    img_title = titles[0] if titles and len(titles) > 0 else Path(image_cube).name
    con_title = titles[1] if titles and len(titles) > 1 else Path(contour_cube).name

    fig = plt.figure(figsize=(9.0, 7.5))
    ax = fig.add_subplot(111, projection=WCS(img_hdr2d))
    im = ax.imshow(img, origin='lower', cmap=colormap, norm=norm, interpolation='nearest')
    cs = ax.contour(con_r, levels=lvls, colors=contour_color,
                    linewidths=contour_linewidth, alpha=0.9)
    if label_contours:
        ax.clabel(cs, inline=True, fontsize=6, fmt='%.3g')

    cbar = fig.colorbar(im, ax=ax, pad=0.02)
    cbar.set_label(f"{img_label}  —  {img_title}")

    ax.set_xlabel('RA (J2000)')
    ax.set_ylabel('Dec (J2000)')
    vr = (f"{velocity_range[0]:.0f}–{velocity_range[1]:.0f} km/s"
          if velocity_range else "full band")
    ax.set_title(f"{mode}, {vr}\ncolour: {img_title}   contours: {con_title}", fontsize=10)
    ax.coords.grid(color='grey', ls=':', alpha=0.5)

    fig.tight_layout()
    if plot_output:
        fig.savefig(plot_output, dpi=150, bbox_inches='tight')
        print(f"  Overlay figure saved: {plot_output}")
    if show:
        plt.show()
    return fig
