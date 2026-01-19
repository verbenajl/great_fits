"""
Core gridding and mapping functions for spectral data using WCS and cygrid.

Provides functionality to create proper WCS-based spatial maps from spectral 
observations using cygrid for optimal gridding with Gaussian kernel weighting.
Falls back to scipy griddata if cygrid is not available.
"""

from typing import Optional, Tuple, Dict, Any
import numpy as np
from astropy.io import fits
from astropy import units as u
from astropy.wcs import WCS
from astropy.table import Table
import matplotlib.pyplot as plt
import warnings

# Try to import cygrid for optimal gridding
HAS_CYGRID = False
try:
    import cygrid
    HAS_CYGRID = True
except (ImportError, ValueError) as e:
    # ValueError can occur due to numpy binary incompatibility
    warnings.warn(f"cygrid not available ({type(e).__name__}) - using scipy griddata instead", UserWarning)


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
    from oi_zeigt.basic_io import get_config
    
    # Load config
    cfg = get_config(config_path)
    gridding_cfg = cfg.get('gridding', {})
    
    # Get beamsize from config (in arcseconds) or use provided value
    if beamsize_deg is None:
        beamsize_arcsec = gridding_cfg.get('beamsize_arcsec', 15.0)  # Default 15 arcsec
        beamsize_deg = beamsize_arcsec / 3600.0  # Convert arcseconds to degrees
    
    # Get pixel size from config (in arcseconds) or calculate it
    if pixsize_deg is None:
        pixsize_arcsec = gridding_cfg.get('pixel_size_arcsec', None)
        if pixsize_arcsec is None or pixsize_arcsec == 0:
            # Calculate as beamsize / 3 (good for Nyquist sampling)
            pixsize_deg = beamsize_deg / 3.0
        else:
            pixsize_deg = pixsize_arcsec / 3600.0  # Convert arcseconds to degrees
    
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
    header['CTYPE3'] = 'VRAD'
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
    header['TELESCOP'] = 'SOFIA/GREAT'
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
    # Get the binary table
    table = hdul[1].data
    
    # Extract reference velocity (should be same for all observations)
    if 'VELOCITY' in table.names:
        velo_ref = float(table['VELOCITY'][0])
    else:
        warnings.warn("VELOCITY column not found - using 0.0 as reference", UserWarning)
        velo_ref = 0.0
    
    # Extract velocity spacing per channel
    if 'DELTAV' in table.names:
        deltav = float(table['DELTAV'][0])
    else:
        warnings.warn("DELTAV column not found - using 1.0 m/s as default", UserWarning)
        deltav = 1.0
    
    # Extract rest frequency
    if 'RESTFREQ' in table.names:
        restfreq = float(table['RESTFREQ'][0])
    else:
        warnings.warn("RESTFREQ column not found - using 1.0 Hz as default", UserWarning)
        restfreq = 1.0
    
    # Extract velocity definition
    if 'VELDEF' in table.names:
        veldef_val = table['VELDEF'][0]
        # Handle bytes vs string
        veldef = veldef_val.decode().strip() if isinstance(veldef_val, bytes) else str(veldef_val).strip()
    else:
        warnings.warn("VELDEF column not found - using 'RADI-LSR' as default", UserWarning)
        veldef = 'RADI-LSR'
    
    # Get number of spectral channels from SPECTRUM
    if 'SPECTRUM' in table.names:
        nchans = table['SPECTRUM'].shape[1]
    else:
        raise ValueError("SPECTRUM column not found in FITS table")
    
    # Extract reference pixel for spectral axis
    # First try to get it from the binary table header (CRPIX1)
    table_header = hdul[1].header
    if 'CRPIX1' in table_header:
        crpix1_spec = float(table_header['CRPIX1'])
    else:
        # Try primary header
        primary_header = hdul[0].header if len(hdul) > 0 else None
        if primary_header and 'CRPIX1' in primary_header:
            crpix1_spec = float(primary_header['CRPIX1'])
        else:
            # Default to first pixel if not found
            warnings.warn("CRPIX1 not found in headers - using 1.0 (first channel) as reference pixel", UserWarning)
            crpix1_spec = 1.0
    
    return velo_ref, deltav, restfreq, veldef, nchans, crpix1_spec


def _get_celestial_coords(hdul: fits.HDUList) -> Tuple[np.ndarray, np.ndarray]:
    """
    Extract celestial RA and Dec coordinates from FITS data.
    
    Uses CDELT2/CRVAL2 for RA and CDELT3/CRVAL3 for Dec, as these represent
    the actual celestial coordinates of the observations. The LONGITUDE/LATITUDE
    columns represent telescope position, not celestial position.
    
    Follows the formula:
        RA = CRVAL2 + CDELT2
        Dec = CRVAL3 + CDELT3
    
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
    
    # Extract delta values from data columns
    # RA = CRVAL2 + CDELT2 (in degrees)
    # Dec = CRVAL3 + CDELT3 (in degrees)
    ras = crval2 + data['CDELT2']
    decs = crval3 + data['CDELT3']
    
    return ras, decs


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
    
    # Grid using cygrid if available, otherwise use scipy
    if HAS_CYGRID:
        kernelsize_sigma = beamsize_deg / 2.355
        kernel_type = 'gauss1d'
        kernel_params = (kernelsize_sigma,)
        kernel_support = 3.0 * kernelsize_sigma
        hpx_maxres = kernelsize_sigma / 2.0
        
        # Create WCS header for cygrid
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
        
        # Create cygrid gridder with WCS header
        gridder = cygrid.WcsGrid(wcs_header_dict)
        gridder.set_kernel(
            kernel_type,
            kernel_params,
            kernel_support,
            hpx_maxres,
        )
        
        # Add data points
        gridder.grid(ras, decs, values)
        
        # Get the gridded data
        grid_map = gridder.get_datacube()  # Returns 2D array (naxis2, naxis1)
    else:
        # Fallback to scipy
        from scipy.interpolate import griddata
        
        ra_grid = np.linspace(ra_min, ra_max, naxis1)
        dec_grid = np.linspace(dec_min, dec_max, naxis2)
        ra_mesh, dec_mesh = np.meshgrid(ra_grid, dec_grid)
        
        points = np.column_stack([ras, decs])
        grid_map = griddata(points, values, (ra_mesh, dec_mesh), method='linear')
    
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
                         show_scatter: bool = False) -> Tuple[np.ndarray, fits.Header, plt.Figure]:
    """
    Create a spatial map from integrated spectral intensity.
    
    Integrates spectra across all channels, then grids with proper WCS.
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
    
    # Compute integrated intensity
    integrated = np.nansum(spectra, axis=1)
    
    # Filter valid points
    valid = np.isfinite(ras) & np.isfinite(decs) & np.isfinite(integrated)
    ras = ras[valid]
    decs = decs[valid]
    integrated = integrated[valid]
    
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
    
    # Grid using cygrid if available, otherwise use scipy
    if HAS_CYGRID:
        kernelsize_sigma = beamsize_deg / 2.355
        kernel_type = 'gauss1d'
        kernel_params = (kernelsize_sigma,)
        kernel_support = 3.0 * kernelsize_sigma
        hpx_maxres = kernelsize_sigma / 2.0
        
        # Create WCS header for cygrid
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
        
        # Create cygrid gridder with WCS header
        gridder = cygrid.WcsGrid(wcs_header_dict)
        gridder.set_kernel(
            kernel_type,
            kernel_params,
            kernel_support,
            hpx_maxres,
        )
        
        # Add data points
        gridder.grid(ras, decs, integrated)
        
        # Get the gridded data
        grid_map = gridder.get_datacube()  # Returns 2D array (naxis2, naxis1)
    else:
        # Fallback to scipy
        from scipy.interpolate import griddata
        
        ra_grid = np.linspace(ra_min, ra_max, naxis1)
        dec_grid = np.linspace(dec_min, dec_max, naxis2)
        ra_mesh, dec_mesh = np.meshgrid(ra_grid, dec_grid)
        
        points = np.column_stack([ras, decs])
        grid_map = griddata(points, integrated, (ra_mesh, dec_mesh), method='linear')
    
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
    ax.set_title(f'Integrated Intensity Map ({gridding_method}, beam={beamsize_deg*3600:.1f}″)',
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
                     spectral_params: Optional[Dict[str, Any]] = None) -> None:
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
    primary_hdu = fits.PrimaryHDU(data=grid_map, header=header)
    
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


def create_spectral_datacube(hdul: fits.HDUList,
                            beamsize_deg: float = 0.25,
                            pixsize: Optional[float] = None,
                            object_filter: Optional[str] = None,
                            figsize: Tuple[float, float] = (12, 5),
                            output_file: Optional[str] = None) -> Tuple[np.ndarray, fits.Header, plt.Figure]:
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
    
    # Extract data
    table = hdul[1].data
    spectra = table['SPECTRUM']  # Shape: (nobs, nvel)
    
    # Extract celestial coordinates (RA/Dec from CDELT2/3 and CRVAL2/3)
    ras, decs = _get_celestial_coords(hdul)
    
    nobs, nvel = spectra.shape
    
    # Filter by object if specified
    if object_filter:
        if 'OBJECT' in table.names:
            # Use substring matching: object_filter should be contained in OBJECT field
            objects = table['OBJECT']
            mask = np.array([
                object_filter.upper() in (obj.decode().strip() if isinstance(obj, bytes) else str(obj).strip()).upper()
                for obj in objects
            ])
            if not np.any(mask):
                raise ValueError(f"No observations found matching object filter '{object_filter}'")
            spectra = spectra[mask]
            ras = ras[mask]
            decs = decs[mask]
            nobs = len(ras)
            print(f"Filtered to {nobs} observations matching '{object_filter}'")
        else:
            warnings.warn("OBJECT column not found in FITS table - not filtering", UserWarning)
    else:
        print(f"Using all {nobs} observations")
    
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
    datacube = np.full((nvel, naxis2, naxis1), np.nan, dtype=np.float64)
    
    # Grid each velocity channel
    print(f"Creating datacube: {nvel} channels × {naxis2} × {naxis1} pixels...")
    
    for ichannel in range(nvel):
        if ichannel % max(1, nvel // 10) == 0:
            print(f"  Channel {ichannel+1}/{nvel}...", end='\r')
        
        channel_data = spectra[:, ichannel]
        
        # Filter valid points
        valid = np.isfinite(ras) & np.isfinite(decs) & np.isfinite(channel_data)
        if not np.any(valid):
            continue
        
        ras_valid = ras[valid]
        decs_valid = decs[valid]
        data_valid = channel_data[valid]
        
        if HAS_CYGRID:
            try:
                # Use cygrid with proper kernel parameters
                kernelsize_sigma = beamsize_deg / 2.355
                kernel_type = 'gauss1d'
                kernel_params = (kernelsize_sigma,)
                kernel_support = 3.0 * kernelsize_sigma
                hpx_maxres = kernelsize_sigma / 2.0
                
                # Create WCS header for cygrid
                wcs_header_dict = create_wcs_header(
                    naxis1=naxis1, naxis2=naxis2, naxis3=1,
                    crval_ra=ra_center, crval_dec=dec_center,
                    crpix1=crpix1, crpix2=crpix2, crpix3=crpix1_spec,
                    cdelt1=-pixsize, cdelt2=pixsize, cdelt3=1.0,
                    beamsize_deg=beamsize_deg,
                    velo_lsr=velo_lsr,
                    restfreq=restfreq,
                    object_name=object_filter or "",
                )
                
                # Create cygrid gridder with WCS header
                gridder = cygrid.WcsGrid(wcs_header_dict)
                gridder.set_kernel(
                    kernel_type,
                    kernel_params,
                    kernel_support,
                    hpx_maxres,
                )
                
                # Grid data
                gridder.grid(ras_valid, decs_valid, data_valid)
                
                # Get the gridded data
                channel_map = gridder.get_datacube()  # Returns 2D array (naxis2, naxis1)
            except Exception:
                # Fallback to scipy
                from scipy.interpolate import griddata
                ra_grid = np.linspace(ra_min, ra_max, naxis1)
                dec_grid = np.linspace(dec_min, dec_max, naxis2)
                ra_mesh, dec_mesh = np.meshgrid(ra_grid, dec_grid)
                points = np.column_stack([ras_valid, decs_valid])
                channel_map = griddata(points, data_valid, (ra_mesh, dec_mesh), method='linear')
        else:
            # Use scipy
            from scipy.interpolate import griddata
            ra_grid = np.linspace(ra_min, ra_max, naxis1)
            dec_grid = np.linspace(dec_min, dec_max, naxis2)
            ra_mesh, dec_mesh = np.meshgrid(ra_grid, dec_grid)
            points = np.column_stack([ras_valid, decs_valid])
            channel_map = griddata(points, data_valid, (ra_mesh, dec_mesh), method='linear')
        
        datacube[ichannel] = channel_map
    
    print(f"  Channel {nvel}/{nvel}... Done!")
    
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
    )
    
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
        save_map_to_fits(datacube, wcs_header, output_file, beam_maj_deg=beamsize_deg,
                        spectral_params=spectral_params)
    
    return datacube, wcs_header, fig

