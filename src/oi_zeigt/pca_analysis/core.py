"""
Core PCA analysis functions for spectral datacubes.

Provides principal component analysis tools for:
- Decomposing spectral datacubes into principal components
- Analyzing spectral variability
- Noise characterization
- Component-to-channel mapping
"""

from typing import Optional, Tuple, Dict, Any
import numpy as np
from astropy.io import fits
import warnings

# Try to import sklearn for PCA
try:
    from sklearn.decomposition import PCA
    from sklearn.preprocessing import StandardScaler
    HAS_SKLEARN = True
except ImportError:
    HAS_SKLEARN = False
    warnings.warn("scikit-learn not available - PCA functionality limited", UserWarning)


def perform_spectral_pca(datacube: np.ndarray, 
                        n_components: Optional[int] = None,
                        standardize: bool = True,
                        mask: Optional[np.ndarray] = None) -> Dict[str, Any]:
    """
    Perform Principal Component Analysis on a spectral datacube.
    
    Decomposes a datacube into principal components, useful for:
    - Understanding spectral variability
    - Dimensionality reduction
    - Noise characterization
    - Identifying dominant spectral features
    
    Parameters
    ----------
    datacube : np.ndarray
        3D datacube array of shape (nvel, ndec, nra)
    n_components : int, optional
        Number of principal components to compute. If None, uses min(nvel, ndec*nra)
    standardize : bool, optional
        Whether to standardize the data before PCA. Default is True.
    mask : np.ndarray, optional
        2D boolean mask for spatial pixels to include. If None, uses all pixels.
    
    Returns
    -------
    dict
        Dictionary containing:
        - 'components': Principal component spectra (n_components × nvel)
        - 'loadings': Component loadings for each pixel (ndec × nra × n_components)
        - 'variance_explained': Variance explained by each component
        - 'cumulative_variance': Cumulative variance explained
        - 'n_components': Number of components computed
        - 'scaler': StandardScaler object (if standardize=True)
    
    Notes
    -----
    This function reshapes the datacube into a 2D array (nvel, npixels) for PCA,
    then reshapes results back to original spatial dimensions.
    
    Examples
    --------
    >>> from astropy.io import fits
    >>> from oi_zeigt.pca_analysis import perform_spectral_pca
    >>> 
    >>> # Load datacube
    >>> hdul = fits.open('datacube.fits')
    >>> datacube = hdul[0].data
    >>> 
    >>> # Perform PCA with 10 components
    >>> results = perform_spectral_pca(datacube, n_components=10)
    >>> 
    >>> # Access results
    >>> components = results['components']           # (10, 1264)
    >>> loadings = results['loadings']               # (12, 11, 10)
    >>> var_exp = results['variance_explained']      # (10,)
    >>> cum_var = results['cumulative_variance']     # (10,)
    """
    if not HAS_SKLEARN:
        raise ImportError("scikit-learn is required for PCA analysis. "
                         "Install with: pip install scikit-learn")
    
    nvel, ndec, nra = datacube.shape
    
    # Apply mask if provided
    if mask is not None:
        if mask.shape != (ndec, nra):
            raise ValueError(f"Mask shape {mask.shape} doesn't match spatial dims ({ndec}, {nra})")
        pixels_to_use = np.where(mask.flatten())[0]
    else:
        pixels_to_use = np.arange(ndec * nra)
        mask = np.ones((ndec, nra), dtype=bool)
    
    # Reshape datacube to 2D: (nvel, npixels)
    datacube_2d = datacube.reshape(nvel, -1)[:, pixels_to_use]
    
    # Determine number of components
    if n_components is None:
        n_components = min(nvel, len(pixels_to_use))
    elif n_components > min(nvel, len(pixels_to_use)):
        warnings.warn(f"Requested {n_components} components but only "
                     f"{min(nvel, len(pixels_to_use))} possible. Using max possible.",
                     UserWarning)
        n_components = min(nvel, len(pixels_to_use))
    
    # Standardize if requested
    scaler = None
    if standardize:
        scaler = StandardScaler()
        datacube_2d = scaler.fit_transform(datacube_2d.T).T
    
    # Perform PCA
    pca = PCA(n_components=n_components)
    loadings_2d = pca.fit_transform(datacube_2d.T)  # (npixels, n_components)
    components = pca.components_  # (n_components, nvel)
    
    # Reshape loadings back to spatial dimensions
    loadings_full = np.zeros((ndec * nra, n_components), dtype=np.float32)
    loadings_full[pixels_to_use, :] = loadings_2d
    loadings = loadings_full.reshape(ndec, nra, n_components)
    
    # Calculate variance explained
    var_explained = pca.explained_variance_ratio_
    cum_var = np.cumsum(var_explained)
    
    results = {
        'components': components,
        'loadings': loadings,
        'variance_explained': var_explained,
        'cumulative_variance': cum_var,
        'n_components': n_components,
        'scaler': scaler,
        'mask': mask,
        'pca_object': pca,
    }
    
    return results


def reconstruct_from_pca(results: Dict[str, Any], 
                        n_components: Optional[int] = None) -> np.ndarray:
    """
    Reconstruct datacube from PCA results.
    
    Parameters
    ----------
    results : dict
        Dictionary returned by perform_spectral_pca()
    n_components : int, optional
        Number of components to use for reconstruction. If None, uses all.
    
    Returns
    -------
    np.ndarray
        Reconstructed datacube of same shape as original
    """
    components = results['components']
    loadings = results['loadings']
    mask = results['mask']
    scaler = results.get('scaler')
    
    if n_components is None:
        n_components = components.shape[0]
    
    # Use subset of components
    components_subset = components[:n_components]
    loadings_subset = loadings[:, :, :n_components]
    
    # Reconstruct: datacube = loadings @ components
    ndec, nra, _ = loadings_subset.shape
    nvel = components_subset.shape[1]
    
    # Reshape loadings to 2D
    loadings_2d = loadings_subset.reshape(-1, n_components)
    
    # Reconstruct
    datacube_2d = loadings_2d @ components_subset
    
    # Inverse standardization if needed
    if scaler is not None:
        datacube_2d = scaler.inverse_transform(datacube_2d.T).T
    
    # Reshape back to 3D
    datacube_reconstructed = datacube_2d.reshape(ndec, nra, nvel)
    
    # Transpose to (nvel, ndec, nra)
    datacube_reconstructed = np.transpose(datacube_reconstructed, (2, 0, 1))
    
    return datacube_reconstructed


def get_component_stats(results: Dict[str, Any]) -> Dict[str, Any]:
    """
    Get statistics about PCA components.
    
    Parameters
    ----------
    results : dict
        Dictionary returned by perform_spectral_pca()
    
    Returns
    -------
    dict
        Dictionary with statistics for each component
    """
    n_components = results['n_components']
    var_explained = results['variance_explained']
    cum_var = results['cumulative_variance']
    components = results['components']
    
    stats = {
        'n_components': n_components,
        'variance_explained_percent': var_explained * 100,
        'cumulative_variance_percent': cum_var * 100,
        'component_amplitudes': np.max(np.abs(components), axis=1),
        'component_rms': np.sqrt(np.mean(components**2, axis=1)),
    }
    
    return stats


def filter_datacube_by_pca(datacube: np.ndarray,
                          n_components: int,
                          standardize: bool = True,
                          mask: Optional[np.ndarray] = None) -> np.ndarray:
    """
    Filter datacube using PCA - removes noise by keeping only top N components.
    
    Parameters
    ----------
    datacube : np.ndarray
        3D datacube array
    n_components : int
        Number of top components to keep
    standardize : bool, optional
        Whether to standardize before PCA
    mask : np.ndarray, optional
        Spatial mask for pixels to include
    
    Returns
    -------
    np.ndarray
        Filtered datacube (same shape as input)
    """
    results = perform_spectral_pca(datacube, n_components=n_components,
                                   standardize=standardize, mask=mask)
    filtered = reconstruct_from_pca(results, n_components=n_components)
    return filtered


__all__ = [
    'perform_spectral_pca',
    'reconstruct_from_pca',
    'get_component_stats',
    'filter_datacube_by_pca',
]
