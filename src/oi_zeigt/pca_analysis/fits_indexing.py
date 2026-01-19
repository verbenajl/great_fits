"""
FITS file indexing and metadata extraction.

Scans FITS files in a directory, extracts metadata from headers, and creates
an indexed DataFrame for filtering and accessing spectral data by mission,
telescope, scan, and source.
"""

import numpy as np
import logging
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Union
from astropy.io import fits
from astropy.table import Table
import pandas as pd

from .errors import FITSError, NoSpectraFoundError, InsufficientDataError

logger = logging.getLogger(__name__)


class FITSIndexer:
    """
    Index and manage FITS spectral files.
    
    Scans a directory for FITS files, extracts metadata from headers,
    and provides convenient access to spectral data with filtering
    by mission, telescope, scan, and source.
    """
    
    def __init__(self, fits_directory: Union[str, Path]):
        """
        Initialize FITS indexer.
        
        Parameters
        ----------
        fits_directory : str or Path
            Directory containing FITS files
        """
        self.fits_directory = Path(fits_directory)
        self.index_df = None
        self.fits_files = []
        
    def scan_directory(self) -> None:
        """
        Scan directory and index all FITS files.
        
        Raises
        ------
        FITSError
            If directory doesn't exist or contains no FITS files
        """
        if not self.fits_directory.exists():
            raise FITSError(f"FITS directory not found: {self.fits_directory}")
        
        # Find all FITS files
        fits_files = list(self.fits_directory.glob("*.fits"))
        fits_files += list(self.fits_directory.glob("*.fits.gz"))
        
        if not fits_files:
            raise FITSError(f"No FITS files found in {self.fits_directory}")
        
        self.fits_files = sorted(fits_files)
        logger.info(f"Found {len(self.fits_files)} FITS files")
        
        # Build index
        self._build_index()
    
    def _build_index(self) -> None:
        """Build pandas DataFrame index from FITS headers."""
        records = []
        
        for filepath in self.fits_files:
            try:
                with fits.open(filepath) as hdul:
                    # Get primary header
                    header = hdul[0].header
                    
                    # Extract standard keywords
                    record = self._extract_header_metadata(header, filepath)
                    records.append(record)
            
            except Exception as e:
                logger.warning(f"Failed to read FITS file {filepath}: {e}")
                continue
        
        if not records:
            raise FITSError("No valid FITS files could be indexed")
        
        self.index_df = pd.DataFrame(records)
        logger.info(f"Built index with {len(self.index_df)} entries")
    
    def _extract_header_metadata(self, header: fits.Header, filepath: Path) -> dict:
        """
        Extract metadata from FITS header.
        
        Parameters
        ----------
        header : fits.Header
            FITS header
        filepath : Path
            Path to FITS file
            
        Returns
        -------
        dict
            Dictionary with extracted metadata
        """
        record = {
            'filepath': str(filepath),
            'filename': filepath.name,
        }
        
        # Standard FITS keywords
        standard_keys = {
            'MISSION': 'mission',
            'TELESCOP': 'telescope',
            'INSTRUME': 'instrument',
            'OBSERVER': 'observer',
            'OBJECT': 'source',
            'SCAN': 'scan',
            'SUBSCAN': 'subscan',
            'DATE-OBS': 'date_obs',
            'TIME-OBS': 'time_obs',
            'NAXIS3': 'n_channels',
            'NAXIS1': 'n_subscan',
            'NAXIS2': 'n_readouts',
            'CRVAL3': 'freq_start',
            'CDELT3': 'freq_step',
        }
        
        for fits_key, df_key in standard_keys.items():
            if fits_key in header:
                record[df_key] = header[fits_key]
            else:
                record[df_key] = None
        
        # Additional useful keywords
        for key in ['EXTVER', 'EXTNAME', 'ORIGIN', 'DATATYPE']:
            if key in header:
                record[key.lower()] = header[key]
            else:
                record[key.lower()] = None
        
        return record
    
    def filter_by_mission(self, mission: str) -> 'FITSIndexer':
        """
        Filter index by mission.
        
        Parameters
        ----------
        mission : str
            Mission name (e.g., 'SOFIA', 'OI')
            
        Returns
        -------
        FITSIndexer
            New indexer with filtered data
        """
        if self.index_df is None:
            raise FITSError("Index not built. Call scan_directory() first.")
        
        filtered = self.index_df[self.index_df['mission'] == mission]
        
        if len(filtered) == 0:
            raise NoSpectraFoundError(f"No FITS files found for mission: {mission}")
        
        new_indexer = FITSIndexer(self.fits_directory)
        new_indexer.index_df = filtered.reset_index(drop=True)
        new_indexer.fits_files = [Path(f) for f in filtered['filepath']]
        
        return new_indexer
    
    def filter_by_telescope(self, telescope: str) -> 'FITSIndexer':
        """
        Filter index by telescope.
        
        Parameters
        ----------
        telescope : str
            Telescope name
            
        Returns
        -------
        FITSIndexer
            New indexer with filtered data
        """
        if self.index_df is None:
            raise FITSError("Index not built. Call scan_directory() first.")
        
        filtered = self.index_df[self.index_df['telescope'] == telescope]
        
        if len(filtered) == 0:
            raise NoSpectraFoundError(f"No FITS files found for telescope: {telescope}")
        
        new_indexer = FITSIndexer(self.fits_directory)
        new_indexer.index_df = filtered.reset_index(drop=True)
        new_indexer.fits_files = [Path(f) for f in filtered['filepath']]
        
        return new_indexer
    
    def filter_by_source(self, source: str) -> 'FITSIndexer':
        """
        Filter index by source.
        
        Parameters
        ----------
        source : str
            Source name (e.g., 'M51CENTER', 'SKYCHOPDIFF')
            
        Returns
        -------
        FITSIndexer
            New indexer with filtered data
        """
        if self.index_df is None:
            raise FITSError("Index not built. Call scan_directory() first.")
        
        filtered = self.index_df[self.index_df['source'] == source]
        
        if len(filtered) == 0:
            raise NoSpectraFoundError(f"No FITS files found for source: {source}")
        
        new_indexer = FITSIndexer(self.fits_directory)
        new_indexer.index_df = filtered.reset_index(drop=True)
        new_indexer.fits_files = [Path(f) for f in filtered['filepath']]
        
        return new_indexer
    
    def filter_by_scan(self, scan: int) -> 'FITSIndexer':
        """
        Filter index by scan number.
        
        Parameters
        ----------
        scan : int
            Scan number
            
        Returns
        -------
        FITSIndexer
            New indexer with filtered data
        """
        if self.index_df is None:
            raise FITSError("Index not built. Call scan_directory() first.")
        
        filtered = self.index_df[self.index_df['scan'] == scan]
        
        if len(filtered) == 0:
            raise NoSpectraFoundError(f"No FITS files found for scan: {scan}")
        
        new_indexer = FITSIndexer(self.fits_directory)
        new_indexer.index_df = filtered.reset_index(drop=True)
        new_indexer.fits_files = [Path(f) for f in filtered['filepath']]
        
        return new_indexer
    
    def get_unique_values(self, column: str) -> List:
        """
        Get unique values in a column.
        
        Parameters
        ----------
        column : str
            Column name
            
        Returns
        -------
        list
            Unique values in column
        """
        if self.index_df is None:
            raise FITSError("Index not built. Call scan_directory() first.")
        
        if column not in self.index_df.columns:
            raise ValueError(f"Column not found: {column}")
        
        return self.index_df[column].dropna().unique().tolist()
    
    def get_fits_data(self, index: int, hdu: int = 0) -> Tuple[np.ndarray, fits.Header]:
        """
        Get data and header from specific FITS file.
        
        Parameters
        ----------
        index : int
            Row index in the DataFrame
        hdu : int, optional
            HDU index (default: 0)
            
        Returns
        -------
        data : np.ndarray
            Data from FITS file
        header : fits.Header
            Header from FITS file
        """
        if self.index_df is None:
            raise FITSError("Index not built. Call scan_directory() first.")
        
        if index < 0 or index >= len(self.index_df):
            raise IndexError(f"Index out of bounds: {index}")
        
        filepath = self.index_df.iloc[index]['filepath']
        
        try:
            with fits.open(filepath) as hdul:
                data = hdul[hdu].data
                header = hdul[hdu].header
        except Exception as e:
            raise FITSError(f"Failed to read FITS data from {filepath}: {e}")
        
        return data, header
    
    def get_all_fits_data(self, hdu: int = 0) -> List[Tuple[np.ndarray, fits.Header, dict]]:
        """
        Get data from all FITS files in index.
        
        Parameters
        ----------
        hdu : int, optional
            HDU index (default: 0)
            
        Returns
        -------
        list of tuples
            Each tuple contains (data, header, metadata_row)
        """
        if self.index_df is None:
            raise FITSError("Index not built. Call scan_directory() first.")
        
        if len(self.index_df) == 0:
            raise NoSpectraFoundError("No FITS files in index")
        
        results = []
        
        for idx in range(len(self.index_df)):
            try:
                data, header = self.get_fits_data(idx, hdu=hdu)
                metadata = self.index_df.iloc[idx].to_dict()
                results.append((data, header, metadata))
            except Exception as e:
                logger.warning(f"Failed to read FITS file at index {idx}: {e}")
                continue
        
        if not results:
            raise FITSError("Failed to read any FITS files")
        
        return results
    
    def summary(self) -> str:
        """
        Get summary of indexed FITS files.
        
        Returns
        -------
        str
            Summary string
        """
        if self.index_df is None:
            return "Index not built"
        
        lines = [
            f"FITS Index Summary:",
            f"  Total files: {len(self.index_df)}",
            f"  Missions: {', '.join(map(str, self.get_unique_values('mission')))}",
            f"  Telescopes: {', '.join(map(str, self.get_unique_values('telescope')))}",
            f"  Sources: {', '.join(map(str, self.get_unique_values('source')))}",
            f"  Scans: {len(self.get_unique_values('scan'))} unique",
        ]
        
        return "\n".join(lines)
    
    def get_dataframe(self) -> pd.DataFrame:
        """
        Get the full index DataFrame.
        
        Returns
        -------
        pd.DataFrame
            Index DataFrame
        """
        if self.index_df is None:
            raise FITSError("Index not built. Call scan_directory() first.")
        
        return self.index_df.copy()


def load_spectral_data(filepath: Union[str, Path], hdu: int = 0) -> Tuple[np.ndarray, fits.Header]:
    """
    Load spectral data from FITS file.
    
    Parameters
    ----------
    filepath : str or Path
        Path to FITS file
    hdu : int, optional
        HDU index (default: 0)
        
    Returns
    -------
    data : np.ndarray
        Spectral data
    header : fits.Header
        FITS header
        
    Raises
    ------
    FITSError
        If file cannot be read
    """
    try:
        filepath = Path(filepath)
        if not filepath.exists():
            raise FileNotFoundError(f"FITS file not found: {filepath}")
        
        with fits.open(filepath) as hdul:
            if hdu >= len(hdul):
                raise IndexError(f"HDU {hdu} not found in {filepath}")
            
            data = hdul[hdu].data
            header = hdul[hdu].header
        
        if data is None:
            raise FITSError(f"No data in HDU {hdu} of {filepath}")
        
        return data, header
    
    except Exception as e:
        raise FITSError(f"Failed to load spectral data from {filepath}: {e}")


def save_spectral_data(data: np.ndarray, filepath: Union[str, Path],
                       header: Optional[fits.Header] = None,
                       overwrite: bool = False) -> None:
    """
    Save spectral data to FITS file.
    
    Parameters
    ----------
    data : np.ndarray
        Spectral data
    filepath : str or Path
        Output FITS file path
    header : fits.Header, optional
        FITS header
    overwrite : bool, optional
        Overwrite existing file
    """
    try:
        filepath = Path(filepath)
        filepath.parent.mkdir(parents=True, exist_ok=True)
        
        if header is None:
            header = fits.Header()
        
        hdu = fits.PrimaryHDU(data=data, header=header)
        hdu.writeto(filepath, overwrite=overwrite)
        
        logger.info(f"Saved spectral data to {filepath}")
    
    except Exception as e:
        raise FITSError(f"Failed to save spectral data to {filepath}: {e}")
