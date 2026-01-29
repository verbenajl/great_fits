#!/usr/bin/env python3
"""
pytest tests for NaN and blank value detection.
"""
import sys
from pathlib import Path
import numpy as np
import pytest
from astropy.io import fits

# Add src to path
sys.path.insert(0, str(Path(__file__).parent / 'src'))

from oi_zeigt.reduction.core import detect_nan_channels, detect_missing_channels


class TestDetectNanChannels:
    """Tests for detect_nan_channels function."""
    
    def test_nan_detection_basic(self):
        """Test basic NaN detection."""
        spectrum = np.array([1.0, 2.0, np.nan, 4.0, np.nan, 6.0])
        mask, frac = detect_nan_channels(spectrum)
        
        assert np.sum(mask) == 2
        assert np.isclose(frac, 2/6)
    
    def test_no_nans(self):
        """Test spectrum with no NaNs."""
        spectrum = np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0])
        mask, frac = detect_nan_channels(spectrum)
        
        assert np.sum(mask) == 0
        assert frac == 0.0
    
    def test_all_nans(self):
        """Test spectrum with all NaNs."""
        spectrum = np.array([np.nan, np.nan, np.nan])
        mask, frac = detect_nan_channels(spectrum)
        
        assert np.sum(mask) == 3
        assert frac == 1.0
    
    def test_single_nan(self):
        """Test spectrum with single NaN."""
        spectrum = np.array([1.0, np.nan, 3.0])
        mask, frac = detect_nan_channels(spectrum)
        
        assert np.sum(mask) == 1
        assert np.isclose(frac, 1/3)


class TestDetectMissingChannels:
    """Tests for detect_missing_channels function."""
    
    def test_nan_only_mode(self):
        """Test NaN-only detection (minimal mode)."""
        spectrum = np.array([1.0, 0.0, np.nan, 4.0])
        mask, frac = detect_missing_channels(spectrum, include_blanks=False, gildas_blank=False)
        
        # Should only detect NaN (1 out of 4)
        assert np.sum(mask) == 1
        assert np.isclose(frac, 0.25)
    
    def test_with_blank_zeros(self):
        """Test detection including zero blanks."""
        spectrum = np.array([1.0, 0.0, np.nan, 4.0])
        mask, frac = detect_missing_channels(spectrum, include_blanks=True, gildas_blank=False)
        
        # Should detect NaN + 0.0 (2 out of 4)
        assert np.sum(mask) == 2
        assert np.isclose(frac, 0.5)
    
    def test_gildas_blanks(self):
        """Test detection of GILDAS-style blanks."""
        spectrum = np.array([1.0, 2.0, -1e31, 4.0, np.nan])
        mask, frac = detect_missing_channels(spectrum, include_blanks=False, gildas_blank=True)
        
        # Should detect NaN + GILDAS blank (2 out of 5)
        assert np.sum(mask) == 2
        assert np.isclose(frac, 0.4)
    
    def test_infinity_detection(self):
        """Test detection of infinity values."""
        spectrum = np.array([1.0, np.inf, np.nan, 4.0, -np.inf])
        mask, frac = detect_missing_channels(spectrum, include_blanks=False, gildas_blank=False)
        
        # Should detect NaN + inf + -inf (3 out of 5)
        assert np.sum(mask) == 3
        assert np.isclose(frac, 0.6)
    
    def test_comprehensive_detection(self):
        """Test full comprehensive detection."""
        spectrum = np.array([100.0, 0.0, np.nan, -1e31, np.inf, 200.0, -np.inf, 300.0])
        mask, frac = detect_missing_channels(spectrum, include_blanks=True, gildas_blank=True)
        
        # Should detect: 0.0, NaN, -1e31, inf, -inf (5 out of 8)
        assert np.sum(mask) == 5
        assert np.isclose(frac, 5/8)
    
    def test_custom_blank_value(self):
        """Test detection with custom blank value."""
        spectrum = np.array([1.0, -999.0, np.nan, 4.0])
        mask, frac = detect_missing_channels(spectrum, 
                                             include_blanks=True, 
                                             blank_value=-999.0,
                                             gildas_blank=False)
        
        # Should detect NaN + -999.0 (2 out of 4)
        assert np.sum(mask) == 2
        assert np.isclose(frac, 0.5)


class TestRealFITSData:
    """Tests using real FITS data if available."""
    
    @pytest.fixture
    def fits_file(self):
        """Fixture providing path to test FITS file."""
        path = Path('/home/diskB/sofiaobsdata/m51/m51_central_region/kalibrated_all/m51_central_tile.fits')
        return path if path.exists() else None
    
    def test_real_data_nan_detection(self, fits_file):
        """Test NaN detection on real FITS data."""
        if fits_file is None:
            pytest.skip("Test FITS file not found")
        
        with fits.open(fits_file) as hdul:
            for hdu in hdul:
                if hasattr(hdu, 'data') and hdu.data is not None:
                    if 'SPECTRUM' in hdu.data.dtype.names:
                        # Test first spectrum
                        spectrum = hdu.data[0]['SPECTRUM']
                        mask, frac = detect_nan_channels(spectrum)
                        
                        # Real data should have some NaNs
                        assert frac > 0.0
                        assert frac < 0.5  # But not too many
                        assert np.sum(mask) > 0
                        return
        
        pytest.fail("No SPECTRUM column found in FITS file")
    
    def test_real_data_no_gildas_blanks(self, fits_file):
        """Verify real data doesn't have GILDAS blanks."""
        if fits_file is None:
            pytest.skip("Test FITS file not found")
        
        with fits.open(fits_file) as hdul:
            for hdu in hdul:
                if hasattr(hdu, 'data') and hdu.data is not None:
                    if 'SPECTRUM' in hdu.data.dtype.names:
                        # Test first few spectra
                        for i in range(min(10, len(hdu.data))):
                            spectrum = hdu.data[i]['SPECTRUM']
                            
                            # Check for GILDAS blanks
                            gildas_mask = np.abs(spectrum) > 1e30
                            assert np.sum(gildas_mask) == 0, f"Unexpected GILDAS blanks in row {i}"
                        return
        
        pytest.fail("No SPECTRUM column found in FITS file")


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
