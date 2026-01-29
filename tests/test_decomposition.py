"""
Unit tests for Phase 2 Decomposition module.

Tests the complete PCA decomposition workflow:
- DecompositionConfig configuration
- PCADecomposer sklearn integration
- DecompositionResult results storage
"""

import pytest
import tempfile
import numpy as np
from pathlib import Path
import pickle

from src.oi_zeigt.pca_analysis import (
    DecompositionConfig,
    PCADecomposer,
    DecompositionResult,
)
from src.oi_zeigt.pca_analysis.errors import DecompositionError


class TestDecompositionConfig:
    """Tests for DecompositionConfig dataclass."""
    
    def test_init_required(self):
        """Test with required n_components."""
        config = DecompositionConfig(n_components=5)
        
        assert config.n_components == 5
        assert config.pca_source == "SKYCHOPDIFF"
        assert config.normalize == True
        assert config.scale == False
    
    def test_init_custom_all(self):
        """Test with custom parameters."""
        config = DecompositionConfig(
            n_components=10,
            pca_source="CUSTOM",
            normalize=False,
            scale=True,
            add_sky_diff=True,
            noise_cutoff=True,
            scramble=True,
        )
        
        assert config.n_components == 10
        assert config.pca_source == "CUSTOM"
        assert config.normalize == False
        assert config.scale == True
        assert config.add_sky_diff == True
        assert config.noise_cutoff == True
        assert config.scramble == True
    
    def test_to_dict(self):
        """Test conversion to dictionary."""
        config = DecompositionConfig(n_components=8, pca_source="TEST")
        config_dict = config.to_dict()
        
        assert isinstance(config_dict, dict)
        assert config_dict['n_components'] == 8
        assert config_dict['pca_source'] == 'TEST'
        assert 'normalize' in config_dict
        assert 'scale' in config_dict
    
    def test_from_dict(self):
        """Test creation from dictionary."""
        config_dict = {
            'n_components': 7,
            'pca_source': 'TEST',
            'normalize': False,
            'scale': True,
        }
        config = DecompositionConfig.from_dict(config_dict)
        
        assert config.n_components == 7
        assert config.pca_source == 'TEST'
        assert config.normalize == False
        assert config.scale == True
    
    def test_roundtrip_serialization(self):
        """Test to_dict -> from_dict roundtrip."""
        original = DecompositionConfig(
            n_components=12,
            pca_source="MISSION",
            normalize=False,
            scale=True,
        )
        
        config_dict = original.to_dict()
        restored = DecompositionConfig.from_dict(config_dict)
        
        assert restored.n_components == original.n_components
        assert restored.pca_source == original.pca_source
        assert restored.normalize == original.normalize
        assert restored.scale == original.scale


class TestPCADecomposer:
    """Tests for PCADecomposer class."""
    
    def setup_method(self):
        """Set up test fixtures."""
        self.decomposer = PCADecomposer(n_components=3)
    
    def test_init(self):
        """Test PCADecomposer initialization."""
        assert self.decomposer.n_components == 3
        assert self.decomposer.fitted == False
        assert self.decomposer.pca_model is None
    
    def test_fit_valid_data(self):
        """Test fitting PCA to valid spectral data."""
        # Create synthetic spectral matrix [n_spectra, n_channels]
        n_spectra = 50
        n_channels = 256
        spectral_matrix = np.random.randn(n_spectra, n_channels)
        
        result = self.decomposer.fit(spectral_matrix)
        
        # fit() returns self
        assert result is self.decomposer
        assert self.decomposer.fitted == True
        assert self.decomposer.pca_model is not None
    
    def test_get_components(self):
        """Test retrieving PCA components."""
        spectral_matrix = np.random.randn(30, 128)
        self.decomposer.fit(spectral_matrix)
        
        components = self.decomposer.get_components()
        
        assert components.shape == (3, 128)
        assert isinstance(components, np.ndarray)
    
    def test_get_explained_variance_ratio(self):
        """Test getting explained variance ratios."""
        spectral_matrix = np.random.randn(40, 256)
        self.decomposer.fit(spectral_matrix)
        
        variance = self.decomposer.get_explained_variance_ratio()
        
        assert variance.shape == (3,)
        assert np.all(variance > 0)
        assert np.all(variance <= 1.0)
    
    def test_get_explained_variance(self):
        """Test getting absolute explained variance."""
        spectral_matrix = np.random.randn(40, 256)
        self.decomposer.fit(spectral_matrix)
        
        variance = self.decomposer.get_explained_variance()
        
        assert variance.shape == (3,)
        assert np.all(variance > 0)
    
    def test_transform(self):
        """Test transforming data to component space."""
        spectral_matrix = np.random.randn(50, 128)
        self.decomposer.fit(spectral_matrix)
        
        transformed = self.decomposer.transform(spectral_matrix)
        
        assert transformed.shape == (50, 3)
        assert isinstance(transformed, np.ndarray)
    
    def test_reconstruct(self):
        """Test reconstructing spectra from coefficients."""
        spectral_matrix = np.random.randn(30, 256)
        self.decomposer.fit(spectral_matrix)
        
        coefficients = np.random.randn(10, 3)
        reconstructed = self.decomposer.reconstruct(coefficients)
        
        assert reconstructed.shape == (10, 256)
        assert isinstance(reconstructed, np.ndarray)
    
    def test_roundtrip_transform_reconstruct(self):
        """Test that transform->reconstruct gives similar data."""
        spectral_matrix = np.random.randn(20, 128)
        self.decomposer.fit(spectral_matrix)
        
        # Transform first few spectra
        subset = spectral_matrix[:5]
        coefficients = self.decomposer.transform(subset)
        reconstructed = self.decomposer.reconstruct(coefficients)
        
        # Reconstruction should be reasonably close
        mse = np.mean((subset - reconstructed) ** 2)
        assert mse < 2.0  # Allow some reconstruction error
    
    def test_fit_error_unfitted(self):
        """Test error when calling methods on unfitted decomposer."""
        with pytest.raises(DecompositionError):
            self.decomposer.get_components()
    
    def test_multiple_fits(self):
        """Test refitting with different data."""
        # First fit
        data1 = np.random.randn(30, 100)
        self.decomposer.fit(data1)
        comp1 = self.decomposer.get_components()
        
        # Second fit with different data
        data2 = np.random.randn(50, 100)
        self.decomposer.fit(data2)
        comp2 = self.decomposer.get_components()
        
        # Components should be different
        assert not np.allclose(comp1, comp2)


class TestDecompositionResult:
    """Tests for DecompositionResult class."""
    
    def setup_method(self):
        """Set up test fixtures."""
        self.components = np.random.randn(3, 256)
        self.explained_variance_ratio = np.array([0.5, 0.3, 0.2])
        self.explained_variance = np.array([100.0, 60.0, 40.0])
        self.mean_spectrum = np.random.randn(256)
        self.config = DecompositionConfig(n_components=3)
        self.metadata = {'mission': 'TEST', 'tag': 'test_run'}
    
    def test_init_minimal(self):
        """Test initialization with required parameters."""
        result = DecompositionResult(
            components=self.components,
            explained_variance_ratio=self.explained_variance_ratio,
            explained_variance=self.explained_variance,
            mean_spectrum=self.mean_spectrum,
            config=self.config,
            metadata=self.metadata,
        )
        
        assert np.array_equal(result.components, self.components)
        assert np.array_equal(result.explained_variance_ratio, self.explained_variance_ratio)
        assert np.array_equal(result.mean_spectrum, self.mean_spectrum)
        assert result.config == self.config
    
    def test_spectrum_metadata_default(self):
        """Test that spectrum_metadata defaults to empty list."""
        result = DecompositionResult(
            components=self.components,
            explained_variance_ratio=self.explained_variance_ratio,
            explained_variance=self.explained_variance,
            mean_spectrum=self.mean_spectrum,
            config=self.config,
            metadata=self.metadata,
        )
        
        assert isinstance(result.spectrum_metadata, list)
        assert len(result.spectrum_metadata) == 0
    
    def test_summary(self):
        """Test summary generation."""
        result = DecompositionResult(
            components=self.components,
            explained_variance_ratio=self.explained_variance_ratio,
            explained_variance=self.explained_variance,
            mean_spectrum=self.mean_spectrum,
            config=self.config,
            metadata={'mission': 'TEST', 'n_spectra': 100},
        )
        
        summary = result.summary()
        
        assert isinstance(summary, str)
        assert len(summary) > 0
        # Check for key terms
        assert any(term in summary.lower() for term in ['component', 'variance', 'summary'])
    
    def test_save_and_load(self):
        """Test saving and loading results."""
        result = DecompositionResult(
            components=self.components,
            explained_variance_ratio=self.explained_variance_ratio,
            explained_variance=self.explained_variance,
            mean_spectrum=self.mean_spectrum,
            config=self.config,
            metadata={'mission': 'TESTMISSION', 'n_spectra': 50},
        )
        
        with tempfile.NamedTemporaryFile(suffix='.pkl', delete=False) as f:
            filepath = f.name
        
        try:
            result.save(filepath)
            
            # Verify file exists
            assert Path(filepath).exists()
            
            # Load and verify
            loaded = DecompositionResult.load(filepath)
            
            assert np.allclose(loaded.components, self.components)
            assert np.allclose(loaded.explained_variance_ratio, self.explained_variance_ratio)
            assert np.allclose(loaded.mean_spectrum, self.mean_spectrum)
        finally:
            Path(filepath).unlink(missing_ok=True)
    
    def test_metadata_preservation(self):
        """Test that metadata is preserved through save/load."""
        metadata = {
            'mission': 'MISSION123',
            'tag': 'test_tag',
            'n_spectra': 150,
            'custom_field': 'custom_value',
            'timestamp': '2024-01-15T10:30:00',
        }
        
        result = DecompositionResult(
            components=self.components,
            explained_variance_ratio=self.explained_variance_ratio,
            explained_variance=self.explained_variance,
            mean_spectrum=self.mean_spectrum,
            config=self.config,
            metadata=metadata,
        )
        
        with tempfile.NamedTemporaryFile(suffix='.pkl', delete=False) as f:
            filepath = f.name
        
        try:
            result.save(filepath)
            loaded = DecompositionResult.load(filepath)
            
            assert loaded.metadata['mission'] == 'MISSION123'
            assert loaded.metadata['n_spectra'] == 150
            assert loaded.metadata['custom_field'] == 'custom_value'
            assert loaded.metadata['timestamp'] == '2024-01-15T10:30:00'
        finally:
            Path(filepath).unlink(missing_ok=True)
    
    def test_spectrum_metadata_preservation(self):
        """Test that spectrum metadata list is preserved."""
        spectrum_metadata = [
            {'spectrum_id': 0, 'quality': 'good'},
            {'spectrum_id': 1, 'quality': 'excellent'},
            {'spectrum_id': 2, 'quality': 'fair'},
        ]
        
        result = DecompositionResult(
            components=self.components,
            explained_variance_ratio=self.explained_variance_ratio,
            explained_variance=self.explained_variance,
            mean_spectrum=self.mean_spectrum,
            config=self.config,
            metadata=self.metadata,
            spectrum_metadata=spectrum_metadata,
        )
        
        with tempfile.NamedTemporaryFile(suffix='.pkl', delete=False) as f:
            filepath = f.name
        
        try:
            result.save(filepath)
            loaded = DecompositionResult.load(filepath)
            
            assert len(loaded.spectrum_metadata) == 3
            assert loaded.spectrum_metadata[0]['spectrum_id'] == 0
            assert loaded.spectrum_metadata[1]['quality'] == 'excellent'
        finally:
            Path(filepath).unlink(missing_ok=True)


class TestIntegration:
    """Integration tests for complete decomposition workflow."""
    
    def test_decomposer_with_synthetic_data(self):
        """Test complete decomposition with synthetic data."""
        config = DecompositionConfig(n_components=4)
        
        # Create synthetic data (more samples than components)
        spectral_matrix = np.random.randn(100, 256)
        
        # Decompose
        decomposer = PCADecomposer(n_components=4)
        decomposer.fit(spectral_matrix)
        components = decomposer.get_components()
        variance_ratio = decomposer.get_explained_variance_ratio()
        
        assert components.shape == (4, 256)
        assert variance_ratio.shape == (4,)
        # With only 4 components from 256 dimensions, ratios will be < 1.0 total
        assert np.all(variance_ratio > 0)
    
    def test_result_roundtrip_complete(self):
        """Test complete save/load cycle with all attributes."""
        config = DecompositionConfig(n_components=5)
        components = np.random.randn(5, 128)
        variance_ratio = np.array([0.4, 0.3, 0.15, 0.1, 0.05])
        variance = np.array([100.0, 75.0, 37.5, 25.0, 12.5])
        mean_spectrum = np.random.randn(128)
        
        spectrum_metadata = [
            {'id': i, 'quality': 'good'} for i in range(50)
        ]
        
        original_result = DecompositionResult(
            components=components,
            explained_variance_ratio=variance_ratio,
            explained_variance=variance,
            mean_spectrum=mean_spectrum,
            config=config,
            metadata={
                'mission': 'TESTMISSION',
                'tag': 'integration_test',
                'n_spectra': 100,
            },
            spectrum_metadata=spectrum_metadata,
        )
        
        with tempfile.NamedTemporaryFile(suffix='.pkl', delete=False) as f:
            filepath = f.name
        
        try:
            original_result.save(filepath)
            loaded_result = DecompositionResult.load(filepath)
            
            # Verify all attributes preserved
            assert np.allclose(loaded_result.components, original_result.components)
            assert np.allclose(loaded_result.explained_variance_ratio, 
                             original_result.explained_variance_ratio)
            assert np.allclose(loaded_result.explained_variance, 
                             original_result.explained_variance)
            assert np.allclose(loaded_result.mean_spectrum, original_result.mean_spectrum)
            assert loaded_result.metadata['mission'] == 'TESTMISSION'
            assert len(loaded_result.spectrum_metadata) == 50
        finally:
            Path(filepath).unlink(missing_ok=True)
    
    def test_pca_with_different_components(self):
        """Test PCA with different numbers of components."""
        spectral_matrix = np.random.randn(100, 256)
        
        for n_comp in [1, 3, 5, 10]:
            decomposer = PCADecomposer(n_components=n_comp)
            decomposer.fit(spectral_matrix)
            
            components = decomposer.get_components()
            assert components.shape == (n_comp, 256)
            
            variance = decomposer.get_explained_variance_ratio()
            assert variance.shape == (n_comp,)


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
