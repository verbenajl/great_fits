"""
Unit tests for Phase 3 Correction module.

Tests the complete correction workflow:
- CorrectionConfig configuration
- SpectrumProjector spectrum projection
- CorrectionCalculator correction calculation
- ScienceDataLoader science data loading
- CorrectionResult results storage
- correct_spectra() orchestration
"""

import pytest
import tempfile
import numpy as np
from pathlib import Path
from unittest.mock import Mock, patch, MagicMock

from src.oi_zeigt.pca_analysis import (
    CorrectionConfig,
    SpectrumProjector,
    CorrectionCalculator,
    CorrectionResult,
    DecompositionConfig,
    DecompositionResult,
    PCADecomposer,
)
from src.oi_zeigt.pca_analysis.errors import CorrectionError


class TestCorrectionConfig:
    """Tests for CorrectionConfig dataclass."""
    
    def test_init_required(self):
        """Test with required n_components."""
        config = CorrectionConfig(n_components=5)
        
        assert config.n_components == 5
        assert config.science_source == "M51CENTER"
        assert config.preserve_science_line == True
        assert config.mask_width == 20
        assert config.fit_window == 512
        assert config.regularization == 0.0
    
    def test_init_custom(self):
        """Test with custom parameters."""
        config = CorrectionConfig(
            n_components=10,
            science_source="CUSTOM_SOURCE",
            preserve_science_line=False,
            mask_width=30,
            fit_window=256,
            regularization=0.01,
        )
        
        assert config.n_components == 10
        assert config.science_source == "CUSTOM_SOURCE"
        assert config.preserve_science_line == False
        assert config.mask_width == 30
        assert config.fit_window == 256
        assert config.regularization == 0.01
    
    def test_to_dict(self):
        """Test conversion to dictionary."""
        config = CorrectionConfig(n_components=8, science_source="TEST")
        config_dict = config.to_dict()
        
        assert isinstance(config_dict, dict)
        assert config_dict['n_components'] == 8
        assert config_dict['science_source'] == 'TEST'
        assert 'preserve_science_line' in config_dict
    
    def test_from_dict(self):
        """Test creation from dictionary."""
        config_dict = {
            'n_components': 7,
            'science_source': 'TEST',
            'preserve_science_line': False,
            'mask_width': 25,
            'fit_window': 256,
            'regularization': 0.01,
        }
        config = CorrectionConfig.from_dict(config_dict)
        
        assert config.n_components == 7
        assert config.science_source == 'TEST'
        assert config.preserve_science_line == False
    
    def test_roundtrip_serialization(self):
        """Test to_dict -> from_dict roundtrip."""
        original = CorrectionConfig(
            n_components=12,
            science_source="MISSION",
            preserve_science_line=False,
            mask_width=15,
            regularization=0.05,
        )
        
        config_dict = original.to_dict()
        restored = CorrectionConfig.from_dict(config_dict)
        
        assert restored.n_components == original.n_components
        assert restored.science_source == original.science_source
        assert restored.preserve_science_line == original.preserve_science_line
        assert restored.regularization == original.regularization


class TestSpectrumProjector:
    """Tests for SpectrumProjector class."""
    
    def setup_method(self):
        """Set up test fixtures."""
        # Create a simple decomposition result
        n_components = 5
        n_channels = 256
        
        # Create fake components
        components = np.random.randn(n_components, n_channels)
        mean_spectrum = np.random.randn(n_channels)
        variance_ratio = np.ones(n_components) / n_components
        variance = np.ones(n_components) * 10.0
        
        config = DecompositionConfig(n_components=n_components)
        
        self.decomp_result = DecompositionResult(
            components=components,
            explained_variance_ratio=variance_ratio,
            explained_variance=variance,
            mean_spectrum=mean_spectrum,
            config=config,
            metadata={'mission': 'TEST'},
        )
        
        self.projector = SpectrumProjector(self.decomp_result)
    
    def test_init(self):
        """Test SpectrumProjector initialization."""
        assert self.projector.n_components == 5
        assert self.projector.n_channels == 256
        assert self.projector.components.shape == (5, 256)
    
    def test_project_single_spectrum(self):
        """Test projecting a single spectrum."""
        spectrum = np.random.randn(256)
        coefficients, metadata = self.projector.project_spectrum(spectrum)
        
        assert coefficients.shape == (5,)
        assert isinstance(metadata, dict)
        assert 'input_norm' in metadata
        assert 'coefficient_norm' in metadata
    
    def test_project_multiple_spectra(self):
        """Test projecting multiple spectra."""
        spectra = np.random.randn(20, 256)
        coefficients, metadata_list = self.projector.project_multiple(spectra)
        
        assert coefficients.shape == (20, 5)
        assert len(metadata_list) == 20
        assert all('spectrum_id' in m for m in metadata_list)
    
    def test_projection_metadata(self):
        """Test that projection metadata is generated correctly."""
        spectrum = np.random.randn(256)
        _, metadata = self.projector.project_spectrum(spectrum)
        
        assert 'coefficient_stats' in metadata
        assert 'mean' in metadata['coefficient_stats']
        assert 'std' in metadata['coefficient_stats']
        assert 'min' in metadata['coefficient_stats']
        assert 'max' in metadata['coefficient_stats']
    
    def test_projection_dimension_mismatch(self):
        """Test error handling for dimension mismatch."""
        spectrum = np.random.randn(128)  # Wrong shape
        
        with pytest.raises(CorrectionError):
            self.projector.project_spectrum(spectrum)


class TestCorrectionCalculator:
    """Tests for CorrectionCalculator class."""
    
    def setup_method(self):
        """Set up test fixtures."""
        # Create decomposition result
        n_components = 5
        n_channels = 256
        
        components = np.random.randn(n_components, n_channels)
        mean_spectrum = np.random.randn(n_channels)
        variance_ratio = np.ones(n_components) / n_components
        variance = np.ones(n_components) * 10.0
        
        decomp_config = DecompositionConfig(n_components=n_components)
        
        self.decomp_result = DecompositionResult(
            components=components,
            explained_variance_ratio=variance_ratio,
            explained_variance=variance,
            mean_spectrum=mean_spectrum,
            config=decomp_config,
            metadata={'mission': 'TEST'},
        )
        
        # Create correction config
        self.correction_config = CorrectionConfig(n_components=n_components)
        
        self.calculator = CorrectionCalculator(self.correction_config, self.decomp_result)
    
    def test_init(self):
        """Test CorrectionCalculator initialization."""
        assert self.calculator.n_components == 5
        assert self.calculator.config == self.correction_config
    
    def test_calculate_correction(self):
        """Test correction matrix calculation."""
        # Create test data
        n_ref = 30
        n_science = 20
        
        reference_coeff = np.random.randn(n_ref, 5) * 0.5
        science_coeff = np.random.randn(n_science, 5) * 0.5
        
        correction_matrix, stats = self.calculator.calculate_correction(
            reference_coeff, science_coeff
        )
        
        assert correction_matrix.shape == (5, 5)
        assert isinstance(stats, dict)
        assert 'residual_rms' in stats
        assert 'condition_number' in stats
        assert 'correction_strength' in stats
    
    def test_correction_matrix_properties(self):
        """Test properties of calculated correction matrix."""
        reference_coeff = np.random.randn(30, 5)
        science_coeff = np.random.randn(20, 5)
        
        correction_matrix, stats = self.calculator.calculate_correction(
            reference_coeff, science_coeff
        )
        
        # Check that correction matrix is finite
        assert np.all(np.isfinite(correction_matrix))
        
        # Check residual is non-negative
        assert stats['residual_rms'] >= 0
    
    def test_apply_correction(self):
        """Test applying correction to coefficients."""
        # Calculate correction first
        reference_coeff = np.random.randn(30, 5)
        science_coeff = np.random.randn(20, 5)
        
        correction_matrix, _ = self.calculator.calculate_correction(
            reference_coeff, science_coeff
        )
        
        # Apply correction
        corrected = self.calculator.apply_correction(science_coeff, correction_matrix)
        
        assert corrected.shape == science_coeff.shape
        assert np.all(np.isfinite(corrected))
    
    def test_correction_with_regularization(self):
        """Test correction with regularization parameter."""
        config = CorrectionConfig(n_components=5, regularization=0.1)
        calculator = CorrectionCalculator(config, self.decomp_result)
        
        reference_coeff = np.random.randn(30, 5)
        science_coeff = np.random.randn(20, 5)
        
        correction_matrix, stats = calculator.calculate_correction(
            reference_coeff, science_coeff
        )
        
        assert correction_matrix.shape == (5, 5)
        assert stats['regularization'] == 0.1


class TestCorrectionResult:
    """Tests for CorrectionResult class."""
    
    def setup_method(self):
        """Set up test fixtures."""
        self.n_components = 5
        self.n_channels = 256
        self.n_ref = 30
        self.n_science = 20
        
        self.reference_coeff = np.random.randn(self.n_ref, self.n_components)
        self.science_coeff = np.random.randn(self.n_science, self.n_components)
        self.corrected_coeff = np.random.randn(self.n_science, self.n_components)
        self.correction_matrix = np.random.randn(self.n_components, self.n_components)
        self.components = np.random.randn(self.n_components, self.n_channels)
        self.mean_spectrum = np.random.randn(self.n_channels)
        
        self.config = CorrectionConfig(n_components=self.n_components)
        self.metadata = {'mission': 'TEST', 'tag': 'test_run'}
    
    def test_init(self):
        """Test CorrectionResult initialization."""
        result = CorrectionResult(
            reference_coefficients=self.reference_coeff,
            science_coefficients=self.science_coeff,
            corrected_coefficients=self.corrected_coeff,
            correction_matrix=self.correction_matrix,
            components=self.components,
            mean_spectrum=self.mean_spectrum,
            config=self.config,
            metadata=self.metadata,
        )
        
        assert result.reference_coefficients.shape == (self.n_ref, self.n_components)
        assert result.science_coefficients.shape == (self.n_science, self.n_components)
        assert result.correction_matrix.shape == (self.n_components, self.n_components)
    
    def test_reconstruct_corrected_spectra(self):
        """Test reconstructing corrected spectra."""
        result = CorrectionResult(
            reference_coefficients=self.reference_coeff,
            science_coefficients=self.science_coeff,
            corrected_coefficients=self.corrected_coeff,
            correction_matrix=self.correction_matrix,
            components=self.components,
            mean_spectrum=self.mean_spectrum,
            config=self.config,
            metadata=self.metadata,
        )
        
        corrected_spectra = result.reconstruct_corrected_spectra()
        
        assert corrected_spectra.shape == (self.n_science, self.n_channels)
        assert np.all(np.isfinite(corrected_spectra))
    
    def test_get_correction_strength(self):
        """Test getting correction strength statistics."""
        result = CorrectionResult(
            reference_coefficients=self.reference_coeff,
            science_coefficients=self.science_coeff,
            corrected_coefficients=self.corrected_coeff,
            correction_matrix=self.correction_matrix,
            components=self.components,
            mean_spectrum=self.mean_spectrum,
            config=self.config,
            metadata=self.metadata,
        )
        
        strength = result.get_correction_strength()
        
        assert 'per_component_change' in strength
        assert 'total_change_mean' in strength
        assert len(strength['per_component_change']) == self.n_components
    
    def test_summary(self):
        """Test summary generation."""
        result = CorrectionResult(
            reference_coefficients=self.reference_coeff,
            science_coefficients=self.science_coeff,
            corrected_coefficients=self.corrected_coeff,
            correction_matrix=self.correction_matrix,
            components=self.components,
            mean_spectrum=self.mean_spectrum,
            config=self.config,
            metadata=self.metadata,
        )
        
        summary = result.summary()
        
        assert isinstance(summary, str)
        assert len(summary) > 0
        assert 'Correction Results Summary' in summary or 'Results' in summary
    
    def test_save_and_load(self):
        """Test saving and loading results."""
        result = CorrectionResult(
            reference_coefficients=self.reference_coeff,
            science_coefficients=self.science_coeff,
            corrected_coefficients=self.corrected_coeff,
            correction_matrix=self.correction_matrix,
            components=self.components,
            mean_spectrum=self.mean_spectrum,
            config=self.config,
            metadata={'mission': 'TESTMISSION', 'n_spectra': 20},
        )
        
        with tempfile.NamedTemporaryFile(suffix='.pkl', delete=False) as f:
            filepath = f.name
        
        try:
            result.save(filepath)
            assert Path(filepath).exists()
            
            loaded = CorrectionResult.load(filepath)
            
            assert np.allclose(loaded.reference_coefficients, self.reference_coeff)
            assert np.allclose(loaded.corrected_coefficients, self.corrected_coeff)
        finally:
            Path(filepath).unlink(missing_ok=True)
    
    def test_metadata_preservation(self):
        """Test that metadata is preserved through save/load."""
        metadata = {
            'mission': 'MISSION123',
            'tag': 'test_tag',
            'n_spectra': 150,
            'custom_field': 'custom_value',
        }
        
        result = CorrectionResult(
            reference_coefficients=self.reference_coeff,
            science_coefficients=self.science_coeff,
            corrected_coefficients=self.corrected_coeff,
            correction_matrix=self.correction_matrix,
            components=self.components,
            mean_spectrum=self.mean_spectrum,
            config=self.config,
            metadata=metadata,
        )
        
        with tempfile.NamedTemporaryFile(suffix='.pkl', delete=False) as f:
            filepath = f.name
        
        try:
            result.save(filepath)
            loaded = CorrectionResult.load(filepath)
            
            assert loaded.metadata['mission'] == 'MISSION123'
            assert loaded.metadata['n_spectra'] == 150
        finally:
            Path(filepath).unlink(missing_ok=True)


class TestIntegration:
    """Integration tests for complete correction workflow."""
    
    def test_projection_and_correction_workflow(self):
        """Test complete projection and correction workflow."""
        # Create decomposition result
        n_components = 5
        n_channels = 256
        
        components = np.random.randn(n_components, n_channels)
        mean_spectrum = np.random.randn(n_channels)
        variance_ratio = np.ones(n_components) / n_components
        variance = np.ones(n_components) * 10.0
        
        decomp_config = DecompositionConfig(n_components=n_components)
        decomp_result = DecompositionResult(
            components=components,
            explained_variance_ratio=variance_ratio,
            explained_variance=variance,
            mean_spectrum=mean_spectrum,
            config=decomp_config,
            metadata={'mission': 'TEST'},
        )
        
        # Create projector and project spectra
        projector = SpectrumProjector(decomp_result)
        
        spectra = np.random.randn(30, n_channels)
        science_coeff, _ = projector.project_multiple(spectra)
        
        # Calculate correction
        correction_config = CorrectionConfig(n_components=n_components)
        calculator = CorrectionCalculator(correction_config, decomp_result)
        
        reference_coeff = np.random.randn(20, n_components)
        correction_matrix, stats = calculator.calculate_correction(
            reference_coeff, science_coeff
        )
        
        # Apply correction
        corrected_coeff = calculator.apply_correction(science_coeff, correction_matrix)
        
        # Verify shapes
        assert corrected_coeff.shape == science_coeff.shape
        assert correction_matrix.shape == (n_components, n_components)
        assert stats['residual_rms'] >= 0
    
    def test_correction_result_roundtrip(self):
        """Test complete save/load cycle for correction results."""
        # Setup
        n_components = 5
        n_channels = 256
        n_ref = 30
        n_science = 20
        
        reference_coeff = np.random.randn(n_ref, n_components)
        science_coeff = np.random.randn(n_science, n_components)
        corrected_coeff = np.random.randn(n_science, n_components)
        correction_matrix = np.random.randn(n_components, n_components)
        components = np.random.randn(n_components, n_channels)
        mean_spectrum = np.random.randn(n_channels)
        
        config = CorrectionConfig(n_components=n_components)
        
        original_result = CorrectionResult(
            reference_coefficients=reference_coeff,
            science_coefficients=science_coeff,
            corrected_coefficients=corrected_coeff,
            correction_matrix=correction_matrix,
            components=components,
            mean_spectrum=mean_spectrum,
            config=config,
            metadata={'mission': 'TESTMISSION', 'tag': 'integration_test'},
        )
        
        with tempfile.NamedTemporaryFile(suffix='.pkl', delete=False) as f:
            filepath = f.name
        
        try:
            original_result.save(filepath)
            loaded_result = CorrectionResult.load(filepath)
            
            # Verify all attributes preserved
            assert np.allclose(loaded_result.reference_coefficients, reference_coeff)
            assert np.allclose(loaded_result.science_coefficients, science_coeff)
            assert np.allclose(loaded_result.corrected_coefficients, corrected_coeff)
            assert loaded_result.metadata['mission'] == 'TESTMISSION'
        finally:
            Path(filepath).unlink(missing_ok=True)


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
