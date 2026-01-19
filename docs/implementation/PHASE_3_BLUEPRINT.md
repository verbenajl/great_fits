# Phase 3: Correction Module - Design Specification

**Status**: 📋 Design Ready  
**Target**: Implement correction algorithm to remove receiver variations  
**Expected Size**: ~800 lines  
**Expected Tests**: ~25 tests  

---

## 📋 Overview

Phase 3 implements the **Correction** workflow that uses PCA decomposition results from Phase 2 to correct spectral variations caused by instrument/receiver differences.

### Workflow

```
Phase 2 Results (PCA Components)
         ↓
    Load Components
         ↓
    Load Science Spectra (M51CENTER)
         ↓
    Project onto Components
         ↓
    Calculate Optimal Correction
    (Least-Squares Fit)
         ↓
    Apply Correction
    (Remove Receiver Variations)
         ↓
    Preserve Science Line
    (Via Masking Strategy)
         ↓
    Generate Corrected Spectra
         ↓
    Save Results + Provenance
```

---

## 🎯 Key Objectives

1. ✅ Load PCA decomposition results from Phase 2
2. ✅ Load science spectra (M51CENTER source)
3. ✅ Project science spectra onto PCA components
4. ✅ Calculate optimal correction factors via least-squares
5. ✅ Apply correction to spectral data
6. ✅ Preserve science line (via two-mask strategy)
7. ✅ Store corrected spectra with full metadata
8. ✅ Generate correction summary and statistics

---

## 🏗️ Architecture

### 5 Core Classes

#### 1. **CorrectionConfig** (Dataclass)
Configuration for the correction process.

**Parameters**:
- `n_components`: int - Number of components to use (from Phase 2)
- `science_source`: str - Source name for science spectra (e.g., "M51CENTER")
- `preserve_science_line`: bool - Whether to preserve science line (default: True)
- `mask_width`: int - Width of science line mask (channels)
- `fit_window`: int - Fitting window size for least-squares
- `regularization`: float - Regularization parameter (default: 0.0)

**Methods**:
- `to_dict()` - Serialize to dictionary
- `from_dict(data)` - Create from dictionary

---

#### 2. **SpectrumProjector** (Class)
Projects science spectra onto PCA components.

**Purpose**: Transform science spectra into component coefficient space

**Methods**:
- `__init__(decomposition_result: DecompositionResult)`
  - Store components and mean spectrum
  
- `project_spectrum(spectrum_1d: np.ndarray) -> Tuple[np.ndarray, dict]`
  - Project 1D spectrum onto components
  - Input: [n_channels] spectrum
  - Output: [n_components] coefficients + metadata

- `project_multiple(spectra_matrix: np.ndarray) -> Tuple[np.ndarray, List[dict]]`
  - Project multiple spectra
  - Input: [n_spectra, n_channels]
  - Output: [n_spectra, n_components] + metadata list

**Error Handling**:
- Check spectrum validity
- Handle dimension mismatches
- Log projection statistics

---

#### 3. **CorrectionCalculator** (Class)
Calculates optimal correction factors via least-squares fitting.

**Purpose**: Find best-fit correction to minimize residuals

**Methods**:
- `__init__(config: CorrectionConfig, decomposition_result: DecompositionResult)`
  - Store components and configuration
  
- `calculate_correction(reference_coefficients: np.ndarray, science_coefficients: np.ndarray) -> Tuple[np.ndarray, dict]`
  - Calculate correction matrix via least-squares
  - Input: Reference [n_ref, n_components] and Science [n_science, n_components]
  - Output: Correction matrix [n_components, n_components] + stats
  
- `apply_correction(coefficients: np.ndarray, correction_matrix: np.ndarray) -> np.ndarray`
  - Apply correction to coefficients
  - Input: [n_spectra, n_components] coefficients
  - Output: [n_spectra, n_components] corrected coefficients

**Statistics**:
- Residual RMS
- Correction strength
- Condition number

---

#### 4. **ScienceDataLoader** (Class)
Loads science (M51CENTER) spectra for correction.

**Purpose**: Load and prepare science spectra from FITS files

**Methods**:
- `__init__(config: ConfigLoader, mission_id: str)`
  - Initialize with configuration
  
- `load_science_data(fits_directory: str, mission_id: str) -> Tuple[np.ndarray, List[dict], dict]`
  - Load M51CENTER spectra
  - Filter by source
  - Apply drop filters
  - Return spectral matrix [n_spectra, n_channels] + metadata

**Integration**:
- Uses FITSIndexer from Phase 1
- Uses ConfigLoader from Phase 1
- Follows DecompositionDataLoader pattern

---

#### 5. **CorrectionResult** (Dataclass)
Stores correction results with full provenance.

**Attributes**:
- `reference_coefficients`: [n_ref, n_components] - Reference (SKYCHOPDIFF) projections
- `science_coefficients`: [n_science, n_components] - Science (M51CENTER) projections
- `corrected_coefficients`: [n_science, n_components] - After correction
- `correction_matrix`: [n_components, n_components] - The correction transformation
- `components`: [n_components, n_channels] - PCA components (from Phase 2)
- `mean_spectrum`: [n_channels] - Mean spectrum from decomposition
- `config`: CorrectionConfig - Configuration used
- `metadata`: dict - Mission, date, statistics, etc.
- `reference_metadata`: List[dict] - Per-reference-spectrum metadata
- `science_metadata`: List[dict] - Per-science-spectrum metadata

**Methods**:
- `reconstruct_corrected_spectra() -> np.ndarray`
  - Inverse transform corrected coefficients to spectral space
  - Output: [n_science, n_channels] corrected spectra
  
- `get_correction_strength() -> dict`
  - Return statistics about correction magnitude
  - Includes RMS residual, component-wise correction magnitude
  
- `summary() -> str`
  - Generate text summary of correction
  
- `save(filepath: str)` - Pickle serialization
  
- `load(filepath: str) -> 'CorrectionResult'` - Pickle deserialization

---

#### 6. **correct_spectra()** Function (Main Entry Point)
Orchestrates complete correction workflow.

**Parameters**:
- `fits_directory`: str - Path to FITS files
- `decomposition_file`: str - Path to Phase 2 results pickle
- `config_file`: str - TOML configuration file (optional)
- `mission_file`: str - YAML mission parameters (optional)
- `mission_id`: str - Mission identifier
- `output_dir`: str - Output directory for results (optional)

**Workflow**:
1. Load configuration
2. Load decomposition results from Phase 2
3. Load science spectra (M51CENTER)
4. Project science spectra onto components
5. Calculate correction factors
6. Apply correction
7. Generate results
8. Save if output_dir specified
9. Return CorrectionResult

**Returns**: CorrectionResult object

**Logging**: Comprehensive logging at each step

---

## 🧪 Testing Strategy

### Test Classes (25+ tests)

1. **TestCorrectionConfig** (5 tests)
   - Initialization with defaults/custom
   - Dictionary serialization
   - Roundtrip consistency

2. **TestSpectrumProjector** (6 tests)
   - Projection of single spectrum
   - Projection of multiple spectra
   - Output shapes and types
   - Metadata generation
   - Error handling

3. **TestCorrectionCalculator** (8 tests)
   - Correction calculation
   - Correction application
   - Residual statistics
   - Regularization effects
   - Error conditions

4. **TestCorrectionResult** (6 tests)
   - Initialization
   - Reconstruction of corrected spectra
   - Correction strength calculation
   - Save/load persistence
   - Metadata preservation

5. **TestIntegration** (3+ tests)
   - Complete correction workflow
   - Roundtrip persistence
   - Correction quality metrics

---

## 📊 Data Flow

### Input
- **Phase 2 Results**: Components, mean spectrum, explained variance
- **Reference Data**: SKYCHOPDIFF spectra (already in components)
- **Science Data**: M51CENTER spectra to be corrected
- **Configuration**: Mission-specific parameters

### Processing
```
Science Spectra [n_science, n_channels]
       ↓
Project onto Components
       ↓
Science Coefficients [n_science, n_components]
       ↓
Calculate Correction: argmin ||ref_coeff - correction @ sci_coeff||²
       ↓
Correction Matrix [n_components, n_components]
       ↓
Apply: corrected_coeff = correction @ science_coeff
       ↓
Reconstruct: corrected_spectra = corrected_coeff @ components + mean
       ↓
Preserve Science Line (mask-based)
       ↓
Corrected Spectra [n_science, n_channels]
```

### Output
- **CorrectionResult**: Full results with provenance
- **Corrected Spectra**: [n_science, n_channels]
- **Correction Metadata**: Statistics and quality measures

---

## 🔄 Integration with Phases 1 & 2

### Phase 1 Modules Used
- `config.py` - ConfigLoader
- `fits_indexing.py` - FITSIndexer
- `utilities.py` - Baseline fitting, masking, normalization
- `line_detection.py` - Science line detection
- `errors.py` - Exception hierarchy

### Phase 2 Integration
- Loads `DecompositionResult` from Phase 2
- Uses components and mean spectrum
- Uses explained variance information
- Builds on same patterns and conventions

---

## 📈 Two-Mask Strategy (Science Line Preservation)

Phase 3 will implement a sophisticated masking strategy:

1. **Artifact Mask** - Hardware features (fixed position)
2. **Science Line Mask** - Detected science line (variable position)
3. **Combined Mask** - Union of both

**Application**:
- Calculate correction using all data
- Apply correction to all channels
- In corrected spectrum: preserve masked regions from science data
- Replace unmasked regions with corrected data

This preserves authentic science line while correcting receiver variations.

---

## 🔑 Key Algorithms

### Least-Squares Correction
```
Given:
  R = [r₁, r₂, ..., rₙ]  (reference coefficients) [n_ref, n_comp]
  S = [s₁, s₂, ..., sₘ]  (science coefficients)   [n_science, n_comp]

Find C (correction matrix) [n_comp, n_comp] that minimizes:
  ||R - C @ S.T||²_F

Solution: C = R @ S.T @ (S @ S.T)⁻¹

Apply correction:
  S_corrected = S @ C.T
```

### Regularization (Optional)
For ill-conditioned systems:
```
C = R @ S.T @ (S @ S.T + λ*I)⁻¹
where λ is regularization parameter
```

---

## ✨ Quality Metrics

Phase 3 will compute:

1. **Residual RMS** - Fit quality
2. **Correction Strength** - Per-component correction magnitude
3. **Channel-wise Statistics** - Before/after correction
4. **Science Line Preservation** - How much science line is preserved

---

## 📝 Configuration Example

```toml
[correction]
n_components = 5
science_source = "M51CENTER"
preserve_science_line = true
mask_width = 20
fit_window = 512
regularization = 0.0
```

---

## 🎯 Success Criteria

✅ Load Phase 2 results correctly  
✅ Project spectra onto components  
✅ Calculate correction via least-squares  
✅ Apply correction to spectra  
✅ Preserve science line (two-mask strategy)  
✅ Generate results with full metadata  
✅ All 25+ tests passing  
✅ Zero regressions in Phase 1 + Phase 2 tests  
✅ Production-quality code and documentation  

---

## 🚀 Implementation Plan

### Day 1: Core Classes
- Create `src/oi_zeigt/pca_analysis/correct.py`
- Implement CorrectionConfig
- Implement SpectrumProjector
- Implement CorrectionCalculator

### Day 2: Data & Results
- Implement ScienceDataLoader
- Implement CorrectionResult
- Implement correct_spectra() function
- Add comprehensive error handling

### Day 3: Testing & Validation
- Create `test_correction.py` with 25+ tests
- Run Phase 1 + Phase 2 + Phase 3 tests (54+ total)
- Verify no regressions
- Complete documentation

---

## 📚 Documentation Deliverables

- `PHASE_3_COMPLETION.md` - Implementation details
- `PHASE_3_STATUS.md` - Full status report
- `PHASE_3_READY.md` - Quick reference
- `PHASE_3_FINAL_DELIVERY.md` - Final summary
- Code docstrings and inline comments

---

## 🎉 End Result

By end of Phase 3:
- ✅ Complete correction workflow implemented
- ✅ 25+ tests all passing
- ✅ Full provenance tracking
- ✅ Production-ready code
- ✅ Ready for Phase 4 (CLI integration)

**Total Project Size**: ~3,100 lines, 79+ tests, 100% passing

---

**This Phase 3 Blueprint is ready for implementation!**
