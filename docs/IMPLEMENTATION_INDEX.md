# PCA Decomposition Tool: Implementation Status & Documentation Index

**Project**: OI-ZEIGT PCA Analysis Module  
**Status**: ✅ Phase 1 Complete, Phase 2 Ready to Begin  
**Last Updated**: January 16, 2026  

## Quick Navigation

### Phase 1: Foundation (✅ COMPLETE)
All foundation modules have been implemented and syntax-validated.

- **[PHASE_1_COMPLETION_SUMMARY.md](PHASE_1_COMPLETION_SUMMARY.md)** - Overview of what was built
- **[PHASE_1_ARCHITECTURE.md](PHASE_1_ARCHITECTURE.md)** - Module organization and dependencies
- **Implementation Files**:
  - `src/oi_zeigt/pca_analysis/errors.py` - Exception hierarchy (14 classes)
  - `src/oi_zeigt/pca_analysis/config.py` - Configuration loading system
  - `src/oi_zeigt/pca_analysis/utilities.py` - Helper functions
  - `src/oi_zeigt/pca_analysis/fits_indexing.py` - FITS file handling
  - `src/oi_zeigt/pca_analysis/line_detection.py` - Line detection & masking

### Phase 2: Decomposition (📋 BLUEPRINT READY)
Detailed specification for PCA decomposition implementation.

- **[PHASE_2_BLUEPRINT.md](PHASE_2_BLUEPRINT.md)** - Full implementation plan for Phase 2
- **Status**: Ready to begin implementation
- **Deliverable**: `src/oi_zeigt/pca_analysis/decompose.py`

### Phase 3: Correction (📋 DESIGN COMPLETE)
See [PCA_CORRECTION_CLARIFICATION.md](#phase-3-correction-reference) for design.

- **Deliverable**: `src/oi_zeigt/pca_analysis/correct.py`

### Phase 4: CLI & Integration (📋 DESIGN COMPLETE)
Integration with main CLI system.

- **Deliverable**: Updates to `src/oi_zeigt/cli.py`

---

## Reference Documents

### Analysis & Design Documents

These documents from earlier analysis phases contain the full design:

- **[PCA_DECOMPOSITION_ANALYSIS.md](PCA_DECOMPOSITION_ANALYSIS.md)** (~600 lines)
  - Original analysis of GILDAS reference implementation
  - Detailed algorithm explanation
  - Component identification and responsibilities

- **[PCA_LINE_DETECTION_CLARIFICATION.md](PCA_LINE_DETECTION_CLARIFICATION.md)** (~300 lines)
  - Clarification on science line detection
  - Two-mask strategy explanation
  - Line detection methodology

- **[PCA_CORRECTION_CLARIFICATION.md](PCA_CORRECTION_CLARIFICATION.md)** (~300 lines)
  - What PCA corrects for (receiver/electronics, not atmosphere)
  - Phase 3 correction algorithm
  - Least-squares fitting methodology

- **[PCA_DECOMPOSITION_TOOL_SUMMARY.md](PCA_DECOMPOSITION_TOOL_SUMMARY.md)** (~600 lines)
  - Complete feature summary
  - Configuration system documentation
  - Testing strategy

- **[PCA_ANALYSIS_INDEX.md](PCA_ANALYSIS_INDEX.md)** (~200 lines)
  - Index of all analysis documents
  - Quick reference to key concepts

### Implementation Documents

- **[PHASE_1_COMPLETION_SUMMARY.md](PHASE_1_COMPLETION_SUMMARY.md)** (~300 lines)
  - What was built in Phase 1
  - Module descriptions
  - Code statistics

- **[PHASE_1_ARCHITECTURE.md](PHASE_1_ARCHITECTURE.md)** (~300 lines)
  - Module dependency graph
  - Data flow diagrams
  - Function call patterns

- **[PHASE_2_BLUEPRINT.md](PHASE_2_BLUEPRINT.md)** (~400 lines)
  - Decomposition phase specification
  - Class and function designs
  - Implementation checklist

---

## Phase Status Details

### ✅ Phase 1: Foundation (COMPLETE)

**5 Core Modules Created** (1,549 lines of code):

1. **errors.py** (63 lines)
   - 14 exception classes in hierarchy
   - Handles all error cases from Phases 1-4
   - ✅ Syntax validated

2. **config.py** (326 lines)
   - ConfigLoader class with 8 methods
   - Supports TOML + YAML configuration
   - Nested key access with dot notation
   - ✅ Syntax validated

3. **utilities.py** (340 lines)
   - 12 core utility functions
   - Baseline fitting, masking, normalization
   - Spectrum validation and statistics
   - ✅ Syntax validated

4. **fits_indexing.py** (430 lines)
   - FITSIndexer class with 11 methods
   - FITS file scanning and indexing
   - Hierarchical filtering (mission → telescope → source)
   - ✅ Syntax validated

5. **line_detection.py** (390 lines)
   - 3 science line detection methods
   - LineDetectionConfig class
   - Masking utilities (artifact + science line)
   - ✅ Syntax validated

**Test Coverage**:
- 40+ unit tests defined in `test_pca_phase1.py`
- All syntax passes Python validation

---

### 📋 Phase 2: Decomposition (READY TO BEGIN)

**What Needs to be Built**:
- `decompose.py` (~600-800 lines estimated)
- Unit tests for decomposition module
- Integration tests with Phase 1

**Key Components** (from PHASE_2_BLUEPRINT.md):
1. DecompositionConfig - Configuration holder
2. SpectrumPreparator - Prepare spectra (baseline, mask, normalize)
3. DecompositionDataLoader - Load FITS data
4. PCADecomposer - Perform PCA fitting
5. DecompositionResult - Store and save results

**Entry Point**:
```python
result = decompose_spectra(
    fits_directory=Path('/data'),
    config_file=Path('config.toml'),
    mission_file=Path('missions.yml')
)
```

---

### 📋 Phase 3: Correction (DESIGNED)

**What Needs to be Built**:
- `correct.py` (~600-800 lines estimated)
- Unit tests for correction module
- Integration tests with Phases 1-2

**Key Components** (from design documents):
1. CorrectionDataLoader - Load science spectra + components
2. CorrectionCalculator - Calculate optimal correction (least-squares)
3. CorrectionApplier - Apply correction while preserving science line
4. CorrectionResult - Store corrected spectra

**Entry Point**:
```python
corrected = correct_spectra(
    decomposition_result=result_from_phase2,
    fits_directory=Path('/data'),
    config_file=Path('config.toml')
)
```

---

### 📋 Phase 4: CLI & Integration (DESIGNED)

**What Needs to be Done**:
- Add `--decompose` and `--correct` subcommands to CLI
- Integrate with configuration system
- Add output reporting and visualization
- Update `src/oi_zeigt/cli.py`

---

## How to Use These Documents

### For Implementation (Phase 2 Start)

1. **Read PHASE_2_BLUEPRINT.md** (10 min)
   - Understand what needs to be built
   - See the class structures
   - Check the implementation checklist

2. **Reference Phase 1 Code** (15 min)
   - Look at `config.py` for ConfigLoader pattern
   - Look at `fits_indexing.py` for data loading pattern
   - Look at `errors.py` for exception handling

3. **Start Implementation**
   - Create `decompose.py`
   - Follow the checklist in PHASE_2_BLUEPRINT.md
   - Use Phase 1 modules as dependencies

4. **Test**
   - Create `test_decomposition.py`
   - Test each class independently
   - Test full pipeline
   - Use Phase 1 test patterns from `test_pca_phase1.py`

### For Understanding the Algorithm

1. **Read PCA_DECOMPOSITION_ANALYSIS.md** (15 min)
   - Overview of what the tool does
   - Algorithm explanation

2. **Read PCA_LINE_DETECTION_CLARIFICATION.md** (10 min)
   - Understand the two-mask strategy
   - Learn about science line detection

3. **Read PCA_CORRECTION_CLARIFICATION.md** (10 min)
   - Understand what gets corrected (receiver, not atmosphere)
   - Learn the correction algorithm

### For Configuration Reference

1. **Read PHASE_2_BLUEPRINT.md** - Configuration section
2. **Look at `config.py`** - How to load and access config
3. **Check sample files**:
   - `src/oi_zeigt/pca_analysis/M51_pca_reduction.toml`
   - `src/oi_zeigt/pca_analysis/mission_id_parameters.yml`

---

## Code Statistics

| Category | Value |
|----------|-------|
| **Phase 1 Lines of Code** | 1,549 |
| **Phase 1 Classes** | 17 |
| **Phase 1 Functions** | 23 |
| **Phase 1 Unit Tests** | 40+ |
| **Exception Classes** | 14 |
| **FITS Indexing Methods** | 11 |
| **Line Detection Methods** | 3 |
| **Baseline Fitting** | 1 (with sigma-clipping) |
| **Syntax Validation** | ✅ 5/5 modules pass |

---

## Key Technical Decisions

1. **Two-Mask Strategy**
   - Artifact mask: Fixed (from config, represents hardware artifact)
   - Science line mask: Auto-detected (from data, represents CII line)
   - Both masks exclude regions during baseline fitting

2. **Configuration System**
   - TOML for main configuration (pca parameters)
   - YAML for mission-specific parameters
   - Dot notation for nested access: `'pca.decompose.number_components'`

3. **FITS Indexing**
   - Pandas DataFrame for flexible filtering
   - Chainable filter methods (mission → telescope → source)
   - Lazy loading (only read FITS when requested)

4. **Line Detection**
   - Three methods: waterfall (default), max_intensity, centroid
   - Confidence metrics for detection quality
   - Validation to ensure detection is reasonable

5. **Error Handling**
   - Hierarchy of exceptions (base → specific)
   - Easy to catch all PCA errors or specific types
   - Informative error messages with context

---

## Dependencies

### Phase 1 Dependencies
- `numpy` - Array operations
- `scipy` - Scientific computing
- `astropy.io.fits` - FITS file handling
- `astropy.table` - Table operations
- `pandas` - DataFrame for indexing
- `toml` (optional) - TOML configuration
- `pyyaml` (optional) - YAML mission parameters

### Phase 2 (Will Add)
- `sklearn.decomposition` - PCA

### Phase 3 (Will Add)
- `scipy.optimize` - Least-squares fitting

---

## Next Steps

### Immediate (Ready Now)
- [ ] Read PHASE_2_BLUEPRINT.md
- [ ] Review Phase 1 code as implementation reference
- [ ] Prepare test FITS files or mock data

### Phase 2 Implementation
- [ ] Create `decompose.py`
- [ ] Implement DecompositionDataLoader
- [ ] Implement SpectrumPreparator
- [ ] Implement PCADecomposer
- [ ] Implement DecompositionResult
- [ ] Create unit tests
- [ ] Create integration tests
- [ ] Document with docstrings

### Phase 2 Completion → Phase 3
- [ ] Verify Phase 2 produces correct components
- [ ] Create PHASE_3_BLUEPRINT.md
- [ ] Implement correction module

---

## Support Files

All original GILDAS reference implementations are preserved in `src/oi_zeigt/pca_analysis/`:
- `pca_decompose.py` - Reference decomposition implementation
- `pca_correct.py` - Reference correction implementation
- `pca_utilities.py` - Reference utilities
- `pca_errors.py` - Reference error handling
- `M51_pca_reduction.toml` - Reference configuration
- `mission_id_parameters.yml` - Reference mission parameters

These can be referenced during implementation if needed.

---

## Questions or Issues?

Refer to the appropriate document:
- **Algorithm questions**: PCA_DECOMPOSITION_ANALYSIS.md
- **Line detection questions**: PCA_LINE_DETECTION_CLARIFICATION.md
- **Correction questions**: PCA_CORRECTION_CLARIFICATION.md
- **Implementation questions**: PHASE_1_ARCHITECTURE.md, PHASE_2_BLUEPRINT.md
- **Configuration questions**: PHASE_2_BLUEPRINT.md (config section)
- **API reference**: Individual module docstrings in `src/oi_zeigt/pca_analysis/`

---

**Ready to proceed with Phase 2! 🚀**
