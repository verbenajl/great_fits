# OI-Zeigt Git Repository Guide

## Repository Overview

**Project**: OI-Zeigt - Optical Interferometry Analysis via PCA  
**Repository Location**: `/home/verbena/software/oi_zeigt`  
**VCS**: Git  
**Status**: ✅ Version controlled with clean history

### Quick Statistics
- **Total Commits**: 9 (organized by phase)
- **Total Files**: 120
- **Total Lines of Code**: ~30,000
- **Test Coverage**: 77/77 tests passing (100%)
- **Phases Completed**: 3 (100% implementation)

---

## Commit History

The repository has been organized with **9 strategic commits**, each focusing on a specific component:

### 1. Phase 1: Baseline Utilities, Configuration, and Line Detection
**Commit**: `4578947`  
**Files**: 13 | **Insertions**: 7,156

Core foundational module with:
- Configuration loader and parameter management
- FITS file detection and catalog utilities
- Line detection for spectral analysis
- Error handling framework
- **Tests**: 26 tests (all passing)

**Key Files**:
- `src/oi_zeigt/pca_analysis/config.py`
- `src/oi_zeigt/pca_analysis/line_detection.py`
- `src/oi_zeigt/pca_analysis/utilities.py`
- `test_pca_phase1.py`

---

### 2. Phase 2: PCA Decomposition Engine
**Commit**: `3109286`  
**Files**: 2 | **Insertions**: 1,175

Principal Component Analysis core:
- Robust decomposition using SVD with automatic rank detection
- Statistical analysis of explained variance
- Component-based spectral representation
- Serialization and persistence
- **Tests**: 24 tests (all passing)

**Key Files**:
- `src/oi_zeigt/pca_analysis/decompose.py`
- `test_decomposition.py`

---

### 3. Phase 3: Spectral Correction via PCA
**Commit**: `9bda989`  
**Files**: 2 | **Insertions**: 1,325

Complete correction framework:
- Mean-centered least-squares alignment algorithm
- SVD-based principal component alignment with automatic scaling
- Science data loading and quality filtering
- Correction result persistence and reconstruction
- **Tests**: 23 tests (all passing)
- **Integration**: All phases working together (77/77 tests passing)

**Key Files**:
- `src/oi_zeigt/pca_analysis/correct.py`
- `test_correction.py`

---

### 4. Documentation and Project Structure
**Commit**: `5417390`  
**Files**: 61 | **Insertions**: 16,196

Comprehensive documentation:
- Phase 1, 2, 3 implementation guides
- Architecture blueprints and design documents
- PCA analysis reference materials
- Project status reports and summaries
- Quick reference guides
- README files for navigation

**Key Folders**:
- `docs/implementation/` - Phase-specific documentation (21 files)
- `docs/reference/` - Analysis and reference materials (10 files)
- `docs/guides/` - User guides and workflows
- `docs/archive/` - Historical documentation

---

### 5. Project Configuration and Setup Files
**Commit**: `f0146c5`  
**Files**: 5 | **Insertions**: 163

Project infrastructure:
- `pyproject.toml` - Modern Python project configuration
- `setup.py` - Package installation script
- `MANIFEST.in` - Distribution manifest
- `config.toml` - Application configuration
- `.gitignore` - Version control exclusions

---

### 6. Core Module Infrastructure and Sample Data
**Commit**: `1cac986`  
**Files**: 6 | **Insertions**: 1,205

Module essentials:
- PCA analysis module initialization and exports
- Error handling and custom exceptions
- FITS file indexing utilities
- Mission ID parameters for spectral calibration
- Sample FITS file list for testing

**Key Files**:
- `src/oi_zeigt/pca_analysis/__init__.py`
- `src/oi_zeigt/pca_analysis/errors.py`
- `src/oi_zeigt/pca_analysis/core.py`
- `src/oi_zeigt/pca_analysis/fits_indexing.py`
- `src/oi_zeigt/pca_analysis/mission_id_parameters.yml`

---

### 7. Development and Diagnostic Utilities
**Commit**: `b4364da`  
**Files**: 10 | **Insertions**: 1,222

Testing and validation tools:
- CRPIX1 preservation tests
- Mapping integration tests
- Spectral axis validation tests
- Velocity axis analysis tests
- FITS inspection tools
- Diagnostic utilities for development

**Key Files**:
- `test_crpix1_preservation.py`
- `test_map_integrated.py`
- `test_mapping_multi_hdu.py`
- `test_spectral_axis_*.py`
- `test_velocity_axis.py`
- `inspect_fits*.py`
- `check_spec.py`

---

### 8. Analysis and Debugging Scripts
**Commit**: `60b7915`  
**Files**: 12 | **Insertions**: 1,266

Development and analysis tools:
- FITS header analysis and inspection
- Reference channel investigation utilities
- Spectral axis and velocity axis analysis
- FITS file combining examples
- Visualization and debugging aids

**Key Files**:
- `analyze_reference_channel.py`
- `check_*.py` (multiple variants)
- `inspect_datacube_velocity.py`
- `investigate_reference_channel.py`
- `search_refchan.py`
- `visualize_velocity_axis.py`
- `example_combine_fits.sh`

---

### 9. Legacy Modules and Sample Data
**Commit**: `67fd254`  
**Files**: 9 | **Insertions**: 4,398

Archival and reference:
- Legacy PCA module implementations
- Broken CLI (for reference)
- Sample FITS data for testing
- Sample visualizations
- Configuration templates

**Key Files**:
- `src/oi_zeigt/pca_analysis/pca_*.py` (legacy implementations)
- `averaged_data.fits` (sample data)
- `spectra_with_blanks.png` (sample visualization)
- `velocity_axis_visualization.png`
- `M51_pca_reduction.toml` (config template)

---

## Viewing Commits

To see details of any commit:

```bash
# View full commit history with graphs
git log --graph --oneline

# View detailed commit information
git log -p <commit_hash>

# View changes in a specific commit
git show <commit_hash>

# View author statistics
git log --shortstat

# View commits for a specific file
git log -- <filepath>
```

---

## Branch Structure

Currently on **main branch**: `master`

To create a feature branch for development:

```bash
git checkout -b feature/<feature-name>
```

---

## Working with the Repository

### Cloning
```bash
git clone <repository-url> oi-zeigt
cd oi-zeigt
```

### Installation
```bash
# Install in development mode
pip install -e .

# Install with test dependencies
pip install -e ".[test]"
```

### Running Tests
```bash
# Run all tests
pytest

# Run specific phase tests
pytest test_pca_phase1.py       # Phase 1 (26 tests)
pytest test_decomposition.py    # Phase 2 (24 tests)
pytest test_correction.py       # Phase 3 (23 tests)

# Run with coverage
pytest --cov=src/oi_zeigt
```

### Making Changes

1. **Create a feature branch**:
   ```bash
   git checkout -b feature/my-feature
   ```

2. **Make your changes** and test them:
   ```bash
   pytest  # Run all tests
   ```

3. **Commit with clear messages**:
   ```bash
   git add <files>
   git commit -m "Brief description of changes"
   ```

4. **Push to remote** (when configured):
   ```bash
   git push origin feature/my-feature
   ```

---

## Repository Structure

```
oi_zeigt/
├── .git/                          # Git directory (version control)
├── .gitignore                     # Git exclusion rules
├── README.md                      # Project overview
├── setup.py                       # Package setup
├── pyproject.toml                 # Project configuration
├── MANIFEST.in                    # Distribution manifest
├── config.toml                    # Application configuration
│
├── src/oi_zeigt/                  # Main source code
│   ├── __init__.py
│   ├── basic_io.py
│   ├── cli.py
│   ├── mapping/                   # Mapping module
│   ├── reduction/                 # Reduction utilities
│   ├── statistics/                # Statistical analysis
│   └── pca_analysis/              # PCA analysis module
│       ├── __init__.py
│       ├── config.py              # Phase 1: Configuration
│       ├── line_detection.py      # Phase 1: Line detection
│       ├── utilities.py           # Phase 1: Utilities
│       ├── errors.py              # Error handling
│       ├── core.py                # Core functionality
│       ├── fits_indexing.py       # FITS utilities
│       ├── decompose.py           # Phase 2: Decomposition
│       ├── correct.py             # Phase 3: Correction
│       └── ...other modules       # Legacy implementations
│
├── tests/
│   ├── test_pca_phase1.py         # Phase 1 tests (26 tests)
│   ├── test_decomposition.py      # Phase 2 tests (24 tests)
│   ├── test_correction.py         # Phase 3 tests (23 tests)
│   └── ...other test files
│
├── docs/                          # Documentation
│   ├── README.md
│   ├── IMPLEMENTATION_INDEX.md
│   ├── implementation/            # Phase implementation docs
│   ├── reference/                 # Reference materials
│   ├── guides/                    # User guides
│   └── archive/                   # Historical docs
│
└── sample_fits_list.txt           # Sample data reference
```

---

## Commit Workflow

Each commit is organized by functionality:

| Commit | Component | Purpose | Test Status |
|--------|-----------|---------|-------------|
| 1 | Phase 1 | Baseline infrastructure | 26/26 ✅ |
| 2 | Phase 2 | Decomposition engine | 24/24 ✅ |
| 3 | Phase 3 | Correction framework | 23/23 ✅ |
| 4 | Documentation | Guides and references | N/A |
| 5 | Configuration | Setup and deployment | N/A |
| 6 | Core Infrastructure | Module essentials | N/A |
| 7 | Development Tools | Utilities and diagnostics | N/A |
| 8 | Analysis Scripts | Investigation tools | N/A |
| 9 | Legacy Code | Reference implementations | N/A |

---

## Adding Remote Repository

To add a remote repository (GitHub, GitLab, etc.):

```bash
# Add remote (replace with your URL)
git remote add origin https://github.com/username/oi-zeigt.git

# Verify remote
git remote -v

# Push to remote
git push -u origin master
```

---

## Useful Git Commands

```bash
# Check status
git status

# View unstaged changes
git diff

# View staged changes
git diff --staged

# View branch history
git log --oneline -10

# Create a tag
git tag -a v1.0.0 -m "Version 1.0.0"

# View tags
git tag -l

# Push tags
git push origin --tags

# List branches
git branch -a

# Switch branches
git checkout <branch-name>

# Merge a branch
git merge <branch-name>

# View who changed what
git blame <filepath>

# Search commit messages
git log --grep="keyword"

# Get commit statistics
git log --stat
```

---

## Next Steps

1. **Add Remote**: Set up GitHub/GitLab remote and push
2. **Create Tags**: Version releases with git tags
3. **Set up CI/CD**: Automate testing with GitHub Actions or similar
4. **Create Branches**: Use branches for feature development
5. **Documentation**: Keep docs/ folder updated with changes

---

## Support and Questions

- See `docs/IMPLEMENTATION_INDEX.md` for navigation guide
- See `docs/reference/` for technical analysis
- See `test_*.py` files for usage examples
- Run tests to validate functionality: `pytest`

---

**Repository Initialized**: January 19, 2026  
**Current Status**: ✅ Clean and ready for development  
**Last Updated**: $(date)
