# oi-zeigt Documentation

Tools to process and analyse the OI 63 µm SOFIA/GREAT archive.

---

## Quick navigation

| Document | Purpose |
|---|---|
| [START_HERE.md](START_HERE.md) | Project overview and orientation |
| [QUICK_REFERENCE.md](QUICK_REFERENCE.md) | Command cheat-sheet |
| [REPOSITORY_GUIDE.md](REPOSITORY_GUIDE.md) | Repository layout and git workflow |
| [VENV_SETUP.md](VENV_SETUP.md) | Virtual environment setup (cygrid, oi) |

---

## Subdirectories

### `pca_correction/`
PCA-based sky subtraction pipeline.

| File | Content |
|---|---|
| [pca_correct_workflow.md](pca_correction/pca_correct_workflow.md) | Detailed step-by-step walkthrough of `pca_correct` |
| [PCA_CORRECTION_WORKFLOW.md](pca_correction/PCA_CORRECTION_WORKFLOW.md) | Config integration and parameter reference |
| [PCA_DECOMPOSE_CONFIG.md](pca_correction/PCA_DECOMPOSE_CONFIG.md) | `pca_decompose` configuration guide |
| [PCA_QUICK_REFERENCE_CARD.md](pca_correction/PCA_QUICK_REFERENCE_CARD.md) | One-page quick reference |
| [MISSION_TELESCOPE_DECOMPOSITION.md](pca_correction/MISSION_TELESCOPE_DECOMPOSITION.md) | Per-mission/telescope decomposition |
| `reference/` | Channel trimming, variance, parameter details |
| `telluric_masking/` | Telluric line masking in PCA |

### `features/`
Documentation for individual pipeline features.

| File / folder | Content |
|---|---|
| [BLANK_VALUE_DETECTION.md](features/BLANK_VALUE_DETECTION.md) | Blank/missing value detection |
| [REDUCE_SPECTRA_METHODS.md](features/REDUCE_SPECTRA_METHODS.md) | Spectral reduction methods |
| `filter_fits/` | `filter_fits` command reference and examples |

### `velocity_axis/`
Velocity axis alignment, reconstruction, and reference.

### `combine_fits/`
Guide for combining multiple FITS files.

### `guides/`
Workflow walkthroughs (complete pipeline, before/after comparisons).

### `reference/`
PCA analysis indices, reduce_spectra options, test results.

### `implementation/`
Internal implementation notes, phase checklists, spectral axis work.

### `archive/`
Superseded status summaries and old delivery notes.
