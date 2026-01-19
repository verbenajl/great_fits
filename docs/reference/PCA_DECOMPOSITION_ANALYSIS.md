# PCA Decomposition Tool - Requirements Analysis

## Executive Summary

This document analyzes the existing GILDAS/CLASS-based PCA tools to guide implementation of a Python-native version for `oi_zeigt`. The tool performs PCA decomposition of atmospheric/instrumental components from SKYCHOPDIFF spectra and correction of M51CENTER spectra.

---

## 1. Overall Architecture & Workflow

### Two-Phase Process

The implementation consists of two main phases:

#### Phase 1: **DECOMPOSITION** (`pca_decompose.py`)
- Input: FITS spectra from GILDAS database (filtered by mission, telescope, scan)
- Process: Analyze SKYCHOPDIFF spectra for principal components
- Output: Pickled PCA models per flight/telescope pair
- Key: Ignore atmospheric line windows during analysis

#### Phase 2: **CORRECTION** (`pca_correct.py`)
- Input: PCA models + M51CENTER science spectra
- Process: Apply component analysis to remove atmospheric contamination
- Output: Corrected M51CENTER spectra in FITS format
- Key: Use noise ratios to selectively apply corrections

---

## 2. Configuration System

### Two-Tier Config Structure

#### **2.1 Main Configuration: `M51_pca_reduction.toml`**

Common parameters (apply to both phases):
- `pca_source`: Data source for PCA decomposition (DEFAULT: "SKYCHOPDIFF")
- `line_window`: Spectral range to ignore when detecting lines [ch_start, ch_end]
- `tag`: Identifier for decomposition run (for matching decompose ↔ correct)
- `limit_to_scan`: Debug option to process single scan only
- `limit_to_telescope`: Debug option to process single telescope only
- `smoothing_kernel_size`: Frequency-space smoothing of components (boxcar, default=3)
- `rolling_noise_window`: Window for rolling RMS estimate of components
- `decomposition_type`: "PCA", "SparsePCA", or "ICA" (default: "PCA")
- `do_scale`: Normalize by std dev before PCA (default: False)
- `z_score_cutoff`: Outlier detection threshold (default: False)

Decomposition-specific (`[decompose]`):
- `number_components`: Components to derive (default: 5)
- `add_sky_diff`: Include SKY-DIFF data in addition to SKYCHOPDIFF
- `noise_cutoff`: Z-score threshold for emission-free channels
- `scramble`: Random phase multiplication (experimental)

Correction-specific (`[correct]`):
- `cutoff`: Min variance explained ratio to use component
- `output_folder`: Diagnostic plots directory (default: "pca_plots")
- `export_components`: Write component spectra to output
- `global_noise_ratio_cutoff`: Noise ratio threshold (default: 15)
- `noise_ratio_cutoff`: Override automatic threshold (default: False)
- `line_kernel_size`: Blur kernel for line detection in waterfall plot (default: 61)
- `dump_exp_variance`: Write variance ratios to file
- `cut_coefficients`: Max coefficient value to use component (experimental)
- `dump_original_spec_array`: ASCII output of original spectra

**STATUS**: Can be merged into main `config.toml` with new `[pca]` section

#### **2.2 Flight Metadata: `mission_id_parameters.yml`**

Per-flight configuration (flight is mission date + flight number):

```yaml
2017-02-10_GR_F373:
  telluric_line_center: 592        # Channel where atmospheric line peaks
  telluric_line_width: 30          # Width of atmospheric line
  drop:                            # Data quality filters
    telescope:                     # Entire telescopes to exclude
      - LFAH_3
    scans:                         # Bad scans
      complete:                    # Exclude entire scan
        - 18529
      telescope:                   # Scan only on specific telescope
        LFAV_0:
          - 18611
  spikes:                          # Spike markers (experimental)
    - LFAH_3
```

**Key Features**:
- **Telluric line parameters**: Define atmospheric contamination per flight
- **Drop filters**: Exclude bad data at multiple levels (flight → telescope → scan → scan+telescope)
- **Flexible filtering**: Can specify quality cuts at different granularities

**STATUS**: Excellent reference; should integrate into main config under `[pca.missions]`

---

## 3. Data Filtering & Selection

### Selection Hierarchy

```
All FITS Files
    ↓
Filter by Flight (mission date + number)
    ↓
Filter by Telescope
    ↓
Filter by Scan
    ↓
Filter by Source (SKYCHOPDIFF vs M51CENTER)
    ↓
Filter by Quality (apply drop rules)
    ↓
Selected Spectra
```

### Implementation in Python

From `pca_utilities.py::create_index()`:
- Uses FITS headers to create pandas DataFrame with columns:
  - `number`, `telescope`, `scan`, `subscan`, `source`, `line`, `boff`, `loff`, `version`, `mission_id`
- Caches index in `.pca_index_*.pkl` to avoid recreating
- Groups by flight, telescope, scan, source for selective loading

**Key insight**: Use FITS structure directly (no GILDAS dependency)

---

## 4. Line Detection & Masking Strategy

### TWO TYPES OF LINES TO MASK

#### 4.1 TELLURIC LINES (Atmospheric contamination)
**What**: Atmospheric features that contaminate the signal
**Example**: Water vapor absorption at ~566 channels in some flights

**Where to find it**: `mission_id_parameters.yml`
```yaml
2017-02-10_GR_F373:
  telluric_line_center: 592        # Channel where atmospheric line peaks
  telluric_line_width: 30          # Width of atmospheric line (full width)
```

**How to mask**: 
- Hardcoded in configuration (pre-calibrated from observations)
- Always known per flight (consistent across observations)
- Mask region: `[telluric_center - width/2, telluric_center + width/2]`

**Why mask it**:
- During DECOMPOSITION: Prevents atmospheric features from entering PCA components
- During CORRECTION: Baseline-fitting doesn't try to fit atmospheric contamination

---

#### 4.2 SCIENCE LINES (CII in M51 case)
**What**: The actual astronomical feature we care about (CII emission from M51)
**Example**: CII 158 μm line in the spectral datacube

**Where to find it**: AUTO-DETECTION using waterfall plot + OpenCV
**Why auto-detect?**
- Different observing conditions → slightly different line positions
- Need to identify exact line channel range dynamically
- Waterfall plot reveals where the signal is concentrated

**How detection works**: From `pca_utilities.py::find_lines()`

**Step 1: Stack spectra vertically (waterfall plot)**
```
Create waterfall plot:
  Rows: Individual SKYCHOPDIFF spectra (to understand noise baseline)
  Cols: Frequency channels
  Values: Intensity (background noise level)
  
  Goal: See which channels have consistent signal across all spectra
```

**Step 2: Blur image to find concentrated regions**
```python
# Gaussian blur with kernel_size=line_kernel_size (default: 61)
# This smooths the waterfall plot → finds contiguous line regions
blurred = cv2.GaussianBlur(img, (kernel_size, kernel_size), 0)
```

**Step 3: Adaptive thresholding**
```python
# Find regions where intensity > mean + N*std
# These regions show where the CII line dominates
threshold = cv2.threshold(blurred, mean + cutoff_std * std, 1, cv.THRESH_BINARY)
```

**Step 4: Iterative refinement** (to remove false positives)
```python
# Mask detected line regions
# Recalculate mean/std using only "noise" regions (outside thresholded areas)
# Repeat thresholding with new statistics
# Continue until statistics converge
```

**Step 5: Contour detection**
```python
# Find spatial extent of detected regions
contours = cv2.findContours(threshold, cv2.RETR_TREE, cv2.CHAIN_APPROX_NONE)

# Result: Precise channel range of CII line in SKYCHOPDIFF spectra
```

**Result**: Binary mask showing which channels contain the CII science line

### Integration with PCA

**DECOMPOSITION Phase**:
```
SKYCHOPDIFF spectra
    ↓
Mask telluric line (from mission_id_parameters.yml)
    ↓
Auto-detect science line (CII) using waterfall + opencv
    ↓
Create combined mask: exclude BOTH telluric AND science lines
    ↓
Baseline-fit only in line-free regions
    ↓
PCA decomposition of atmospheric components only
    (science line signal is NOT in the PCA components)
```

**CORRECTION Phase**:
```
M51CENTER spectra
    ↓
Mask telluric line (from mission_id_parameters.yml)
    ↓
Use same science line mask detected during decomposition
    ↓
Apply PCA correction (subtract atmospheric components)
    BUT: Don't touch the science line channels (already masked)
    ↓
Result: CII line preserved, atmospheric lines removed
```

**For Python implementation**:
- **Telluric line**: Always from `mission_id_parameters.yml` (config-based)
- **Science line**: Auto-detected using waterfall plot + opencv (data-driven)
- Both masks applied during both decomposition and correction
- Configuration parameter: `enable_science_line_detection` = true/false

---

## 5. PCA Decomposition Process

### Data Preparation

From `pca_decompose.py::InputSpectra` and `create_components()`:

```
1. Filter spectra by mission/telescope/source
2. For each spectrum:
   a. Baseline (fit polynomial, exclude line window)
   b. Apply line mask (zero out atmospheric regions)
   c. Normalize if do_scale=True (divide by std dev)
3. Stack into 2D array: (n_spectra, n_channels)
```

### Decomposition

From `pca_decompose.py::PCADecomposition.compute_PCA()`:

```python
from sklearn.decomposition import PCA

# Input: normalized spectra array
pca = PCA(n_components=number_components - 1)
pca.fit(spectra_array)

# Output: 
#   pca.components_      (n_components-1, n_channels) - eigenspectra
#   pca.explained_variance_ratio_ (n_components-1,) - variance per component
#   pca.mean_            (n_channels,) - mean spectrum
```

### Post-Processing

**Smoothing** (`smooth_pca_components()`):
```python
# Frequency-space smoothing with boxcar filter
kernel = np.ones(smoothing_kernel_size) / smoothing_kernel_size
smoothed_components = [np.convolve(comp, kernel, mode='same') 
                       for comp in components]
```

**Global structure evaluation** (`evaluate_global_structure()`):
```python
# For each component, compute signal-to-noise ratio
for component in components:
    # RMS noise from rolling window
    rolling_rms = np.std(rolling_window(component, window_size), axis=1)
    mean_rolling_rms = rolling_rms.mean()
    
    # Global noise (overall std)
    global_noise = component.std()
    
    # Structure metric
    structure_metric = global_noise / mean_rolling_rms
    global_structure.append(structure_metric)
```

### Output Format

Pickle file: `<tag>___<telescope>___<mission_id>.pkl`

Contains: `PCADecomposition` object with:
- `pca`: sklearn.decomposition.PCA instance
- `pca_comp`: (n_components, n_channels) including mean
- `global_structure`: (n_components-1,) structure metrics
- Config dictionary

---

## 6. PCA Correction Process

### Phase 1: Component Selection

From `pca_correct.py` logic:

**For each component**, calculate:
```python
# Noise ratio: intrinsic noise of spectrum vs component amplitude
noise_ratio = spectrum_noise / component_amplitude

# Scaled noise ratio: after correcting for component strength
noise_ratio_scaled = noise_ratio * coefficient_of_fit
```

**Selection criteria**:
1. `cutoff`: Min variance explained ratio
2. `global_noise_ratio_cutoff`: Max noise ratio (default: 15)
   - KDE of noise ratios across all spectra
   - Find first peak of distribution
   - Components above threshold are used
3. `noise_ratio_cutoff`: Manual override
4. `cut_coefficients`: Max fit coefficient (experimental)

### Phase 2: Spectrum Correction

For each M51CENTER spectrum:

```python
1. Baseline fit (excluding line window)
2. For each selected component:
   a. Calculate least-squares fit coefficient
   b. Compute noise ratio
   c. If noise_ratio < threshold:
      - Subtract fitted_component from spectrum
      - Keep spectrum
   d. Else:
      - Skip this component (don't correct)
      - Mark as "skipped" in diagnostic output
3. Return corrected spectrum
```

### Diagnostics & Plots

From `plot_pca_decomposition()` and `plot_example_correction()`:

**1. Component plot**: Stack eigenspectra with:
   - Mean spectrum (red)
   - Each component colored by selected/excluded status
   - Variance explained per component
   - Global structure metric per component

**2. Correction example**: Show:
   - Original spectrum (black)
   - Each correction step (blue for component)
   - Final corrected spectrum (gray)
   - Diagnostic info: coefficient, noise ratios, skip status

---

## 7. Error Handling & Edge Cases

### Custom Exception Classes

From `pca_errors.py`:

```python
class InvalidOption(Exception):           # Bad config value
class MissingMandatoryOption(Exception):  # Missing required param
class LinePresent(Exception):             # Line detected in spectrum
class NoDataFound(Exception):              # No matching spectra
class NoScienceSourceFound(Exception):     # No M51CENTER in scan
class MultipleScienceSourcesFound(Exception):  # Multiple sources
```

### Handling Edge Cases

**Insufficient data**:
- If `n_components` > number of spectra available
- Action: Reduce n_components, warn user, write empty pickle

**No data for flight/telescope**:
- Raise `NoDataFound`, skip to next group

**Multiple science sources in scan**:
- Raise `MultipleScienceSourcesFound`, investigate

**Line detection ambiguity**:
- Use configuration parameters as fallback
- Can override with `telluric_line_center` ± `telluric_line_width`

---

## 8. Key Implementation Decisions for Python Port

### What to Keep (From GILDAS)

1. **Two-phase architecture**: Decompose → Correct
2. **Multi-level data filtering**: Mission → Telescope → Scan
3. **Flexible line masking**: Config + auto-detection
4. **Noise ratio-based component selection**
5. **Diagnostic plotting**

### What to Change (For Python/FITS)

| GILDAS | Python/FITS |
|--------|-----------|
| pyclass FITS access | astropy.io.fits |
| pgutils utilities | Custom Python utilities |
| GILDAS baselines | scipy.signal for polynomial fits |
| CLASS variables | FITS headers + pandas DataFrames |
| sicparse command parser | argparse + config files |
| Jinja2 templates | Simple config file generation |

### What to Add

1. **Direct FITS file handling**: Read spectra from local FITS files
2. **Output format flexibility**: FITS or ASCII
3. **Visualization**: matplotlib for diagnostic plots
4. **Modularity**: Separate utilities, decomposition, correction into distinct modules

---

## 9. Configuration Migration Plan

### Current Structure (GILDAS)

```toml
[common]
pca_source = "SKYCHOPDIFF"
line_window = [450, 550]
tag = "snakemake_decomposition"
...

[decompose]
number_components = 5
...

[correct]
output_folder = "pca_plots"
...
```

### Proposed Structure (Main config.toml)

```toml
[pca.common]
enabled = true
pca_source = "SKYCHOPDIFF"
line_window = [450, 550]
tag = "decomposition_v1"
smoothing_kernel_size = 3
rolling_noise_window = 11
decomposition_type = "PCA"

[pca.decompose]
number_components = 5
add_sky_diff = false
noise_cutoff = false
scramble = false

[pca.correct]
cutoff = false
output_folder = "pca_plots"
export_components = false
global_noise_ratio_cutoff = 15
line_kernel_size = 61

[pca.missions]
# Import mission parameters from mission_id_parameters.yml
# OR define here inline
"2017-02-10_GR_F373" = {
  telluric_line_center = 592
  telluric_line_width = 30
  drop = { telescope = ["LFAH_3"] }
}
```

---

## 10. Implementation Roadmap

### Phase 1: Foundation
1. ✅ Error classes
2. ✅ Configuration system (config.toml integration)
3. ✅ FITS utilities (reading, header parsing, indexing)
4. ✅ Baseline fitting (polynomial regression)

### Phase 2: Decomposition
1. Data filtering by mission/telescope/scan
2. Line masking (config-based + opencv auto-detection)
3. PCA decomposition (sklearn wrapper)
4. Component smoothing & evaluation

### Phase 3: Correction
1. Component selection (noise ratio thresholding)
2. Spectrum correction (least-squares fitting)
3. Diagnostic output (plots, diagnostics)
4. Output formatting (FITS + ASCII)

### Phase 4: Integration
1. CLI commands for decompose/correct
2. Workflow automation (Snakemake integration)
3. Caching system (pickle storage)
4. Documentation & examples

---

## 11. Dependencies

### Required
- `numpy`: Array operations
- `scipy`: Polynomial fitting, statistics
- `scikit-learn`: PCA decomposition
- `astropy`: FITS I/O
- `pandas`: DataFrame indexing
- `opencv-python`: Line detection
- `pyyaml`: Mission metadata
- `toml`: Configuration files
- `matplotlib`: Plotting

### Optional
- `scikit-learn[sparse]`: Sparse PCA variant
- `joblib`: Parallel processing
- `tqdm`: Progress bars

---

## 12. Data Flow Diagram

```
┌─────────────────────────────────────────────────────────┐
│                  FITS Data Files                         │
│            (mission date + telescope + scan)             │
└──────────────────────────┬──────────────────────────────┘
                           │
                           ▼
        ┌──────────────────────────────────────┐
        │     Create Index (pandas DF)          │
        │  - Extract headers (telescope, scan)  │
        │  - Identify missions                  │
        │  - Cache in .pkl                      │
        └──────┬───────────────────────────────┘
               │
        ┌──────▼────────────────────────────────┐
        │  DECOMPOSITION PHASE                  │
        │  (Process SKYCHOPDIFF spectra)        │
        ├──────────────────────────────────────┤
        │ 1. Filter by flight/telescope        │
        │ 2. Load SKYCHOPDIFF spectra          │
        │ 3. Detect line window                │
        │ 4. Mask atmospheric lines            │
        │ 5. Stack spectra vertically          │
        │ 6. Run PCA decomposition             │
        │ 7. Smooth components                 │
        │ 8. Evaluate global structure         │
        │ 9. Save PCA model (pickle)           │
        └──────┬───────────────────────────────┘
               │
        ┌──────▼────────────────────────────────┐
        │  CORRECTION PHASE                     │
        │  (Apply to M51CENTER spectra)        │
        ├──────────────────────────────────────┤
        │ 1. Load PCA model for flight         │
        │ 2. Filter by flight/telescope        │
        │ 3. Load M51CENTER spectra            │
        │ 4. For each spectrum:                │
        │    a. Baseline (no line window)      │
        │    b. For each component:            │
        │       - Fit component to spectrum    │
        │       - Calculate noise ratio        │
        │       - If ratio < threshold         │
        │         subtract component           │
        │ 5. Generate diagnostics              │
        │ 6. Save corrected spectra (FITS)     │
        └──────┬───────────────────────────────┘
               │
               ▼
    ┌─────────────────────────────────┐
    │  Corrected M51CENTER Spectra    │
    │  (atmospheric lines removed)     │
    └─────────────────────────────────┘
```

---

## 13. Summary: Key Takeaways

✅ **Two-phase architecture**: Decompose atmospheric components, then apply to science data

✅ **Flexible filtering**: Multi-level selection (flight → telescope → scan → source)

✅ **Smart line masking**: Config-based + opencv waterfall analysis

✅ **Noise-ratio driven**: Selective component application based on signal-to-noise

✅ **Comprehensive diagnostics**: Plots, variance tracking, correction examples

✅ **Modular design**: Separate concerns into decomposition, correction, utilities

✅ **Configuration-driven**: All parameters in TOML files, no hardcoding

---

## 14. Questions for Clarification

1. **Auto line detection**: Always use opencv waterfall method, or always use config?
   - Recommendation: Try both, use opencv if no config, override with config if provided

2. **Correction scope**: Only M51CENTER, or apply to other sources too?
   - Current code: M51CENTER only, could generalize

3. **Output format**: FITS only, or also ASCII spectra?
   - Current code: Both options available

4. **Parallel processing**: Process multiple flights/telescopes in parallel?
   - Current code: Sequential, could add joblib

5. **Caching**: Keep pickle files for future use?
   - Current code: Yes, cached in `.pca_pickled_objects/` folder

---

This analysis is ready for implementation! Let me know if you want to proceed with any specific phase.
