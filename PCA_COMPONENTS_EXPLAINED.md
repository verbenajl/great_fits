# PCA Components Determination - Comprehensive Guide

## Overview

This document explains how PCA (Principal Component Analysis) components are determined in the OI-Zeigt project, with detailed code references and mathematical background.

**TL;DR:** PCA finds the main directions of variation in 13,986 sky spectra, extracting 5 principal directions (components) that explain 99.76% of the variance.

---

## Quick Summary

### The Algorithm in 3 Steps

#### 1. **COLLECT & CENTER**
```
Input: 13,986 spectra × 1229 channels

mean_spectrum = average of all 13,986 spectra
centered = each spectrum - mean_spectrum

Result: Data centered around zero
```

#### 2. **MEASURE CORRELATION**
```
covariance_matrix = how much do channels co-vary?

Shape: 1229 × 1229

Tells us: which wavelengths change together?
```

#### 3. **FIND PRINCIPAL DIRECTIONS**
```
Eigendecomposition of covariance matrix

Eigenvectors = principal directions (your components)
Eigenvalues = how much variance in each direction

Select top 5 eigenvectors

Result: 5 components × 1229 channels each
```

---

## What's in Your Pickle File

The file `decomposition_2017-02-01_GR_F_2017021_components.pkl` contains:

### ✓ components [5 × 1229]
- The 5 principal directions found by PCA
- Each row is one eigenvector (direction)
- Each column is the "weight" at that wavelength
- Range: -0.31 to +0.56
- **Memory:** 24 KB

### ✓ mean_spectrum [1229]
- Average of all 13,986 reference spectra
- **Essential for:** Centering future spectra, reconstructing spectra
- Range: 6.6 to 31.1 (intensity units)
- **Memory:** 4.8 KB

### ✓ explained_variance_ratio [5]
- Fraction of total variance explained by each component
- **Values:** [0.9924, 0.0034, 0.0009, 0.0007, 0.0003]
- **Interpretation:**
  - Component 1: 99.24% of total variance
  - Component 2: 0.34% of total variance
  - Component 3: 0.09% of total variance
  - Component 4: 0.07% of total variance
  - Component 5: 0.03% of total variance
  - **TOTAL: 99.76%** ✓

### ✓ explained_variance [5]
- Raw eigenvalues from eigendecomposition
- Absolute magnitude of variation in each direction
- **Values:** [7.92e+06, 2.69e+04, 6.98e+03, 5.33e+03, 2.59e+03]
- Component 1's eigenvalue is ~300× larger than Component 2's

### ✓ config (dictionary)
Settings used during decomposition:
- `n_components: 5`
- `pca_source: SKYCHOPDIFF`
- `normalize: True`
- `scale: False`
- `add_sky_diff: False`
- `noise_cutoff: False`
- `scramble: False`

### ✓ metadata (dictionary)
Information about the input data:
- `n_reference_spectra: 13986`
- `n_channels: 1229`
- `source: SKYCHOPDIFF`
- `mission_id: 2017-02-01_GR_F`

### ✓ spectrum_metadata (list)
Per-spectrum metadata (currently empty, length: 0)

---

## The Mathematical Process

### Step 1: Data Loading and Preprocessing

**Location:** `src/oi_zeigt/pca_analysis/decompose.py`, lines 753-804 in `main_cli()`

```python
# Load SKYCHOPDIFF spectra from FITS
with fits_io.open(fits_file) as hdul:
    matrix_hdu = hdul['MATRIX']
    mask = matrix_hdu.data['OBJECT'] == 'SKYCHOPDIFF'
    spectra = np.array([row['SPECTRUM'] for row in matrix_hdu.data[mask]])

# Result: spectral_matrix [13986, 1229]
# Each row = one SKYCHOPDIFF spectrum
# Each column = one wavelength channel

# Find valid channels and fill NaNs
valid_channels = ~np.all(np.isnan(spectra), axis=0)
first_valid = np.where(valid_channels)[0][0]
last_valid = np.where(valid_channels)[0][-1]
spectra = spectra[:, first_valid:last_valid+1]

for ch in range(spectra.shape[1]):
    channel_data = spectra[:, ch]
    if np.any(np.isnan(channel_data)):
        mean_val = np.nanmean(channel_data)
        spectra[np.isnan(channel_data), ch] = mean_val
```

**Result:** Clean data [13986 × 1229] with no NaNs

### Step 2: Create PCADecomposer and Fit

**Location:** `src/oi_zeigt/pca_analysis/decompose.py`, lines 806-809 in `main_cli()`

```python
decomposer = PCADecomposer(n_components=5, scale=False)
decomposer.fit(spectra)
```

This calls the `PCADecomposer.fit()` method (lines 400-448).

### Step 3: sklearn.PCA Performs the Core Decomposition

**Location:** `src/oi_zeigt/pca_analysis/decompose.py`, lines 420-427 in `PCADecomposer.fit()`

```python
from sklearn.decomposition import PCA

self.pca_model = PCA(
    n_components=5,
    whiten=False,
    random_state=42
)

self.pca_model.fit(spectral_matrix)
```

### Step 4: What sklearn.PCA.fit() Does Internally

#### 4a. **CENTER THE DATA**
```python
mean_spectrum = np.mean(spectra, axis=0)  # Average each channel → [1229]
centered = spectra - mean_spectrum        # Subtract mean from each spectrum → [13986, 1229]
```

**This is your mean_spectrum saved in the pickle!**

#### 4b. **COMPUTE COVARIANCE MATRIX**
```python
covariance = (centered.T @ centered) / (n_spectra - 1)
# Shape: [1229, 1229]
```

Shows correlation between each pair of wavelength channels.

#### 4c. **EIGENDECOMPOSITION (SVD)**
```python
eigenvalues, eigenvectors = np.linalg.svd(centered, full_matrices=False)
# OR equivalently: eigenvalues, eigenvectors = np.linalg.eigh(covariance)
```

For 1229 channels: produces 1229 eigenvectors and 1229 eigenvalues

#### 4d. **SELECT TOP 5 COMPONENTS**
```python
components = eigenvectors[:5]  # Take first 5 (ranked by eigenvalue)
# Shape: [5, 1229]
```

Sorted by eigenvalue (largest first) = ranked by importance

#### 4e. **COMPUTE EXPLAINED VARIANCE RATIO**
```python
explained_variance_ratio = eigenvalues[:5] / sum(all_eigenvalues)
# Result: [0.9924, 0.0034, 0.0009, 0.0007, 0.0003]
# Sum: 0.9976 (99.76%)
```

### Step 5: Extract Results

**Location:** `src/oi_zeigt/pca_analysis/decompose.py`, lines 811-828 in `main_cli()`

```python
mean_spectrum = decomposer.pca_model.mean_
components = decomposer.get_components()
explained_variance = decomposer.get_explained_variance()
explained_variance_ratio = decomposer.get_explained_variance_ratio()

result = DecompositionResult(
    mean_spectrum=mean_spectrum,
    components=components,
    explained_variance=explained_variance,
    explained_variance_ratio=explained_variance_ratio,
    config=DecompositionConfig(n_components=5),
    metadata={
        'n_reference_spectra': 13986,
        'n_channels': 1229,
        'source': 'SKYCHOPDIFF',
        'mission_id': '2017-02-01_GR_F',
    }
)
```

### Step 6: Save to Pickle

**Location:** `src/oi_zeigt/pca_analysis/decompose.py`, lines 830-835 in `main_cli()`

```python
output_file = "decomposition_2017-02-01_GR_F_2017021_components.pkl"
result.save(output_file)
```

Pickle contains all 7 components shown above.

---

## Mathematical Intuition

### What Are Components?

Imagine you have 13,986 points in 1229-dimensional space (one dimension per wavelength channel).

**PCA finds the directions of MAXIMUM SPREAD:**

```
Component 1 Direction: 99.24% of variation
  │
  │     ╱╱╱╱╱╱╱        Most spectra
  │    ╱╱╱╱╱╱╱╱╱╱╱    spread out along
  │   ╱╱╱╱╱╱╱╱╱╱╱╱╱   this axis (overall
  │  ╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱  brightness/
  ├─────────────────→   continuum level)
  │
Component 2 Direction: 0.34% of variation (much smaller spread)
  │
  │ │ │ │ │ │
  ├─│─│─│─│─│─→  Much tighter variation
  │ │ │ │ │ │     (sky lines, residuals)
```

**Components are ORTHOGONAL** (perpendicular) to each other
→ They describe independent types of variation

### Key Mathematical Facts

| Aspect | Definition | Your Values |
|--------|-----------|-------------|
| **Eigenvectors** | Principal directions of variation | components [5, 1229] |
| **Eigenvalues** | Amount of variance in each direction | [7.92e+06, 2.69e+04, ...] |
| **Variance Ratio** | Eigenvalue / Total variance | [99.24%, 0.34%, ...] |
| **Mean** | Center point of data | mean_spectrum [1229] |

---

## Why Component 1 Is So Dominant (99.24%)

All 13,986 sky spectra from the same flight vary mostly in:
- **Overall brightness/continuum level**
- This is the DOMINANT source of variation

**Component 1 captures this.**

**Components 2-5 capture:**
- Minor variations (sky line residuals, instrument artifacts, etc.)
- Much smaller effects

This makes sense: sky observations have similar overall brightness patterns, with much smaller variations in fine details.

---

## Practical Applications

### 1. **COMPRESS**
```python
# Store any spectrum as just 5 numbers instead of 1229
coefficients = np.dot(components, centered_spectrum)  # [5] instead of [1229]
```

### 2. **CORRECT**
```python
# Remove component 1 to subtract sky continuum
reconstructed_without_comp1 = (
    mean_spectrum + np.dot(components[1:].T, coefficients[1:])
)
```

### 3. **ANALYZE**
```python
# Compare spectra by their component coefficients
coefficients = np.dot(components, centered_spectrum)
# Spectra with similar coefficients have similar properties
```

### 4. **FILTER**
```python
# Reconstruct without noisy components 4-5
clean_spectrum = mean + np.dot(components[:3].T, coefficients[:3])
# Keep only the 3 most important components
```

### 5. **UNDERSTAND**
```python
# See which components matter most
print(explained_variance_ratio)  # [99.24%, 0.34%, 0.09%, 0.07%, 0.03%]
# Component 1 dominates! (~300× more important than Component 2)
```

---

## Example: Projecting a New Spectrum

```python
import pickle
import numpy as np

# Load the decomposition
with open('decomposition_2017-02-01_GR_F_2017021_components.pkl', 'rb') as f:
    decomp = pickle.load(f)

# Your new spectrum (1229 channels)
new_spectrum = ...  # shape (1229,)

# Step 1: Center it
centered = new_spectrum - decomp['mean_spectrum']

# Step 2: Project onto components
coefficients = np.dot(decomp['components'], centered)  # shape (5,)

# Step 3: These 5 coefficients describe your spectrum!

# Step 4: Reconstruct (with some error due to variance loss)
reconstructed = decomp['mean_spectrum'] + np.dot(decomp['components'].T, coefficients)
error = new_spectrum - reconstructed
print(f"Reconstruction error: {np.sqrt(np.mean(error**2)):.4f}")
```

---

## Code Architecture

### Class Hierarchy

```
decompose.py
├── DecompositionConfig [dataclass, lines 47-57]
│   └─ Configuration parameters for PCA
│
├── SpectrumPreparator [lines 62-370]
│   └─ Prepare individual spectra for PCA
│
├── PCADecomposer [lines 378-512]  ← CORE PCA CLASS
│   ├─ __init__() [line 384]
│   ├─ fit() [line 400]              ← Calls sklearn.PCA.fit()
│   ├─ get_components() [line 461]
│   ├─ get_explained_variance() [line 475]
│   ├─ get_explained_variance_ratio() [line 482]
│   ├─ transform() [line 495]
│   └─ reconstruct() [line 514]
│
├── DecompositionResult [dataclass, lines 530-570]
│   └─ Container for results
│
└── main_cli() [lines 728-851]  ← COMMAND-LINE INTERFACE
    ├─ Load FITS [lines 753-758]
    ├─ Extract SKYCHOPDIFF [lines 770-781]
    ├─ Preprocess [lines 795-804]
    ├─ Fit PCA [lines 806-809]
    ├─ Extract results [lines 811-828]
    └─ Save pickle [lines 830-835]
```

### Code Flow Diagram

```
main_cli() [line 728]
  │
  ├─→ Load FITS, extract SKYCHOPDIFF spectra [lines 753-758]
  │   Result: spectral_matrix [13986, 1229]
  │
  ├─→ Preprocess: fill NaNs [lines 772-781]
  │   Result: clean [13986, 1229] matrix
  │
  ├─→ Create PCADecomposer(n_components=5) [lines 783-787]
  │   └─→ PCADecomposer.__init__() [line 384]
  │       └─→ self.pca_model = PCA(n_components=5) [line 424]
  │
  ├─→ decomposer.fit(spectral_matrix) [line 786]
  │   └─→ PCADecomposer.fit() [line 400]
  │       └─→ self.pca_model.fit(spectral_matrix) [line 426]
  │           └─→ [SKLEARN MAGIC HAPPENS HERE]
  │               • Center data
  │               • Compute covariance
  │               • Eigendecomposition
  │               • Sort by eigenvalue
  │               • Keep top 5
  │
  ├─→ Extract results [lines 791-799]
  │   ├─→ components = decomposer.get_components() [line 795]
  │   ├─→ explained_variance = decomposer.get_explained_variance() [line 796]
  │   ├─→ explained_variance_ratio = [line 797]
  │   └─→ mean_spectrum = decomposer.pca_model.mean_ [line 794]
  │
  └─→ Save to pickle [lines 809-811]
      └─→ result.save() → Creates .pkl file with all components
```

---

## Key Code Locations

| Component | Location | Lines |
|-----------|----------|-------|
| **PCADecomposer class** | `decompose.py` | 378-512 |
| **PCADecomposer.__init__()** | `decompose.py` | 384-398 |
| **PCADecomposer.fit()** | `decompose.py` | 400-448 |
| **sklearn.PCA integration** | `decompose.py` | 420-427 |
| **PCADecomposer.get_components()** | `decompose.py` | 461-470 |
| **PCADecomposer.get_explained_variance()** | `decompose.py` | 475-488 |
| **PCADecomposer.get_explained_variance_ratio()** | `decompose.py` | 482-494 |
| **DecompositionResult class** | `decompose.py` | 530-570 |
| **main_cli() function** | `decompose.py` | 728-851 |
| **Load FITS** | `decompose.py` | 753-758 |
| **Preprocess** | `decompose.py` | 772-781 |
| **Fit PCA** | `decompose.py` | 806-809 |
| **Extract results** | `decompose.py` | 811-828 |
| **Save pickle** | `decompose.py` | 830-835 |

---

## Summary of Results

**Your decomposition successfully:**
- ✅ Loaded 13,986 SKYCHOPDIFF reference spectra
- ✅ Processed 1229 wavelength channels
- ✅ Extracted 5 principal components
- ✅ Explained 99.76% of total variance
- ✅ Saved decomposition to pickle file

**Component breakdown:**
| Component | Variance Explained |
|-----------|-------------------|
| 1 | 99.24% |
| 2 | 0.34% |
| 3 | 0.09% |
| 4 | 0.07% |
| 5 | 0.03% |
| **TOTAL** | **99.76%** |

---

## Next Steps

With this decomposition, you can:

1. **Apply to Science Data**
   - Project science spectra onto these components
   - Compare their coefficients

2. **Correct Sky Contamination**
   - Remove or reduce component 1 (sky continuum)
   - Improve science signal

3. **Analyze Patterns**
   - Use component coefficients for classification
   - Identify systematic effects

4. **Compress Data**
   - Store spectra as 5-dimensional coefficient vectors
   - Significant data reduction (1229 → 5 values per spectrum)

---

## References

- **Source Code:** `src/oi_zeigt/pca_analysis/decompose.py`
- **Pickle File:** `output/pca_components/decomposition_2017-02-01_GR_F_2017021_components.pkl`
- **Configuration:** `config.toml`
- **sklearn.decomposition.PCA:** [scikit-learn documentation](https://scikit-learn.org/stable/modules/generated/sklearn.decomposition.PCA.html)

---

**Created:** January 19, 2026
**Project:** OI-Zeigt Spectral Analysis Suite
**Version:** 1.0
