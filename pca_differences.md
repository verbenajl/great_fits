# PCA Implementation Differences: Original pyclass vs New FITS Code

## Status legend
- [ ] Not yet addressed
- [~] Partially addressed / under investigation
- [x] Fixed

---

## 1. [x] Normalisation before PCA decomposition

**Original**: Scaler is commented out in `InputSpectra.normalized_array` — raw spectra passed to `PCA.fit()` with no normalisation.

**New code** (`decompose.py` `SpectrumPreparator.prepare`): Each SKYCHOPDIFF spectrum is normalised to zero mean / unit std before PCA.

**Impact**: Components in the new code capture variance patterns of normalised spectra. In the original, components reflect raw-amplitude patterns. This is likely the largest single difference affecting output quality.

**Files**: `src/oi_zeigt/pca_analysis/decompose.py` (~line 247)

---

## 2. [x] Component smoothing

**Original** (`pca_decompose.py` `smooth_pca_components`): Components are smoothed with a boxcar kernel **at decomposition time** and stored already-smoothed. The subtracted pattern is always smooth.

**New code**: Components are NOT smoothed before subtraction. The `smoothing_kernel_size` parameter only affects the `spectrum_std` used in the noise ratio calculation — NOT what gets subtracted.

**Impact**: The original produces a smoother correction; the new code subtracts sharp, unsmoothed components which could introduce spectral artefacts.

**Files**: `src/oi_zeigt/pca_analysis/decompose.py` (missing smoothing step), `src/oi_zeigt/pca_analysis/pca_correct_fits.py` (~line 402)

---

## 3. [x] Component std in noise ratio

**Original** (`derive_noise_ratio_and_windows` line 928):
```python
noise_in_scaled_component = scaled_component.std()   # ALL channels
```

**New code** (`_pca_apply_correction` line 401):
```python
component_std = np.nanstd(coeff[i] * comp_masked)   # good_channels only
```

**Impact**: When a bright science line is present, the original std includes line channels and is larger → lower noise ratio → more components applied. The new code restricts to baseline channels → different ratio.

**Files**: `src/oi_zeigt/pca_analysis/pca_correct_fits.py` (~line 401)

---

## 4. [x] Noise ratio cutoff mechanism

**Original**: Complex KDE (Kernel Density Estimate) of noise ratio distribution across all spectra per component. Fits a Gaussian to the first peak. Cutoff = `center + 3*sigma`. Also has `global_noise_ratio_cutoff` which rejects entire component if distribution peak exceeds a value. Uses `noise_ratio_scaled` (smoothed spectrum std).

**New code**: Simple per-spectrum threshold — skip if `noise_ratio > cutoff_noise_ratio`. Distribution-based auto-cutoff completely absent.

**Impact**: The original adapts to the actual distribution of noise ratios in the data. The new code uses a fixed threshold regardless of data properties.

**Files**: `src/oi_zeigt/pca_analysis/pca_correct_fits.py` (~line 432), `src/oi_zeigt/pca_analysis/pca_utilities.py` (`derive_noise_ratio_cutoff`, `should_component_be_used`)

---

## 5. [x] Missing `cut_coefficients` parameter

**Original** (line 1422–1425): Component skipped if `|coeff| < cut_coefficients`. Guards against subtracting components that barely project onto the spectrum.

**New code**: Not implemented.

**Impact**: Components with negligible coefficients may be subtracted unnecessarily.

**Files**: `src/oi_zeigt/pca_analysis/pca_correct_fits.py` (missing)

---

## 6. [x] Coefficient re-fitting after component selection

**Original** (lines 1505–1514): After filtering to the components that pass all cutoffs, coefficients are **re-fitted using only those selected components**:
```python
coeff = ma.dot(components_used_for_correction_list[:, good_channels],
               scaled_fit_spec[good_channels])
```

**New code**: No re-fitting. Original coefficient from full fit used even when some components are excluded.

**Impact**: The re-fitted coefficients give a more accurate projection onto the subset actually being subtracted.

**Files**: `src/oi_zeigt/pca_analysis/pca_correct_fits.py` (~line 437)

---

## 7. [ ] StandardScaler during correction

**Original**: Spectrum scaled with `sklearn StandardScaler` before fitting, `inverse_transform` after. Default `do_scale=True`.

**New code**: No scaling during correction.

**Impact**: Scaling normalises spectrum amplitude before projecting onto components. Without scaling, spectra with very different overall amplitudes may produce inconsistent coefficients.

**Files**: `src/oi_zeigt/pca_analysis/pca_correct_fits.py` (missing)

---

## 8. [ ] n_components off-by-one

**Original**: `PCA(self.number_components_used - 1)` — one fewer component than requested.

**New code**: `PCA(n_components=self.n_components)` — exact number requested.

**Impact**: Minor, but results in one extra component compared to original for the same config value.

**Files**: `src/oi_zeigt/pca_analysis/decompose.py` (~line 492)

---

## 9. [x] Baseline subtraction in decomposition

**Original** (`prepare_spectrum`): Order-0 baseline (mean subtraction) applied to SKYCHOPDIFF spectra before PCA.

**New code** (`SpectrumPreparator.prepare`): Polynomial baseline of configurable order (default changed to 0 to match original). Still configurable via `pca.common.baseline_order` in config.

**Files**: `src/oi_zeigt/pca_analysis/decompose.py` (~line 116)

---

## 10. [x] Variance cutoff: inconsistency between two code paths

**Original**: `<= cutoff` (skip component if variance ratio ≤ cutoff).

**Fixed** in `_pca_apply_correction` (line 389): `<= cutoff_variance` ✓ — this is the CLI path.

**Still wrong** in `PCACorrector.apply_correction` (line 981): `< cutoff_variance` ✗ — strict less-than, off by one edge case.

**Impact**: Minor edge case only (only matters when var_ratio == cutoff exactly), but worth fixing for consistency.

**Files**: `src/oi_zeigt/pca_analysis/pca_correct_fits.py` (~line 981)
