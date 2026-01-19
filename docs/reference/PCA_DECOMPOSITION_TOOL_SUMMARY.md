# PCA Decomposition Tool - Final Clarification Summary

## ✅ CORRECTED UNDERSTANDING

### What PCA Actually Does

**Corrects for:**
- ✅ Receiver gain fluctuations (amplifier varying over time)
- ✅ Baseline oscillations (electronics artifacts)
- ✅ Temperature-dependent receiver variations
- ✅ Phase shifts in signal chain
- ✅ Gain imbalances and drifts

**Does NOT correct for:**
- ❌ Atmospheric effects (not part of receiver/electronics)
- ❌ Fixed hardware artifacts (they're constants, not learnable patterns)
- ❌ Random noise (PCA handles correlated patterns, not noise)
- ❌ Science line (intentionally masked to preserve it)

---

## Two Masks, Two Purposes

### Mask #1: Fixed Instrumental Artifact (from config)

**From**: `mission_id_parameters.yml` - `telluric_line_center/width`

**Why called "telluric"**: Historically named but actually means "instrumental artifact"

**What it represents**: A KNOWN receiver problem at specific channels
- Example: "Receiver has 20% lower gain at channels 577-607 on this date"

**Why mask it**:
- It's FIXED (same for all observations that night)
- PCA learns PATTERNS, not constants
- Masking prevents PCA from learning a constant offset
- We accept this artifact (hardware issue, can't fix)

**In practice**:
```yaml
2017-02-10_GR_F373:
  telluric_line_center: 592
  telluric_line_width: 30
  # Interpretation: Mask channels [577-607] - known artifact zone
```

---

### Mask #2: Science Line (CII) - Auto-detected

**From**: Waterfall plot + OpenCV (automatic detection from SKYCHOPDIFF)

**What it represents**: The astronomical feature we're measuring (CII emission)

**Why mask it**:
- We want to PRESERVE it, not remove it
- If we don't mask it: PCA learns "CII pattern" → removes it during correction ❌
- If we do mask it: PCA learns receiver variations only → keeps CII ✅

**In practice**:
```python
# Auto-detect from SKYCHOPDIFF waterfall
science_line_mask = find_lines(waterfall, kernel_size=61)
# Mask channels ~[450-550] or whatever the CII line is at
```

---

## Complete Workflow

### DECOMPOSITION Phase:

```
Goal: Learn receiver variation patterns from SKYCHOPDIFF reference spectra

Steps:
1. Load all SKYCHOPDIFF spectra for a flight/telescope pair
2. Apply Mask #1 (fixed artifact): Zero out channels [577-607]
3. Apply Mask #2 (science line): Zero out channels [450-550]
4. Baseline fit in non-masked regions
5. Run PCA → Learn receiver variation components
6. Save PCA model

Result: Components that describe:
  ✓ How receiver gain varies between observations
  ✓ How baseline oscillates
  ✓ How noise characteristics change
  ✓ Other receiver variations
  
Does NOT include:
  ✗ The fixed artifact
  ✗ The CII line
```

### CORRECTION Phase:

```
Goal: Apply learned receiver patterns to M51CENTER science spectra

Steps:
1. Load M51CENTER spectra for same flight/telescope
2. Load saved PCA model (has learned components)
3. For each spectrum:
   a. Baseline fit (excluding both masks)
   b. For each PCA component:
      - Fit component to spectrum (least-squares)
      - Check if fit is good (noise ratio test)
      - If good: subtract fitted component
      - If bad: skip (too noisy)

Result: M51CENTER spectra with:
  ✓ Receiver variations removed
  ✓ CII signal preserved
  ✓ Fixed artifact untouched (we can't fix hardware issues)
```

---

## Configuration Changes Needed

### Add to main `config.toml`:

```toml
[pca.common]
enabled = true
pca_source = "SKYCHOPDIFF"
smoothing_kernel_size = 3
rolling_noise_window = 11
decomposition_type = "PCA"  # or "SparsePCA", "ICA"
do_scale = false

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
```

### Keep `mission_id_parameters.yml` as-is:

```yaml
2017-02-10_GR_F373:
  telluric_line_center: 592       # Fixed artifact location
  telluric_line_width: 30         # Fixed artifact width
  drop:                           # Data quality filters
    telescope: ["LFAH_3"]         # Exclude bad telescopes
    scans:
      complete: [17563, 17604]    # Exclude bad scans
```

---

## Implementation Pseudocode

```python
class PCADecomposition:
    """Learn receiver variation patterns from SKYCHOPDIFF reference spectra."""
    
    def __init__(self, mission_id, telescope, config, fits_dir):
        # Load config
        self.config = config
        
        # Load mission metadata
        self.mission_params = load_yaml('mission_id_parameters.yml')
        
        # Create Mask #1: Fixed instrumental artifact
        center = self.mission_params[mission_id]['telluric_line_center']
        width = self.mission_params[mission_id]['telluric_line_width']
        self.artifact_mask = create_mask(center, width, n_channels=1024)
        
        # Load SKYCHOPDIFF spectra
        skychopdiff = load_spectra(fits_dir, mission_id, telescope, source='SKYCHOPDIFF')
        
        # Create Mask #2: Science line (auto-detected from waterfall)
        waterfall = np.vstack([s.flux for s in skychopdiff])
        self.science_mask = auto_detect_line(waterfall, kernel_size=61)
        
        # Combined mask
        self.combined_mask = self.artifact_mask | self.science_mask
        
        # Prepare spectra for PCA
        prepared_spectra = []
        for spectrum in skychopdiff:
            spec = spectrum.flux.copy()
            spec[self.combined_mask] = 0      # Zero masked regions
            baseline = fit_baseline(spec, ~self.combined_mask)
            spec -= baseline                   # Subtract baseline
            prepared_spectra.append(spec)
        
        # Run PCA decomposition
        self.pca_model = PCA(n_components=5)
        self.pca_model.fit(np.array(prepared_spectra))
        
        # Save model
        save_pickle(self.pca_model, f'{mission_id}_{telescope}.pkl')


class PCACorrection:
    """Apply learned receiver patterns to M51CENTER science spectra."""
    
    def __init__(self, mission_id, telescope, pca_model, masks):
        self.pca_model = pca_model
        self.artifact_mask = masks['artifact']
        self.science_mask = masks['science']
        self.combined_mask = masks['combined']
    
    def correct_spectrum(self, spectrum):
        """Remove receiver variations from a spectrum."""
        spec = spectrum.copy()
        
        # Baseline fit (excluding masked regions)
        baseline = fit_baseline(spec, ~self.combined_mask)
        spec -= baseline
        
        # Apply PCA corrections
        for i, component in enumerate(self.pca_model.components_):
            # Fit component using only non-masked channels
            coeff, noise_ratio = fit_component(spec, component, ~self.combined_mask)
            
            # Check if component is good to use
            if noise_ratio < self.noise_threshold:
                spec -= coeff * component
        
        return spec
```

---

## Key Points for Implementation

1. **Fixed artifact mask** (from config):
   - Always applied, never changes for a given flight
   - Prevents PCA from learning constant offsets
   - Kept untouched during correction

2. **Science line mask** (auto-detected):
   - Computed once per flight/telescope pair
   - Protects CII signal from being learned/removed by PCA
   - Same mask used in decomposition AND correction

3. **Component selection** (noise ratio based):
   - Not all components are good to use
   - Noisy components might amplify noise instead of correcting
   - Noise ratio threshold filters out unreliable components

4. **Output**:
   - Corrected spectra with receiver variations removed
   - CII signal intact
   - Fixed artifacts unchanged (accepted as hardware issue)

---

## Testing Strategy

Once implemented, verify with:

```python
# Load corrected spectrum
corrected = load_corrected_spectrum()

# Check 1: CII line is still there
assert np.max(corrected[450:550]) > threshold, "CII should be present"

# Check 2: Baseline is smoother than original
assert smoothness(corrected) > smoothness(original), "Baseline should be smoother"

# Check 3: Fixed artifact region is untouched
assert np.allclose(corrected[577:607], original[577:607]), "Artifact region unchanged"

# Check 4: Variance reduced
assert np.std(corrected) < np.std(original), "Overall variation should decrease"
```

---

## Files Created for Reference

1. **PCA_DECOMPOSITION_ANALYSIS.md** - Complete technical analysis
2. **PCA_LINE_DETECTION_CLARIFICATION.md** - Two types of lines explained
3. **PCA_CORRECTION_CLARIFICATION.md** - What PCA actually corrects for
4. **PCA_DECOMPOSITION_TOOL_SUMMARY.md** - This file

All in `/home/verbena/software/oi_zeigt/`

Ready to start implementation! 🚀
