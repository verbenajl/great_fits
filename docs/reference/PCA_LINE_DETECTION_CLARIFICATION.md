# PCA Line Detection - CLARIFICATION

## The Key Distinction

There are **TWO DIFFERENT LINES** that need to be handled in the PCA decomposition/correction workflow:

### 1️⃣ TELLURIC LINE (Instrumental contamination - FIXED)

**What it is:**
- NOT atmospheric absorption, but **instrumental features** from the receivers/electronics
- Instrumental artifacts, gain fluctuations, electronic noise patterns
- Contaminates all observations on the same night/flight
- Examples: Baseline oscillations from electronics, receiver gain variations, amplifier artifacts, etc.

**Where it comes from:**
- **SOURCE**: `mission_id_parameters.yml` 
- Pre-calibrated from previous observations (where this contamination pattern was identified)
- Known and fixed per flight (same electronics/receiver state for all observations)

**Example from mission file:**
```yaml
2017-02-10_GR_F373:
  telluric_line_center: 592        # Peak channel of instrumental feature
  telluric_line_width: 30          # Full width of instrumental feature
```

**How to mask it:**
- Simple calculation: `mask_region = [center - width/2, center + width/2]`
- Mask channels in this range during both decomposition and correction
- Set to zero during baseline fitting to prevent fitting the instrumental artifact

**Why mask it in DECOMPOSITION:**
- Instrumental features would contaminate the PCA components
- We want components that describe the INTRINSIC instrumental VARIATIONS (gain changes, electronics drift)
- NOT the fixed baseline instrumental artifact
- We're trying to learn the VARIABLE part of the instrument, not the fixed part

**Why mask it in CORRECTION:**
- When baseline-fitting M51CENTER spectra, don't try to fit the fixed instrumental artifact
- Polynomial baseline should fit only the true spectral baseline and variable instrument effects
- Instrumental artifact region stays untouched (it's a known artifact we accept)

---

### 2️⃣ SCIENCE LINE (CII - what we care about - AUTO-DETECTED)

**What it is:**
- The actual astronomical feature we're trying to measure
- In M51 case: CII emission line at 158 μm
- The WHOLE POINT of the observation

**Where it comes from:**
- **SOURCE**: AUTO-DETECTION from SKYCHOPDIFF spectra
- Uses waterfall plot + OpenCV image processing
- Computed dynamically per flight

**Why auto-detect?**
```
Question: Why not just put it in config like telluric line?

Answer: 
- Telluric line is INSTRUMENTAL (same night = same position)
- Science line is SPECTRAL FEATURE (depends on exact observation setup)
- If spectrograph was tweaked or calibration changed → science line moves
- Need automatic detection to handle these variations
- Waterfall plot shows concentration of signal → reveals true line position
```

**How detection works:**

```
GOAL: Identify which channels contain the CII science line signal
      in the SKYCHOPDIFF reference spectra

ALGORITHM (from pca_utilities.py::find_lines()):

Step 1: Create waterfall plot
        Stack all SKYCHOPDIFF spectra vertically
        Shape: (n_spectra, n_channels)
        Shows intensity distribution across all spectra

Step 2: Convert to grayscale image
        Normalize intensity range to 0-255
        This is a standard image processing step

Step 3: Apply Gaussian blur
        Kernel size: line_kernel_size (default: 61)
        Purpose: Smooth the image to find contiguous regions
        (Removes pixel-level noise, reveals broader features)

Step 4: Adaptive threshold (iterative)
        Find regions where intensity > mean + cutoff_std * std
        
        First iteration:
        - Calculate mean/std from entire image
        - Create threshold mask
        
        Iterative refinement:
        - Mask detected line regions (set to 0)
        - Recalculate mean/std using only non-line regions
        - Re-threshold with new statistics
        - Repeat until mean/std stabilize
        
        Purpose: Remove false positives from initial threshold

Step 5: Find contours
        Use cv2.findContours() to identify spatial extent
        Result: Precise channel range [ch_start, ch_end] of CII line

OUTPUT: Binary mask showing which channels contain science line
        Shape: (n_channels,)
        Values: 0 (noise) or 1 (science line)
```

**Why this works so well:**
- CII line is BRIGHT in SKYCHOPDIFF (strong emission)
- Baseline regions are UNIFORM (consistent noise level)
- Gaussian blur + threshold naturally separates bright signal from noise
- Iterative refinement handles varying baseline noise

**Why mask it in DECOMPOSITION:**
```
Without masking CII line:
├─ PCA would find a dominant component = the CII line itself
├─ This component would be very strong (it's the main signal)
├─ When we subtract this component during correction:
│  └─ We'd remove the CII line! ❌ (That's the science!)

With masking CII line:
├─ PCA finds components = instrumental/receiver variations ONLY
├─ CII signal is protected from PCA analysis
├─ When we subtract components during correction:
│  └─ We only remove instrumental variations, CII stays intact ✅
```

**Why mask it in CORRECTION:**
```
Even though we masked it during decomposition,
we also mask it during correction to be EXTRA SAFE:

M51CENTER spectrum (what we want to correct):
├─ Contains: Instrumental noise + instrumental artifacts + CII line (signal)
│
├─ Baseline fit (excluding CII): 
│  └─ Polynomial baseline fitted only to non-CII channels
│
├─ Apply PCA corrections (excluding CII):
│  └─ Subtract instrumental/receiver variations
│  └─ But don't touch CII channels
│
└─ Result: Instrumental variations removed, CII line intact ✅
```

---

## Combined Workflow

### During DECOMPOSITION:

```
Load SKYCHOPDIFF spectra for flight F373
    ↓
For each spectrum:
    │
    ├─ Get instrumental artifact location from config:
    │  └─ mask_instrumental = [592-15, 592+15] = [577, 607]
    │
    ├─ Auto-detect science line (CII) from waterfall:
    │  └─ mask_science = detect_cii_line(all_spectra)
    │     └─ Example result: [450, 550]
    │
    ├─ Create combined mask:
    │  └─ mask_combined = mask_instrumental ∪ mask_science
    │                   = [450, 550] ∪ [577, 607]
    │
    ├─ Zero out masked regions:
    │  └─ spectrum[450:550] = 0
    │  └─ spectrum[577:607] = 0
    │
    ├─ Baseline fit (only in non-masked regions):
    │  └─ Fit polynomial using channels outside both masks
    │
    ├─ Subtract baseline:
    │  └─ spectrum -= baseline
    │
    └─ Append to input matrix for PCA
        (spectrum has: receiver variations, electronics drifts, noise, but NO CII!)

PCA.fit(input_matrix)
    ↓
Extract components (these are ONLY instrumental/receiver variations)
```

### During CORRECTION:

```
Load M51CENTER spectra for same flight F373
    ↓
Use same masks: instrumental artifact + science line (auto-detected earlier)
    ↓
For each M51CENTER spectrum:
    │
    ├─ Baseline fit (excluding both mask regions)
    │  └─ Polynomial fitted to clean channels only
    │
    ├─ Subtract baseline
    │
    ├─ For each PCA component:
    │  │
    │  ├─ Calculate fit coefficient (least-squares)
    │  │  BUT exclude masked channels from fitting
    │  │  (Don't let CII or instrumental artifact affect the fit)
    │  │
    │  ├─ Calculate noise ratio (only for non-masked channels)
    │  │
    │  ├─ If noise_ratio < threshold:
    │  │  └─ Subtract component (but NOT from masked channels!)
    │  │
    │  └─ Else:
    │     └─ Skip component
    │
    └─ Return corrected spectrum
        (CII and instrumental artifact regions untouched, receiver variations removed)
```

---

## Configuration Summary

### telluric_line parameters (from `mission_id_parameters.yml`):
```yaml
2017-02-10_GR_F373:
  telluric_line_center: 592         # HARDCODED (always use this)
  telluric_line_width: 30           # HARDCODED (always use this)
```

### science_line parameters (to add to `M51_pca_reduction.toml`):
```toml
[common]
# Auto-detect science line (CII) from waterfall plot
enable_science_line_detection = true

# Settings for auto-detection
line_kernel_size = 61               # Gaussian blur kernel
cutoff_std = 2                      # Threshold: mean + N*std
max_iterations = 10                 # Iterative refinement iterations
```

---

## Python Implementation Pseudocode

```python
class PCADecomposition:
    def __init__(self, mission_id, telescope, config, fits_files):
        self.config = config
        self.mission_id = mission_id
        self.telescope = telescope
        
        # Load telluric line from mission config (ALWAYS)
        self.telluric_mask = self._create_telluric_mask()
        
        # Auto-detect science line (CII) from SKYCHOPDIFF
        self.science_mask = self._auto_detect_science_line(fits_files)
        
        # Combine masks
        self.combined_mask = self.telluric_mask | self.science_mask
        
    def _create_telluric_mask(self):
        """Create mask from mission_id_parameters.yml"""
        center = mission_params[self.mission_id]['telluric_line_center']
        width = mission_params[self.mission_id]['telluric_line_width']
        mask = np.zeros(n_channels, dtype=bool)
        mask[center - width//2 : center + width//2] = True
        return mask
    
    def _auto_detect_science_line(self, fits_files):
        """Detect CII line from SKYCHOPDIFF waterfall plot"""
        # Load all SKYCHOPDIFF spectra
        spectra = [load_spectrum(f) for f in fits_files if f.source == 'SKYCHOPDIFF']
        
        # Create waterfall plot
        waterfall = np.vstack(spectra)  # Shape: (n_spectra, n_channels)
        
        # Use find_lines() from pca_utilities
        threshold, _ = find_lines(
            input_data=waterfall,
            kernel_size=self.config['line_kernel_size'],
            cutoff_std=self.config['cutoff_std']
        )
        
        # Extract channel ranges from threshold
        science_mask = threshold > 0
        return science_mask
    
    def prepare_spectrum(self, spectrum):
        """Prepare spectrum for PCA (mask both line types)"""
        # Zero out masked regions
        spectrum = spectrum.copy()
        spectrum[self.combined_mask] = 0
        
        # Baseline fit (only non-masked channels)
        clean_channels = ~self.combined_mask
        baseline = fit_baseline(spectrum[clean_channels])
        
        # Subtract baseline
        spectrum -= baseline
        
        return spectrum
```

---

## Key Takeaway

✅ **Instrumental artifact mask**: From config (`mission_id_parameters.yml`) - always excluded
✅ **Science line (CII)**: Auto-detected from data (waterfall + opencv) - always excluded
✅ **Both excluded**: During decomposition AND correction
✅ **Result**: PCA learns only instrumental/receiver variations, preserves CII signal
