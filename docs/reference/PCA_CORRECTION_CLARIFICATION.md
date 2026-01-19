# CRITICAL CORRECTION: What PCA Actually Corrects For

## The Essential Distinction

**You clarified**: "We don't want to correct for atmospheric features, but rather instrumental features due to the electronics in the receivers, or gain changes, and such."

This is **fundamentally important** because it changes what we're doing:

### NOT: Correcting for the Earth's atmosphere
### YES: Correcting for receiver/electronics instrumental variations

---

## What the PCA Actually Learns & Removes

### The PCA Model captures:

**Receiver variations** (what PCA DOES learn):
- Gain fluctuations: Receiver amplifier gain varying spectrum-to-spectrum
- Baseline oscillations: Electronic artifacts in the receiver backend
- Temperature effects: Receiver sensitivity drifting with temperature
- Phase shifts: Electronic phase variations in the signal chain
- Noise distribution changes: How receiver noise varies over time

**Examples**:
```
Component 1: "Receiver gain increasing over time"
  └─ All channels scale up uniformly

Component 2: "Baseline oscillation pattern"
  └─ Regular ripples in baseline due to electronics

Component 3: "Noise getting worse in channels 500-600"
  └─ Receiver sensitivity dropping in specific frequency range

Component 4: "Gain imbalance between polarizations"
  └─ Left/right receiver behaving differently
```

### What the PCA model does NOT include:

```
❌ Science line (CII):    MASKED before PCA (auto-detected from waterfall)
❌ Instrumental artifact:  MASKED before PCA (from mission_id_parameters.yml)
❌ Atmospheric effects:    Not relevant here (they're not in receiver variations!)
```

---

## The Fixed Instrumental Artifact (Masked in config)

This is the **telluric_line_center/width** in `mission_id_parameters.yml`:

**What it is:**
- A KNOWN, FIXED instrumental problem with the receiver
- Always appears at the same channels for a given night
- Example: "Receiver gain is 20% lower at channels 577-607 on this date"
- This is a KNOWN ARTIFACT that we accept and ignore

**Why we mask it:**
- It's FIXED, not variable → doesn't belong in PCA
- PCA is for learning VARIABLE patterns, not fixed offsets
- If we included it in PCA, we'd be asking PCA to learn a constant offset
- We want PCA to learn VARIATIONS, not constants

**Example scenario:**
```
All observations on 2017-02-10_GR_F373:
├─ Receiver has permanent 20% gain loss at channels 577-607
├─ This is NOT something we correct (it's a hardware issue)
├─ But we MASK it during PCA training
│  └─ Reason: We want PCA to learn how gain varies between scans,
│             not the fixed gain loss that's always there
├─ During correction: Leave this region untouched
│  └─ We can't fix hardware artifacts with PCA anyway
```

---

## The Science Line (Masked via auto-detection)

This is the **CII line we care about**:

**What it is:**
- The astronomical feature: CII emission from M51
- The SIGNAL we're trying to measure
- Strong and consistent across observations

**Why we mask it:**
- If we DON'T mask it: PCA learns "CII line pattern" as a principal component
  - Then when we apply PCA to M51CENTER, we subtract the CII!
  - We'd remove our own target signal! ❌
- If we DO mask it: PCA learns only receiver variations
  - When we apply PCA, we keep the CII intact ✅

**Example:**
```
SKYCHOPDIFF spectrum (reference, no science):
├─ Receiver variations: "gain is high in this scan"
├─ Fixed artifact: [577-607] channels
├─ CII line: Appears here too! (in the reference spectrum)

Without masking CII:
├─ PCA component: "CII line appears strongly"
├─ In correction: Subtract CII component from M51CENTER
├─ Result: CII signal removed! ❌ (disaster!)

With masking CII:
├─ PCA component: "gain is high in this scan"
├─ In correction: Subtract gain variation from M51CENTER
├─ Result: Receiver variations corrected, CII intact ✅
```

---

## The Complete Picture

### What gets removed during correction:

```
M51CENTER spectrum (before correction):
├─ Receiver gain variations ────────────┐
├─ Baseline oscillations ──────────────┤─ THESE GET REMOVED ✓
├─ Temperature drifts in receiver ────┤  (via PCA subtraction)
├─ Other receiver artifacts ──────────┘
├─ CII line (science signal) ──────────── STAYS INTACT ✓ (masked)
├─ Fixed instrumental artifact ───────── UNCHANGED (masked, can't fix)
└─ Noise ────────────────────────────────  UNCHANGED (PCA doesn't reduce noise)

After correction:
├─ CII line (clean, without receiver variations!)
├─ Fixed instrumental artifact (still there, but we accept it)
└─ Noise (same as before)
```

### Configuration meaning:

```yaml
# mission_id_parameters.yml
2017-02-10_GR_F373:
  telluric_line_center: 592
  telluric_line_width: 30

# INTERPRETATION:
# "On this flight, there's a fixed receiver artifact at channels 577-607.
#  We know about it, we accept it, and we don't try to correct it.
#  We'll mask it during PCA training to avoid confusing PCA with a constant."
```

---

## Implementation Priority

### What MUST be in config (fixed instrumental artifacts):

```toml
# mission_id_parameters.yml
2017-02-10_GR_F373:
  telluric_line_center: 592         # Channel of fixed artifact
  telluric_line_width: 30            # Width of fixed artifact
```

**Why**:
- These are KNOWN problems specific to each observing run
- They're documented in the mission metadata
- They're instrument-specific, not data-dependent

### What MUST be auto-detected (science line):

```python
# In pca_decompose.py
science_line_mask = auto_detect_from_waterfall(skychopdiff_spectra)
```

**Why**:
- CII line position might vary slightly with different setups
- We need to detect the ACTUAL signal location in the data
- Waterfall plot shows exactly where the signal concentrates

---

## Final Clarification Summary

| Item | Source | Type | Purpose |
|------|--------|------|---------|
| **Receiver gain variations** | Inherent in spectra | VARIABLE | Learn & remove via PCA |
| **Baseline oscillations** | Inherent in spectra | VARIABLE | Learn & remove via PCA |
| **Temperature drifts** | Inherent in spectra | VARIABLE | Learn & remove via PCA |
| **Fixed artifact (channels 577-607)** | mission_id_parameters.yml | FIXED | Mask & ignore (can't fix) |
| **CII science line** | Auto-detect from waterfall | VARIABLE | Mask & preserve (don't remove!) |

### PCA removes:
✓ Receiver gain fluctuations
✓ Electronic baseline variations
✓ Instrumental drifts
✓ Temperature-dependent variations
✓ Anything that VARIES and contaminates the science

### PCA does NOT remove:
✗ Fixed receiver artifacts (they're constants, not components)
✗ Atmospheric effects (not part of receiver variations!)
✗ Science line (we masked it!)
✗ Random noise (PCA handles correlated variations, not random noise)

---

## Example Correction Scenario

**Flight**: 2017-02-10_GR_F373
**Receiver problem**: Gain fluctuates ±5% between scans

```
Scan 1 (M51CENTER):  CII + 1.0×receiver_gain + fixed_artifact + noise
Scan 2 (M51CENTER):  CII + 1.05×receiver_gain + fixed_artifact + noise
Scan 3 (M51CENTER):  CII + 0.95×receiver_gain + fixed_artifact + noise

SKYCHOPDIFF (reference):
Scan 1:  reference + 1.0×receiver_gain + fixed_artifact + noise
Scan 2:  reference + 1.05×receiver_gain + fixed_artifact + noise
Scan 3:  reference + 0.95×receiver_gain + fixed_artifact + noise

PCA learned:
├─ Component 1: "gain variation pattern" (detects the ±5% fluctuation)
├─ NOT including: Fixed artifact (masked)
├─ NOT including: CII line (masked)

After correction:
├─ Scan 1: CII + 1.0×receiver_gain (cleaned!)
├─ Scan 2: CII + 1.05×receiver_gain (cleaned!)
├─ Scan 3: CII + 0.95×receiver_gain (cleaned!)
└─ All scans now have receiver gain variations removed ✓
```

---

## Bottom Line

🎯 **PCA corrects for**: Receiver/electronics variations (gain changes, drift, oscillations)
🎯 **PCA does NOT correct for**: Atmospheric effects
🎯 **We protect from PCA**: Science line (CII) via auto-detection
🎯 **We ignore in PCA**: Fixed artifacts (they're constants) via config masking
🎯 **Result**: Clean spectra with receiver variations removed, science signal intact
