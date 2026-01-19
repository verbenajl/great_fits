# PCA Decomposition Tool - Analysis Complete ✅

## 📊 What Was Accomplished

### Analysis Summary
- **Total lines of GILDAS code analyzed**: ~3,500 lines
- **Total documentation created**: 1,740 lines across 5 files
- **Files analyzed**: 6 reference scripts + configs
- **Time investment**: Complete deep-dive analysis

### Files Analyzed
1. ✅ `pca_decompose.py` (462 lines) - Core decomposition algorithm
2. ✅ `pca_correct.py` (1809 lines) - Core correction algorithm  
3. ✅ `pca_utilities.py` (821 lines) - Utilities + line detection
4. ✅ `pca_errors.py` (22 lines) - Exception handling
5. ✅ `M51_pca_reduction.toml` (78 lines) - Parameter config
6. ✅ `mission_id_parameters.yml` (290+ lines) - Flight metadata

### Documentation Created
1. ✅ **PCA_DECOMPOSITION_ANALYSIS.md** (21 KB) - Complete technical deep-dive
2. ✅ **PCA_LINE_DETECTION_CLARIFICATION.md** (11 KB) - Two-mask explanation
3. ✅ **PCA_CORRECTION_CLARIFICATION.md** (8.3 KB) - What PCA corrects for
4. ✅ **PCA_DECOMPOSITION_TOOL_SUMMARY.md** (8.7 KB) - Quick reference + pseudocode
5. ✅ **PCA_ANALYSIS_INDEX.md** (7.8 KB) - Navigation guide

---

## 🎯 Critical Clarifications Made

### Clarification #1: What are we correcting for?
```
ANSWER: Receiver/electronics variations
├─ ✅ Gain fluctuations (amplifier varying)
├─ ✅ Baseline oscillations (electronics artifacts)
├─ ✅ Temperature drifts (receiver sensitivity)
├─ ✅ Phase shifts (signal chain variations)
└─ ❌ NOT atmospheric effects (those aren't receiver variations!)
```

### Clarification #2: Two different masks
```
MASK #1: Fixed Instrumental Artifact (from config)
├─ What: Known receiver problem at specific channels
├─ Why: It's FIXED (constant), PCA learns PATTERNS not constants
├─ Action: Zero it out, then ignore it (we can't fix hardware)

MASK #2: Science Line CII (auto-detected from data)
├─ What: Astronomical feature we're measuring
├─ Why: We want to PRESERVE it, not remove it  
├─ Action: Zero it out during decomposition, protect it during correction
```

### Clarification #3: The workflow
```
DECOMPOSITION: Learn receiver variation patterns
├─ Input: SKYCHOPDIFF reference spectra
├─ Apply BOTH masks (artifact + science)
├─ Run PCA → Learn receiver variations only
└─ Output: PCA model

CORRECTION: Remove receiver variations from science data
├─ Input: M51CENTER science spectra
├─ Load PCA model
├─ Apply learned patterns (with BOTH masks for protection)
└─ Output: Corrected spectra with CII preserved
```

---

## 📋 Key Implementation Insights

### Architecture
- **Two-phase design**: Decomposition → Correction (clean separation)
- **Hierarchical filtering**: Mission → Telescope → Scan → Source
- **Intelligent masking**: Two masks with complementary purposes
- **Flexible configuration**: TOML + YAML, no hardcoding

### Technology Stack
| Component | Current (GILDAS) | Proposed (Python) |
|-----------|---|---|
| FITS I/O | pyclass | astropy.io.fits |
| Data indexing | GILDAS db | pandas DataFrame |
| PCA decomposition | Custom | sklearn.decomposition.PCA |
| Baseline fitting | GILDAS "bas" | scipy.signal |
| Line detection | Custom | opencv (cv2) |
| Config parsing | sicparse | TOML + YAML |

### Configuration Strategy
- **Main config**: `config.toml` with `[pca.common]`, `[pca.decompose]`, `[pca.correct]`
- **Flight metadata**: `mission_id_parameters.yml` with per-flight parameters
- **All configurable**: Line locations, kernel sizes, thresholds, output options

---

## 🚀 Implementation Roadmap

### Phase 1: Foundation (2-3 weeks)
Foundation utilities that other phases depend on
- Error handling
- Configuration system
- FITS indexing
- Baseline fitting
- Masking utilities

### Phase 2: Decomposition (2-3 weeks)
Learning receiver variation patterns
- Data loading & filtering
- Spectrum preparation (apply masks)
- PCA decomposition
- Component smoothing
- Model persistence

### Phase 3: Correction (2 weeks)
Applying learned patterns to science data
- Component selection (noise ratio filtering)
- Spectrum correction (subtraction)
- Diagnostics generation
- FITS output

### Phase 4: Integration (1-2 weeks)
Making it production-ready
- CLI commands
- Workflow automation
- Comprehensive testing
- Documentation

**Total estimated effort**: 4-6 weeks

---

## ✅ Success Criteria

The implementation is complete when:

```
Data I/O:
  ✅ Can load FITS files organized by flight/telescope/scan
  ✅ Can index FITS headers into pandas DataFrame
  ✅ Can filter by mission/telescope/source

Decomposition:
  ✅ Can auto-detect CII line from SKYCHOPDIFF waterfall plot
  ✅ Can apply both masks (artifact + science line)
  ✅ Can run PCA decomposition
  ✅ Can save/load PCA models (pickle)

Correction:
  ✅ Can apply PCA models to M51CENTER spectra
  ✅ Can perform component selection (noise ratio test)
  ✅ Can subtract components safely

Quality Control:
  ✅ CII signal is PRESERVED (not removed)
  ✅ Receiver variations are REDUCED
  ✅ Fixed artifacts are UNCHANGED (accepted)

Output:
  ✅ Can generate diagnostic plots
  ✅ Can write corrected FITS files
  ✅ Can produce detailed logs

Configuration:
  ✅ TOML + YAML configuration system working
  ✅ All parameters externalized
  ✅ Mission metadata properly integrated

Code Quality:
  ✅ Modular design (separate concerns)
  ✅ Comprehensive error handling
  ✅ Well-documented code
  ✅ Unit tests for each component
```

---

## 📚 How to Use the Documentation

### For Quick Start
1. Read **PCA_DECOMPOSITION_TOOL_SUMMARY.md** (15 min)
2. Look at "Implementation Pseudocode" section
3. Start Phase 1 implementation

### For Complete Understanding
1. Read **PCA_DECOMPOSITION_ANALYSIS.md** (45 min)
2. Read **PCA_CORRECTION_CLARIFICATION.md** (20 min)
3. Read **PCA_LINE_DETECTION_CLARIFICATION.md** (20 min)
4. Reference **PCA_DECOMPOSITION_TOOL_SUMMARY.md** while coding

### For Specific Questions
- "How does decomposition work?" → Section 5 of Analysis
- "How does correction work?" → Section 6 of Analysis
- "What gets masked?" → Line Detection Clarification
- "What does PCA correct for?" → Correction Clarification
- "What should I code first?" → Tool Summary pseudocode

---

## 💡 Key Takeaways

1. **This is NOT atmospheric correction**
   - It's receiver/electronics variation removal
   - The name "telluric_line" is historical, means "instrumental artifact"

2. **Two masks, two purposes**
   - Fixed artifact: Masked because it's a constant (PCA learns patterns)
   - Science line: Masked because we want to preserve it

3. **Straightforward two-phase design**
   - Phase 1: Learn variation patterns from reference spectra
   - Phase 2: Apply learned patterns to science spectra

4. **Standard libraries make this easy**
   - sklearn.PCA is well-tested
   - scipy.signal for baseline fitting
   - opencv for line detection
   - astropy for FITS handling

5. **Configuration is flexible**
   - All parameters externalized (TOML + YAML)
   - Per-flight customization possible
   - No hardcoding

---

## 🎉 Status: Ready for Phase 1

You now have:
- ✅ Complete understanding of the algorithm
- ✅ Clear knowledge of masking strategy
- ✅ Detailed pseudocode
- ✅ Configuration structure documented
- ✅ Testing strategy defined
- ✅ 1,740 lines of reference documentation

**All analysis is complete. Ready to start Phase 1 (Foundation) implementation!** 🚀

---

## 📄 Reference Files

All located in `/home/verbena/software/oi_zeigt/`:

1. **PCA_ANALYSIS_INDEX.md** (7.8 KB) - THIS FILE
2. **PCA_DECOMPOSITION_ANALYSIS.md** (21 KB) - Main technical analysis
3. **PCA_LINE_DETECTION_CLARIFICATION.md** (11 KB) - Mask explanations
4. **PCA_CORRECTION_CLARIFICATION.md** (8.3 KB) - What PCA corrects
5. **PCA_DECOMPOSITION_TOOL_SUMMARY.md** (8.7 KB) - Quick reference

---

## Next Action

**Ready to proceed with Phase 1 (Foundation)?**

Phase 1 requires implementing:
- Error classes
- Configuration loader (TOML + YAML)
- FITS indexing (DataFrame creation)
- Baseline fitting utilities
- Masking utilities

All with comprehensive tests! Let me know when you're ready to start. 🚀
