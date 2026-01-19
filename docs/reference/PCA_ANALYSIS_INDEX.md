# PCA Decomposition Tool - Complete Analysis Index

## 📄 Documentation Files

### 1. **PCA_DECOMPOSITION_ANALYSIS.md** (21 KB)
**Complete technical analysis of the GILDAS/CLASS reference code**

Contents:
- Executive summary of two-phase PCA approach
- Complete architecture & workflow
- Configuration system (main config + mission-specific parameters)
- Data filtering hierarchy (flight → telescope → scan → source)
- **Line detection & masking strategy** (with code examples)
- PCA decomposition process (sklearn details)
- PCA correction process (component selection & application)
- Error handling & edge cases
- Implementation decisions for Python port
- Configuration migration plan (TOML + YAML)
- Complete implementation roadmap
- Dependencies list
- Data flow diagram

**Read this for**: Complete understanding of what the tool does and how it works

---

### 2. **PCA_LINE_DETECTION_CLARIFICATION.md** (11 KB)
**Detailed explanation of the two types of lines and masking strategies**

Contents:
- **Section 1: Fixed Instrumental Artifact**
  - What it is (receiver/electronics issues)
  - Where it comes from (mission_id_parameters.yml)
  - How to mask it (simple channel range)
  - Why mask it in decomposition & correction
  - Configuration example

- **Section 2: Science Line (CII)**
  - What it is (astronomical feature we measure)
  - Auto-detection from waterfall plot
  - 5-step detection algorithm (blur → threshold → iterate → contours)
  - Why auto-detect instead of config
  - Why mask it in decomposition & correction

- Complete combined workflow (decomposition + correction phases)
- Configuration parameters for both lines
- Python implementation pseudocode
- Key takeaway summary

**Read this for**: Understanding why we have two masks and what each one does

---

### 3. **PCA_CORRECTION_CLARIFICATION.md** (8.3 KB)
**What PCA actually corrects for (CRITICAL CLARIFICATION)**

Contents:
- ✅ What PCA learns & removes (receiver variations)
  - Gain fluctuations
  - Baseline oscillations
  - Temperature effects
  - Electronics drift
  
- ❌ What PCA does NOT learn (and why)
  - Atmospheric effects (not receiver variations!)
  - Fixed hardware artifacts (they're constants)
  - Science line (intentionally masked)
  - Random noise (not a pattern)

- The fixed artifact explained
- The science line protection
- Complete picture (before/after correction)
- Implementation priority
- Final clarification summary (table)
- Example correction scenario
- Bottom line summary

**Read this for**: Understanding that we're correcting RECEIVER VARIATIONS, not atmospheric effects

---

### 4. **PCA_DECOMPOSITION_TOOL_SUMMARY.md** (8.7 KB)
**Quick reference guide with implementation details**

Contents:
- Corrected understanding (what PCA does)
- Two masks explained (with purposes)
- Complete workflow (decomposition → correction)
- Configuration changes needed (main config.toml)
- Keep mission_id_parameters.yml as-is
- Implementation pseudocode (2 classes with complete logic)
- Key points for implementation (4 critical details)
- Testing strategy (verification checks)
- Files created for reference

**Read this for**: Quick reference before implementation + pseudocode structure

---

## 🎯 Quick Navigation by Question

### "What does the tool do?"
→ Read **PCA_DECOMPOSITION_ANALYSIS.md** sections 1-3

### "What are we correcting for?"
→ Read **PCA_CORRECTION_CLARIFICATION.md** (entire)

### "What are the two masks?"
→ Read **PCA_LINE_DETECTION_CLARIFICATION.md** sections 1-2

### "How do decomposition & correction work?"
→ Read **PCA_DECOMPOSITION_ANALYSIS.md** sections 5-6

### "How should I implement this?"
→ Read **PCA_DECOMPOSITION_TOOL_SUMMARY.md** sections 3-4

### "What configuration changes are needed?"
→ Read **PCA_DECOMPOSITION_TOOL_SUMMARY.md** section "Configuration Changes Needed"

### "How do I test it?"
→ Read **PCA_DECOMPOSITION_TOOL_SUMMARY.md** section "Testing Strategy"

---

## 📋 Key Concepts Summary

### The Masks

| Aspect | Fixed Artifact | Science Line (CII) |
|--------|---|---|
| **Where** | mission_id_parameters.yml | Auto-detected from waterfall |
| **What it is** | Receiver/electronics problem at specific channels | Astronomical feature we're measuring |
| **Why mask** | It's fixed (constant), PCA learns patterns not constants | We want to preserve it, not remove it |
| **In decomposition** | Zero it out (don't learn the constant) | Zero it out (don't learn the signal) |
| **In correction** | Leave it alone (we accept the artifact) | Don't touch it (protect the science) |

### The PCA Model

**Learns:**
- ✅ Receiver gain fluctuations
- ✅ Baseline oscillations
- ✅ Temperature drifts
- ✅ Electronics variations

**Does NOT learn:**
- ❌ Fixed artifacts (masked)
- ❌ Science line (masked)
- ❌ Atmospheric effects (not in receiver variations!)

### The Correction

**Removes:**
- ✅ Receiver gain variations
- ✅ Baseline oscillations
- ✅ Electronics drifts

**Preserves:**
- ✅ Science signal (CII) - was masked
- ✅ Fixed artifacts - we accept them

---

## 🚀 Implementation Order

1. **Read** all 4 documents in this order:
   1. PCA_DECOMPOSITION_ANALYSIS.md (understand the problem)
   2. PCA_CORRECTION_CLARIFICATION.md (critical clarification)
   3. PCA_LINE_DETECTION_CLARIFICATION.md (understand the masks)
   4. PCA_DECOMPOSITION_TOOL_SUMMARY.md (quick reference)

2. **Design** the module structure:
   - `pca_analysis/decompose.py` - Decomposition logic
   - `pca_analysis/correct.py` - Correction logic
   - `pca_analysis/utilities.py` - Helper functions

3. **Implement** in phases:
   - Phase 1: Foundation (errors, config, utilities)
   - Phase 2: Decomposition
   - Phase 3: Correction
   - Phase 4: Integration

4. **Test** with:
   - Unit tests for each component
   - Integration tests with real FITS files
   - Verification that science line is preserved
   - Verification that receiver variations are reduced

---

## 📊 File Statistics

| File | Size | Lines | Purpose |
|------|------|-------|---------|
| PCA_DECOMPOSITION_ANALYSIS.md | 21 KB | 642 | Complete technical analysis |
| PCA_LINE_DETECTION_CLARIFICATION.md | 11 KB | 320 | Two masks explained |
| PCA_CORRECTION_CLARIFICATION.md | 8.3 KB | 240 | What PCA corrects for |
| PCA_DECOMPOSITION_TOOL_SUMMARY.md | 8.7 KB | 310 | Quick reference + pseudocode |

**Total: 49 KB of documentation**

---

## ✅ Analysis Complete!

All aspects of the PCA decomposition tool have been analyzed and documented:

- ✅ Architecture (two-phase: decompose → correct)
- ✅ Configuration (main + mission-specific)
- ✅ Data filtering (hierarchical by mission/telescope/scan)
- ✅ Line masking (two masks, two purposes)
- ✅ PCA algorithm (sklearn implementation)
- ✅ Correction process (component selection + subtraction)
- ✅ Error handling (edge cases documented)
- ✅ Implementation strategy (Python port approach)
- ✅ Testing approach (verification strategy)

**Ready to proceed with implementation!** 🎉

---

## 💬 Questions Clarified

1. **"What is PCA correcting for?"**
   → Receiver/electronics variations (gain, baseline, temperature drift)
   → NOT atmospheric effects

2. **"What are the two line types?"**
   → Fixed instrumental artifact (from config, hardware issue)
   → Science line (auto-detected, preserve it!)

3. **"Why mask both lines?"**
   → Artifact: It's a constant, not a learnable pattern
   → Science: We want to protect it, not remove it

4. **"How do we prevent PCA from removing CII?"**
   → Mask the CII channels during decomposition
   → Don't include CII in what PCA learns
   → Only apply PCA to non-CII channels during correction

5. **"Is this about atmospheric correction?"**
   → No! It's about receiver/electronics correction
   → Atmospheric effects aren't part of receiver variations

---

Next steps: Start Phase 1 implementation (foundation) 🚀
