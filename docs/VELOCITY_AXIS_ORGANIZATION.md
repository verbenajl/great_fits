# Documentation Organization & Navigation Guide

## Overview

Documentation is organized into logical folders for easy navigation. All files related to the recent velocity axis review and implementation are properly categorized.

---

## Folder Structure

```
docs/
├── README.md                           ← Main docs hub
├── REPOSITORY_GUIDE.md
├── PROJECT_STATUS.md
├── PROJECT_STATUS_SUMMARY.md
├── IMPLEMENTATION_INDEX.md
│
├── velocity_axis/                      ← ✨ NEW: Velocity axis review (10 files)
│   ├── VELOCITY_AXIS_QUICK_CHECK.md   ← ONE-PAGE SUMMARY
│   ├── VELOCITY_AXIS_SUMMARY.md
│   ├── VELOCITY_AXIS_CODE_REVIEW.md
│   ├── VELOCITY_AXIS_VISUAL_COMPARISON.md
│   ├── VELOCITY_AXIS_CALL_HIERARCHY.md
│   ├── VELOCITY_AXIS_DOCUMENTATION_INDEX.md
│   ├── VELOCITY_AXIS_QUICK_REFERENCE.md
│   ├── VELOCITY_AXIS_IMPLEMENTATION.md
│   ├── VELOCITY_AXIS_IMPLEMENTATION_SUMMARY.md
│   └── VELOCITY_AXIS_REVIEW_COMPLETE.txt
│
├── features/                           ← Feature documentation
│   └── PLOT_VELOCITY_AXIS_FEATURE.md  ← Plotting with velocity axis
│
├── implementation/                     ← Implementation details (30+ files)
│   ├── EXTRACTION_VELOCITY_IMPLEMENTATION.md
│   ├── EXTRACTION_SUMMARY.md
│   ├── EXTRACTION_CLI_IMPLEMENTATION.md
│   ├── VELOCITY_AXIS_CHANGES.md
│   ├── VELOCITY_AXIS_RECONSTRUCTION.md
│   ├── CRPIX1_FIX.md
│   ├── CRPIX1_PRESERVATION_VERIFIED.md
│   ├── SPECTRAL_AXIS_TABLE_IMPLEMENTATION.md
│   └── ... (Phase 1, 2, 3 documentation)
│
├── guides/                             ← User guides
│   ├── COMPLETE_WORKFLOW.md
│   ├── FITS_COMBINE_GUIDE.md
│   ├── VELOCITY_AXIS_QUICK_SUMMARY.md
│   └── BEFORE_AND_AFTER.md
│
├── reference/                          ← Reference documentation
│   ├── QUICK_REFERENCE.md
│   ├── DOCUMENTATION_INDEX.md
│   ├── REDUCE_SPECTRA_CLEAN_FLAG.md
│   ├── SINGLE_HDU_OPTION.md
│   ├── PCA_ANALYSIS_INDEX.md
│   └── ... (PCA, FIX, STATUS documentation)
│
└── archive/                            ← Old documentation
    └── (Historical docs)
```

---

## ✨ NEW: Velocity Axis Review Documentation

**Location:** `docs/velocity_axis/`

### Quick Navigation

**One-page overview:**
- `VELOCITY_AXIS_QUICK_CHECK.md` ⭐ **START HERE** - 2.5 KB

**Executive summaries:**
- `VELOCITY_AXIS_SUMMARY.md` - Complete findings (8.7 KB)
- `VELOCITY_AXIS_REVIEW_COMPLETE.txt` - Final review status (7.4 KB)

**Detailed analysis:**
- `VELOCITY_AXIS_CODE_REVIEW.md` - Code comparison (11 KB)
- `VELOCITY_AXIS_VISUAL_COMPARISON.md` - Side-by-side comparison (8.4 KB)
- `VELOCITY_AXIS_CALL_HIERARCHY.md` - Integration and flows (12 KB)

**Reference material:**
- `VELOCITY_AXIS_QUICK_REFERENCE.md` - Formulas and quick lookup (6.9 KB)
- `VELOCITY_AXIS_IMPLEMENTATION.md` - Tutorial style (4.2 KB)
- `VELOCITY_AXIS_DOCUMENTATION_INDEX.md` - Navigation guide (8.9 KB)
- `VELOCITY_AXIS_IMPLEMENTATION_SUMMARY.md` - Overview (5.1 KB)

### Status
✓ **APPROVED FOR PRODUCTION USE**

---

## How to Navigate

### By Task

#### "I want a quick overview"
→ `docs/velocity_axis/VELOCITY_AXIS_QUICK_CHECK.md`

#### "I need to understand the code"
→ `docs/velocity_axis/VELOCITY_AXIS_CODE_REVIEW.md`

#### "I'm comparing reduction vs gridding"
→ `docs/velocity_axis/VELOCITY_AXIS_VISUAL_COMPARISON.md`

#### "I need to understand integration"
→ `docs/velocity_axis/VELOCITY_AXIS_CALL_HIERARCHY.md`

#### "I want to use the plotting feature"
→ `docs/features/PLOT_VELOCITY_AXIS_FEATURE.md`

#### "I'm implementing spectral extraction"
→ `docs/implementation/EXTRACTION_VELOCITY_IMPLEMENTATION.md`

---

### By User Type

#### End Users
1. `docs/README.md` - Documentation hub
2. `docs/guides/COMPLETE_WORKFLOW.md` - How to use the package
3. `docs/guides/VELOCITY_AXIS_QUICK_SUMMARY.md` - Quick reference
4. `docs/velocity_axis/VELOCITY_AXIS_QUICK_CHECK.md` - Technical overview

#### Developers
1. `docs/velocity_axis/VELOCITY_AXIS_SUMMARY.md` - Status and findings
2. `docs/velocity_axis/VELOCITY_AXIS_CODE_REVIEW.md` - Code details
3. `docs/implementation/EXTRACTION_VELOCITY_IMPLEMENTATION.md` - Implementation guide
4. `docs/velocity_axis/VELOCITY_AXIS_CALL_HIERARCHY.md` - Integration details

#### Maintainers
1. `docs/PROJECT_STATUS.md` - Overall status
2. `docs/velocity_axis/VELOCITY_AXIS_REVIEW_COMPLETE.txt` - Review status
3. `docs/velocity_axis/VELOCITY_AXIS_CODE_REVIEW.md` - Detailed code review
4. `docs/REPOSITORY_GUIDE.md` - Repository structure

#### Learning/Training
1. `docs/velocity_axis/VELOCITY_AXIS_IMPLEMENTATION.md` - Tutorial
2. `docs/velocity_axis/VELOCITY_AXIS_QUICK_REFERENCE.md` - Quick reference
3. `docs/guides/COMPLETE_WORKFLOW.md` - Full workflow
4. `docs/velocity_axis/VELOCITY_AXIS_DOCUMENTATION_INDEX.md` - Detailed navigation

---

## File Organization Summary

### Velocity Axis (10 files, 98 KB total)
- **Review status documents:** 3 files
- **Detailed analysis:** 3 files
- **Reference materials:** 4 files
- **Status:** ✓ COMPLETE & APPROVED

### Features (1 file, 3.3 KB)
- **Plotting integration:** 1 file
- **Status:** ✓ COMPLETE

### Implementation (30+ files)
- **Extraction guides:** 3 files
- **Baseline documentation:** Various
- **CLI/config documentation:** Various
- **Phase-based implementation:** 3 phases documented
- **Status:** ✓ COMPLETE

### Guides (4 files)
- **Workflow documentation:** 1 file
- **FITS operations:** 1 file
- **Before/after analysis:** 1 file
- **Quick summaries:** 2 files
- **Status:** ✓ COMPLETE

### Reference (15+ files)
- **Quick references:** Multiple
- **Technical specifications:** PCA, CRPIX1, etc.
- **Test results:** Verification files
- **Status:** ✓ COMPLETE

---

## Finding Specific Information

### Velocity Axis Formula
- `docs/velocity_axis/VELOCITY_AXIS_QUICK_REFERENCE.md` (Section: Key Formulas)
- `docs/velocity_axis/VELOCITY_AXIS_CODE_REVIEW.md` (Section: Formula Verification)

### Spectral Extraction
- `docs/implementation/EXTRACTION_VELOCITY_IMPLEMENTATION.md`
- `docs/guides/COMPLETE_WORKFLOW.md`

### Baseline Subtraction
- `docs/implementation/EXTRACTION_VELOCITY_IMPLEMENTATION.md`
- `docs/guides/COMPLETE_WORKFLOW.md`

### Plotting with Velocity Axis
- `docs/features/PLOT_VELOCITY_AXIS_FEATURE.md`
- `docs/guides/COMPLETE_WORKFLOW.md`

### Gridding and Datacubes
- `docs/implementation/SPECTRAL_AXIS_TABLE_IMPLEMENTATION.md`
- `docs/implementation/VELOCITY_AXIS_RECONSTRUCTION.md`

### FITS and File Handling
- `docs/guides/FITS_COMBINE_GUIDE.md`
- `docs/reference/SINGLE_HDU_OPTION.md`

---

## Key Documentation Files

### Status & Overview
- `docs/PROJECT_STATUS.md` - Current project status
- `docs/PROJECT_STATUS_SUMMARY.md` - Quick status
- `docs/velocity_axis/VELOCITY_AXIS_REVIEW_COMPLETE.txt` - Review approval
- `docs/REPOSITORY_GUIDE.md` - Repository structure

### Velocity Axis (Complete Review)
- `docs/velocity_axis/VELOCITY_AXIS_QUICK_CHECK.md` - 1-page summary
- `docs/velocity_axis/VELOCITY_AXIS_SUMMARY.md` - Executive summary
- `docs/velocity_axis/VELOCITY_AXIS_CODE_REVIEW.md` - Detailed review

### Implementation Details
- `docs/implementation/EXTRACTION_VELOCITY_IMPLEMENTATION.md` - Extraction guide
- `docs/implementation/VELOCITY_AXIS_RECONSTRUCTION.md` - WCS info
- `docs/implementation/CRPIX1_PRESERVATION_VERIFIED.md` - Reference pixel handling

### User Guides
- `docs/guides/COMPLETE_WORKFLOW.md` - Full pipeline
- `docs/guides/VELOCITY_AXIS_QUICK_SUMMARY.md` - Quick reference
- `docs/guides/FITS_COMBINE_GUIDE.md` - FITS operations

---

## Statistics

### Total Documentation
- **Files:** 70+ files
- **Total size:** ~500+ KB
- **Organization:** 7 main folders
- **Status:** ✓ Fully organized and indexed

### Velocity Axis Specifically
- **New files:** 10 files
- **Total size:** 98 KB
- **Status:** ✓ APPROVED FOR PRODUCTION USE
- **Last updated:** January 20, 2026

---

## Version Control

**Most recent additions:** January 20, 2026

**Latest files in `/velocity_axis/`:**
- All 10 files created January 20, 2026
- Reviewed and approved for production use

---

## Access Tips

### Via Command Line
```bash
# Browse velocity axis documentation
ls -lh docs/velocity_axis/

# Read a specific file
cat docs/velocity_axis/VELOCITY_AXIS_QUICK_CHECK.md

# Find files about a topic
grep -r "extraction" docs/ | grep -E "\.md:" | head -10
```

### Via File Explorer
1. Navigate to `docs/` folder
2. Choose topic folder (velocity_axis, features, implementation, etc.)
3. Pick specific file

### Using the Documentation Index
- `docs/velocity_axis/VELOCITY_AXIS_DOCUMENTATION_INDEX.md` - Velocity axis topics
- `docs/IMPLEMENTATION_INDEX.md` - Implementation topics
- `docs/reference/DOCUMENTATION_INDEX.md` - Reference topics

---

## Next Steps

1. **Review velocity axis status:** `docs/velocity_axis/VELOCITY_AXIS_QUICK_CHECK.md`
2. **Understand the implementation:** `docs/velocity_axis/VELOCITY_AXIS_CODE_REVIEW.md`
3. **Use the features:** `docs/guides/COMPLETE_WORKFLOW.md`
4. **Get quick references:** `docs/velocity_axis/VELOCITY_AXIS_QUICK_REFERENCE.md`

---

## Need Help?

- **Quick answer:** Check `docs/velocity_axis/VELOCITY_AXIS_DOCUMENTATION_INDEX.md`
- **Understanding code:** Read `docs/velocity_axis/VELOCITY_AXIS_CODE_REVIEW.md`
- **Using features:** See `docs/guides/` folder
- **Technical details:** Check `docs/implementation/` folder

---

**Last Updated:** January 20, 2026  
**Status:** ✓ PRODUCTION READY
