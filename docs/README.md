# OI-ZEIGT Documentation

Welcome to the OI-ZEIGT documentation hub. All documentation is organized by topic for easy navigation.

---

## 📚 Quick Navigation

### 🚀 Getting Started
- **[Main Project README](../README.md)** - Overview and installation
- **[Complete Workflow](guides/COMPLETE_WORKFLOW.md)** - End-to-end pipeline explanation
- **[FITS Combine Guide](guides/FITS_COMBINE_GUIDE.md)** - How to combine FITS files

### 📖 User Guides
Located in [`guides/`](guides/)
- **FITS Combine Guide** - FITS file combining procedure
- **Complete Workflow** - Full pipeline walkthrough
- **Before and After** - Visual comparison of improvements
- **Velocity Axis Quick Summary** - Quick reference for velocity axis
- **README Velocity Axis Fix** - Explanation of velocity axis fixes

### 🛠️ Implementation Details
Located in [`implementation/`](implementation/)

#### Spectral Axis & Velocity Reconstruction
- **[Spectral Axis Table Implementation](implementation/SPECTRAL_AXIS_TABLE_IMPLEMENTATION.md)** ⭐ START HERE
  - How velocity table is stored in datacubes
  - How users access velocity information
  - Direct access methods (no WCS needed)

- **[Velocity Axis Reconstruction](implementation/VELOCITY_AXIS_RECONSTRUCTION.md)**
  - Velocity axis calculation formula
  - Parameter extraction from FITS
  - WCS compliance details

- **[Velocity Axis Changes](implementation/VELOCITY_AXIS_CHANGES.md)**
  - Summary of modifications to datacube creation
  - Before/after comparison
  - Key improvements

#### CRPIX1 Preservation
- **[CRPIX1 Fix](implementation/CRPIX1_FIX.md)**
  - Why CRPIX1 matters (reference pixel for velocity)
  - How CRPIX1 is preserved in single-HDU files
  - Impact on velocity axis accuracy

- **[CRPIX1 Preservation Verified](implementation/CRPIX1_PRESERVATION_VERIFIED.md)**
  - Test results confirming CRPIX1 preservation
  - Real M51 data verification

#### Comprehensive Summaries
- **[Spectral Axis Implementation Complete](implementation/SPECTRAL_AXIS_IMPLEMENTATION_COMPLETE.md)**
  - Full implementation summary
  - Expected output specifications

- **[Spectral Axis Reconstruction Success](implementation/SPECTRAL_AXIS_RECONSTRUCTION_SUCCESS.md)**
  - Verification of correct implementation
  - Test results with real data

- **[Implementation Checklist](implementation/IMPLEMENTATION_CHECKLIST.md)**
  - Step-by-step implementation tracking
  - Verification items

### 📋 Reference
Located in [`reference/`](reference/)
- **Status Spectral Axis** - Current implementation status
- **Documentation Index** - Detailed documentation inventory
- **FIX Reduce Spectra Size** - Information about spectra reduction
- **Reduce Spectra Clean Flag** - Clean flag documentation
- **Single HDU Option** - Details about single-HDU FITS format

### 📦 Archive
Located in [`archive/`](archive/)
- Old implementation notes
- Superseded documentation
- Historical records

---

## 🎯 Common Tasks

### "I want to use the datacube with proper velocity information"
→ Read: **[Spectral Axis Table Implementation](implementation/SPECTRAL_AXIS_TABLE_IMPLEMENTATION.md)**

### "I want to understand the velocity axis"
→ Read: **[Velocity Axis Reconstruction](implementation/VELOCITY_AXIS_RECONSTRUCTION.md)**

### "I want to know what changed"
→ Read: **[Velocity Axis Changes](implementation/VELOCITY_AXIS_CHANGES.md)**

### "I want the full picture"
→ Read: **[Complete Workflow](guides/COMPLETE_WORKFLOW.md)**

### "I want technical details about CRPIX1"
→ Read: **[CRPIX1 Fix](implementation/CRPIX1_FIX.md)**

---

## 📊 Documentation Structure

```
docs/
├── README.md                          ← You are here
├── implementation/                    ← Technical implementation details
│   ├── SPECTRAL_AXIS_TABLE_IMPLEMENTATION.md       ⭐ START HERE
│   ├── VELOCITY_AXIS_RECONSTRUCTION.md
│   ├── VELOCITY_AXIS_CHANGES.md
│   ├── CRPIX1_FIX.md
│   ├── CRPIX1_PRESERVATION_VERIFIED.md
│   ├── SPECTRAL_AXIS_IMPLEMENTATION_COMPLETE.md
│   ├── SPECTRAL_AXIS_RECONSTRUCTION_SUCCESS.md
│   ├── IMPLEMENTATION_CHECKLIST.md
│   ├── IMPLEMENTATION_SUMMARY.md
│   ├── IMPLEMENTATION_COMPLETE.md
│   ├── IMPLEMENTATION_PLAN_SPECTRAL_AXIS_TABLE.md
│   ├── SPECTRAL_AXIS_DATACUBE_OPTIONS.md
│   ├── SPECTRAL_AXIS_RECONSTRUCTION_VERIFICATION.md
│   ├── SPECTRAL_AXIS_BEFORE_AFTER.md
│   ├── VELOCITY_AXIS_RECONSTRUCTION_STATUS.md
│   └── VELOCITY_AXIS_WITH_VS_WITHOUT_CRPIX1.md
├── guides/                            ← User guides and tutorials
│   ├── FITS_COMBINE_GUIDE.md
│   ├── COMPLETE_WORKFLOW.md
│   ├── BEFORE_AND_AFTER.md
│   ├── README_VELOCITY_AXIS_FIX.md
│   └── VELOCITY_AXIS_QUICK_SUMMARY.md
├── reference/                         ← Technical reference
│   ├── STATUS_SPECTRAL_AXIS.md
│   ├── DOCUMENTATION_INDEX.md
│   ├── FIX_REDUCE_SPECTRA_SIZE.md
│   ├── REDUCE_SPECTRA_CLEAN_FLAG.md
│   └── SINGLE_HDU_OPTION.md
└── archive/                           ← Old/superseded documentation
    ├── 00_READ_ME_FIRST.md
    └── FINAL_SUMMARY.md
```

---

## 🔑 Key Concepts

### Spectral Axis Table
The datacube now includes a binary table (SPECTRUM HDU) with:
- **CHANNEL column**: Channel indices (0-1263)
- **VELOCITY column**: Velocity in m/s for each channel
- **Metadata**: CRPIX1, RESTFRQ, VELDEF

This allows users to access velocity information directly without WCS libraries.

### CRPIX1 (Reference Pixel)
The pixel coordinate where the reference velocity is located. For M51:
- CRPIX1 ≈ 505.85 (FITS 1-indexed)
- Reference velocity: 470 km/s
- Reference channel: ~505 (0-indexed)

### Velocity Axis Calibration
The velocity axis is now properly calibrated from FITS parameters:
- **VELOCITY** (from table column): 470,000 m/s reference
- **DELTAV** (from table column): 500 m/s per channel
- **RESTFREQ** (from table column): Rest frequency in Hz
- **VELDEF** (from table column): Velocity definition

---

## 🧪 Quick Test

To verify your datacube has the spectral axis table:

```python
from astropy.io import fits

hdul = fits.open('datacube.fits')
print(f"HDUs: {[hdu.name for hdu in hdul]}")  # Should show PRIMARY, SPECTRUM

# Access velocity information
velocities = hdul['SPECTRUM'].data['VELOCITY']
print(f"Velocity at channel 505: {velocities[505]:.0f} m/s")
# Should print: Velocity at channel 505: 470075 m/s
```

---

## 📖 Reading Guide

**If you're new to OI-ZEIGT:**
1. Start with [Main README](../README.md)
2. Read [Complete Workflow](guides/COMPLETE_WORKFLOW.md)
3. Check [Spectral Axis Table](implementation/SPECTRAL_AXIS_TABLE_IMPLEMENTATION.md)

**If you're implementing new features:**
1. Review [Implementation Checklist](implementation/IMPLEMENTATION_CHECKLIST.md)
2. Check [Status](reference/STATUS_SPECTRAL_AXIS.md)
3. Refer to relevant implementation docs

**If you're troubleshooting:**
1. Check [Spectral Axis Table](implementation/SPECTRAL_AXIS_TABLE_IMPLEMENTATION.md) (troubleshooting section)
2. Review [Status](reference/STATUS_SPECTRAL_AXIS.md)
3. Look in [Archive](archive/) for historical context

---

## ✅ Implementation Status

**Latest Changes:**
- ✅ Spectral axis table added to datacubes
- ✅ CRPIX1 preserved in single-HDU files
- ✅ Velocity axis properly calibrated
- ✅ All documentation organized
- ✅ Tests passing with real M51 data

**Version:** January 2026
**Status:** Production Ready

---

## 💡 Tips

- Use **Ctrl+F** to search across documentation
- Cross-references point to other relevant sections
- Code examples are marked with syntax highlighting
- Technical terms are linked to explanations

---

## 🔗 External Resources

- [FITS Standard](https://fits.gsfc.nasa.gov/)
- [WCS (World Coordinate System)](https://www.atnf.csiro.au/people/mcalabre/WCS/)
- [astropy.io.fits Documentation](https://docs.astropy.org/en/stable/io/fits/)
- [spectral-cube Library](https://spectral-cube.readthedocs.io/)

---

## 📝 Contributing to Documentation

When adding new documentation:
1. Choose appropriate folder (implementation, guides, or reference)
2. Use clear, descriptive filenames
3. Include a table of contents for longer documents
4. Add cross-references to related documents
5. Update this README with the new document

---

**Last Updated:** January 2026  
**Maintainer:** OI-ZEIGT Development Team
