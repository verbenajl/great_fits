# Velocity Axis Documentation: Complete Collection

**Date**: January 20, 2026  
**Status**: Comprehensive Documentation Suite Complete  
**Focus**: X-axis is ALWAYS velocity - read from VELOCITY_AXIS header

---

## 📚 Complete Velocity Axis Documentation Suite

We have created **three complementary documents** that together provide complete, redundant coverage of velocity axis handling:

### 1. **VELOCITY_AXIS_DEFINITIVE.md** (12 KB) - COMPREHENSIVE REFERENCE
**File**: `/docs/VELOCITY_AXIS_DEFINITIVE.md`

**Purpose**: Definitive, complete reference for velocity axis handling

**Contains**:
- The fundamental rule (X-axis is ALWAYS velocity)
- Complete specification (source, units, type)
- Conversion flow (step-by-step m/s → km/s)
- Channel mapping algorithm (nearest-neighbor)
- Real-world implementation code flow
- Verification methods (4 detailed checks)
- Common mistakes and solutions (5 mistakes)
- Integration checklist (10 items)
- Code patterns to follow
- Golden rules summary (8 rules)

**Best For**:
- Understanding complete velocity axis system
- Implementing new velocity-dependent features
- Debugging velocity-related issues
- Verification and validation
- Training new developers

### 2. **VELOCITY_AXIS_REFERENCE.md** (6.3 KB) - IMPLEMENTATION FOCUS
**File**: `/docs/VELOCITY_AXIS_REFERENCE.md`

**Purpose**: Implementation-focused reference with code examples

**Contains**:
- Data source specification
- Config file integration
- Conversion flow
- Channel mapping function documentation
- Usage in PCA correction (lines of code)
- Testing and verification methods
- Common issues and solutions
- Integration checklist

**Best For**:
- Quick implementation reference
- Code examples and patterns
- Troubleshooting specific issues
- Integration requirements
- Fast lookup

### 3. **VELOCITY_AXIS_ORGANIZATION.md** (9.5 KB) - EXISTING ORGANIZATIONAL REFERENCE
**File**: `/docs/VELOCITY_AXIS_ORGANIZATION.md`

**Purpose**: Organizational and contextual reference

**Contains**:
- How velocity axis fits in overall architecture
- Organizational patterns and structure
- Configuration integration
- Data flow through pipeline
- Existing documentation references

**Best For**:
- Understanding system architecture
- Seeing where velocity axis fits
- Organizational context
- Cross-referencing other docs

---

## 🎯 Which Document to Read?

### I want to understand EVERYTHING about velocity axis
→ Read **VELOCITY_AXIS_DEFINITIVE.md** (start to finish)

### I need to implement a velocity-dependent feature
→ Read **VELOCITY_AXIS_DEFINITIVE.md** (implementation section) + code examples

### I'm troubleshooting a velocity issue
→ Read **VELOCITY_AXIS_REFERENCE.md** (problems/solutions section)

### I need quick code examples
→ Read **VELOCITY_AXIS_REFERENCE.md** (data flow section)

### I need architectural context
→ Read **VELOCITY_AXIS_ORGANIZATION.md**

### I need complete redundancy/confirmation
→ Read all three (they're complementary, not redundant)

---

## 📋 The Core Message (All Three Documents Agree)

### ⚠️ The Fundamental Rule

**THE X-AXIS IS ALWAYS A VELOCITY AXIS**

### How It Works
1. **Read** from VELOCITY_AXIS column in FITS file
2. **Units** are meters per second (m/s) - NOT km/s!
3. **Convert** to km/s by dividing by 1000
4. **Map** velocity windows to channel indices
5. **Protect** science lines during correction

### Critical Values
```
FITS VELOCITY_AXIS column:       m/s (e.g., 6,000,000)
User parameters (config.toml):   km/s (e.g., [450, 500])
Conversion rule:                 velocity_kms = velocity_ms / 1000.0
Channel mapping:                 argmin(abs(velocity_axis - target))
```

---

## 🔗 Document Interconnection

```
VELOCITY_AXIS_DEFINITIVE.md (Comprehensive)
    ├─→ For complete understanding
    ├─→ For implementing new features
    ├─→ For debugging
    └─→ For training

VELOCITY_AXIS_REFERENCE.md (Implementation-focused)
    ├─→ For quick lookup
    ├─→ For code examples
    ├─→ For troubleshooting
    └─→ For integration

VELOCITY_AXIS_ORGANIZATION.md (Architectural context)
    ├─→ For understanding overall system
    ├─→ For cross-referencing
    ├─→ For organizational context
    └─→ For seeing relationships

PCA_QUICK_REFERENCE_CARD.md (Quick summary)
    └─→ Printable summary of key concepts
```

---

## ✅ Coverage Verification

### Fundamental Rule
- ✓ DEFINITIVE: Section 1 - "🚨 The Fundamental Rule"
- ✓ REFERENCE: Opening statement
- ✓ ORGANIZATION: How velocity axis fits in system

### Complete Specification
- ✓ DEFINITIVE: Section "📋 The Velocity Axis: Complete Specification"
- ✓ REFERENCE: "Data Source" section
- ✓ ORGANIZATION: "Velocity Axis Structure"

### Conversion Flow
- ✓ DEFINITIVE: Section "🔄 The Conversion Flow" (step-by-step)
- ✓ REFERENCE: "Conversion Flow" section
- ✓ ORGANIZATION: "How Velocity Data Flows"

### Channel Mapping
- ✓ DEFINITIVE: Section "🔍 Channel Mapping: The Critical Function"
- ✓ REFERENCE: "Channel Mapping Function"
- ✓ ORGANIZATION: "Velocity to Channel Conversion"

### Real-World Implementation
- ✓ DEFINITIVE: Section "📊 Real-World Flow in pca_correct_fits.py"
- ✓ REFERENCE: "Usage in PCA correction"
- ✓ ORGANIZATION: "Integration in Pipeline"

### Verification Methods
- ✓ DEFINITIVE: Section "✅ Verification: How to Check Velocity Axis"
- ✓ REFERENCE: "Testing & Verification"
- ✓ ORGANIZATION: "Validation Points"

### Common Mistakes
- ✓ DEFINITIVE: Section "⚠️ Common Mistakes to Avoid"
- ✓ REFERENCE: "Common Issues"
- ✓ ORGANIZATION: "Potential Problems"

### Golden Rules
- ✓ DEFINITIVE: Section "🎯 Summary: The Golden Rules"
- ✓ REFERENCE: "Golden Rules"
- ✓ ORGANIZATION: "Key Principles"

---

## 🚀 Quick Start: Where to Begin

### For a New Developer
1. Read: **VELOCITY_AXIS_DEFINITIVE.md** sections 1-3 (15 min)
2. Review: Code examples in section "📊 Real-World Flow"
3. Run: Verification checks from section "✅ Verification"
4. Bookmark: For future reference

### For Integration/Implementation
1. Read: **VELOCITY_AXIS_DEFINITIVE.md** "🔗 Integration Checklist"
2. Review: "📝 Code Pattern to Follow"
3. Check: "⚠️ Common Mistakes to Avoid"
4. Follow: Template implementation provided

### For Troubleshooting
1. Read: **VELOCITY_AXIS_REFERENCE.md** "Common Issues"
2. Check: **VELOCITY_AXIS_DEFINITIVE.md** "Verification Methods"
3. Review: Diagnostic output with logging
4. Test: Channel mapping with real data

---

## 📍 File Locations

### Velocity Axis Documentation
- `/docs/VELOCITY_AXIS_DEFINITIVE.md` (12 KB) - **START HERE**
- `/docs/VELOCITY_AXIS_REFERENCE.md` (6.3 KB) - Quick reference
- `/docs/VELOCITY_AXIS_ORGANIZATION.md` (9.5 KB) - Context

### Related Documentation
- `/docs/PCA_QUICK_REFERENCE_CARD.md` - Velocity axis summary table
- `/docs/pca_correction/README.md` - PCA correction hub
- `/docs/INDEX.md` - Master documentation index
- `/docs/PCA_DOCUMENTATION_MAP.md` - Documentation structure

### Implementation
- `/src/oi_zeigt/pca_analysis/pca_correct_fits.py` - Source code
  - Lines 419-426: Load velocity axis
  - Lines 450-465: Use with config window
  - Lines 123-159: Channel mapping function

---

## 🎓 Learning Path

### Path 1: Quick Learner (30 minutes)
1. Read: VELOCITY_AXIS_DEFINITIVE.md "🚨 The Fundamental Rule" (5 min)
2. Read: VELOCITY_AXIS_DEFINITIVE.md "🔄 The Conversion Flow" (5 min)
3. Review: Code examples in "📊 Real-World Flow" (10 min)
4. Check: Verification methods in "✅ Verification" (10 min)

### Path 2: Complete Understanding (90 minutes)
1. Read: VELOCITY_AXIS_DEFINITIVE.md (60 min)
2. Review: VELOCITY_AXIS_REFERENCE.md (20 min)
3. Check: Implementation in pca_correct_fits.py (10 min)

### Path 3: Implementation Focus (60 minutes)
1. Read: VELOCITY_AXIS_DEFINITIVE.md sections 1-4 (30 min)
2. Review: "🔗 Integration Checklist" (10 min)
3. Review: "📝 Code Pattern to Follow" (10 min)
4. Check: Implementation examples (10 min)

### Path 4: Troubleshooting (45 minutes)
1. Check: Specific error in VELOCITY_AXIS_DEFINITIVE.md "⚠️ Common Mistakes"
2. Review: "✅ Verification: How to Check"
3. Test: With actual FITS file
4. Review: Diagnostic output

---

## 🔐 Key Principles (All Three Documents Reinforce)

1. **X-axis is velocity** - Not wavelength, not channels
2. **Read from FITS** - Don't assume or calculate values
3. **FITS uses m/s** - Convert to km/s by dividing by 1000
4. **User params use km/s** - No conversion needed
5. **Nearest-neighbor mapping** - `argmin(abs_difference)`
6. **Science lines protected** - Via velocity window conversion
7. **Verify with real data** - Test on actual FITS files
8. **Document assumptions** - Make velocity handling clear

---

## ✨ Why Three Documents?

**One document might be missed or incomplete.**

Three complementary documents ensure:
- ✓ **Comprehensive coverage** - Multiple angles and depths
- ✓ **Redundancy** - Key points reinforced
- ✓ **Cross-checking** - Verify consistency
- ✓ **Multiple learning styles** - Different perspectives
- ✓ **Quick reference** - Find what you need fast
- ✓ **Training material** - For new developers
- ✓ **Implementation guide** - For coding
- ✓ **Troubleshooting help** - For debugging

---

## 📊 Documentation Statistics

| Document | Size | Focus | Best For |
|----------|------|-------|----------|
| VELOCITY_AXIS_DEFINITIVE.md | 12 KB | Comprehensive | Complete understanding |
| VELOCITY_AXIS_REFERENCE.md | 6.3 KB | Implementation | Quick reference |
| VELOCITY_AXIS_ORGANIZATION.md | 9.5 KB | Context | Architectural overview |
| **Total Velocity Axis Docs** | **27.8 KB** | **Complete Suite** | **Everything** |

---

## 🎯 The Bottom Line

### If You Remember Nothing Else, Remember This:

```
The X-axis is ALWAYS a velocity axis.

It is read from the VELOCITY_AXIS column in the FITS binary table.

The units are METERS per SECOND (m/s), not kilometers per second.

For user-facing parameters (config.toml, CLI), divide by 1000 to get km/s.

Use the converted velocity array to map velocity windows to channel indices
via nearest-neighbor matching: argmin(abs(velocity_axis - target_velocity))

This determines which channels contain science lines and must be protected.
```

---

## ✅ Status

**Documentation Suite**: Complete ✓  
**Coverage**: Comprehensive ✓  
**Implementation**: Verified ✓  
**Redundancy**: Confirmed ✓  
**Quality**: Production Ready ✓

---

**Version**: 1.0  
**Last Updated**: January 20, 2026  
**Maintainer**: PCA Analysis Module  
**Status**: Ready for Production Use

Start with: **VELOCITY_AXIS_DEFINITIVE.md**

