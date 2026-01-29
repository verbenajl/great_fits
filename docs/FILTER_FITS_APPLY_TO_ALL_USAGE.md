# filter_fits --apply-to-all Usage Guide

## Quick Start

### Default Behavior (Target Object Only)
```bash
filter_fits --config config.toml
```
- Filters **M51** spectra by NaN threshold
- Keeps **ALL other objects** regardless of NaN content

### New: Filter All Spectra
```bash
filter_fits --config config.toml --apply-to-all
```
- Filters **ALL spectra** by NaN threshold
- Non-target objects also removed if they exceed NaN threshold

---

## Common Use Cases

### 1. Quality Check - Keep Everything Except Known Bad
```bash
# Standard filtering (M51 only affected)
filter_fits --config config.toml
```

**Result:** 
- M51: ~80% (< 20% NaN)
- NGC1234: 100% (no filtering)
- Other objects: 100% (no filtering)

---

### 2. Strict Quality - Filter Everything
```bash
# Apply NaN threshold to ALL spectra
filter_fits --config config.toml --apply-to-all
```

**Result:**
- M51: ~80% (< 20% NaN)
- NGC1234: ~75% (< 20% NaN)
- Other objects: ~70% (< 20% NaN)
- **Overall cleaner dataset, fewer spectra**

---

### 3. Very Strict Quality
```bash
# Combine: apply-to-all + stricter threshold
filter_fits --config config.toml --apply-to-all \
    --nan-threshold 0.05
```

**Result:**
- ALL objects: Only spectra with < 5% NaNs kept
- ~50-60% of data retained
- **Very high quality, but reduced statistics**

---

### 4. Mix Strategies - Remove Bad AND Filter
```bash
# Remove problematic observations AND apply global NaN filter
filter_fits --config config.toml --apply-to-all \
    --remove AOR_ID \
    --remove-values bad_obs_1 --remove-values bad_obs_2
```

**Processing:**
1. Remove rows with specific AOR_IDs → rejected file
2. Apply NaN threshold to ALL remaining spectra
3. Output: Clean (passed NaN filter) + Rejected (AOR_ID + high NaN)

---

## Decision Tree

```
Do you want to filter ALL spectra equally?
│
├─ NO (default): filter_fits --config config.toml
│   └─ Only target object affected by NaN threshold
│
└─ YES: filter_fits --config config.toml --apply-to-all
    └─ All objects affected by NaN threshold equally
```

---

## Output File Differences

### With `--apply-to-all=False` (default)
```
Input: 16,000 total spectra
  - M51: 10,000
  - Other: 6,000

Output:
  Clean:     15,000
    - M51 (filtered):  8,000 (80% kept)
    - Other (all):     7,000 (100% kept)
  
  Rejected:  1,000
    - M51 only:        2,000 (20% rejected)
```

### With `--apply-to-all=True`
```
Input: 16,000 total spectra
  - M51: 10,000
  - Other: 6,000

Output:
  Clean:     11,600
    - M51 (filtered):  8,000 (80% kept)
    - Other (filtered): 3,600 (60% kept)
  
  Rejected:  4,400
    - M51:  2,000 (20% rejected)
    - Other: 2,400 (40% rejected)
```

---

## Tips

**✓ Use `--apply-to-all` when:**
- Making a uniform-quality dataset
- All objects should meet same standards
- Publishing results with consistent quality
- Creating training/testing sets
- Combining observations from multiple sources

**✗ Don't use `--apply-to-all` when:**
- Target object is important, others are secondary
- Other objects have fundamentally different data quality
- You want to preserve all supporting observations
- You only care about one object's quality

---

## Examples with Real Data

### Example 1: M51 Study (Target-Only Filtering)
```bash
# M51 is target, use all available supporting observations
filter_fits --config config.toml
# M51 filtered, but NGC1234 kept even if lower quality
```

### Example 2: Galaxy Survey (Global Filtering)
```bash
# Multiple galaxies, all treated equally
filter_fits --config config.toml --apply-to-all
# All galaxies filtered by same NaN threshold
```

### Example 3: Publication-Ready Data (Very Strict)
```bash
# For high-impact paper: very strict quality
filter_fits --config config.toml --apply-to-all \
    --nan-threshold 0.05 \
    --remove AOR_ID \
    --remove-values known_bad_1 \
    --remove-values known_bad_2

# Result: High-quality subset suitable for publication
```

---

## FAQ

**Q: What's the default behavior?**
A: `apply_to_all=False` - only target object filtered by NaN

**Q: Why would I want to apply to all?**
A: For uniform quality across all objects in dataset

**Q: Can I combine with `--remove`?**
A: Yes! Both filters apply, removal happens first

**Q: Does this change the output format?**
A: No, just different spectra in clean/rejected files

**Q: Can I undo this?**
A: Yes, run again with opposite flag, but regenerated files will overwrite

---

## See Also

- `FILTER_FITS_APPLY_TO_ALL.md` - Technical implementation details
- `FILTER_FITS_COMMAND.md` - Full command reference
- `FILTER_FITS_CODE_WALKTHROUGH.md` - How it works internally

