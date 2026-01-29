# filter_fits Quick Reference - Updated Defaults

## TL;DR

**NEW DEFAULT:** NaN filtering applies to **ALL spectra**

```bash
filter_fits --config config.toml
# Filters ALL objects by NaN threshold
```

**To filter ONLY target object (old default behavior):**

```bash
filter_fits --config config.toml --apply-only-to-object
# Filters only M51, keeps all other objects
```

---

## Behavior Matrix

| Command | M51 < 20% NaN | M51 >= 20% NaN | Other < 20% NaN | Other >= 20% NaN |
|---------|---------------|----------------|-----------------|------------------|
| `filter_fits --config config.toml` | Clean ✓ | Rejected ✗ | Clean ✓ | **Rejected ✗** |
| `--apply-only-to-object` | Clean ✓ | Rejected ✗ | **Clean ✓** | **Clean ✓** |

---

## Common Scenarios

### Scenario 1: I want high-quality data
```bash
filter_fits --config config.toml
# Default: All objects filtered, only best spectra kept
```

### Scenario 2: I want target data with reference objects
```bash
filter_fits --config config.toml --apply-only-to-object
# M51 filtered, other objects kept for reference
```

### Scenario 3: Very strict quality control
```bash
filter_fits --config config.toml --nan-threshold 0.05
# All objects must be pristine (< 5% NaNs)
```

### Scenario 4: Target-focused with reference
```bash
filter_fits --config config.toml \
    --nan-threshold 0.10 \
    --apply-only-to-object
# M51 very strict, others kept regardless
```

---

## What Changed

| Aspect | Before | After |
|--------|--------|-------|
| **Default filtering** | Target object only | ALL spectra |
| **Flag name** | `--apply-to-all` | `--apply-only-to-object` |
| **Flag behavior** | Enable all-filtering | Disable all-filtering |
| **Most common case** | Add flag | Use default |

---

## Migration Guide

**If you have existing scripts:**

**Old way (still works, but different default):**
```bash
filter_fits --config config.toml
# OLD: Filtered M51 only
# NEW: Filters ALL objects
```

**To keep old behavior, add flag:**
```bash
filter_fits --config config.toml --apply-only-to-object
# Filters M51 only (same as old default)
```

---

## Remember

✓ **Default (no flag)** = Filter all spectra equally
✗ **With flag** = Filter target object only, keep everything else

