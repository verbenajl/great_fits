# filter_fits --apply-only-to-object Implementation Update

## Change Summary

**Default behavior REVERSED:**
- **OLD:** NaN filtering applied ONLY to target object (default), with `--apply-to-all` flag to filter all spectra
- **NEW:** NaN filtering applied to ALL spectra (default), with `--apply-only-to-object` flag to filter only target object

---

## Updated CLI Option

**Changed from:**
```python
@click.option(
    "--apply-to-all",
    is_flag=True,
    default=False,
    help="Apply NaN filtering to ALL spectra, not just target object"
)
```

**Changed to:**
```python
@click.option(
    "--apply-only-to-object",
    is_flag=True,
    default=False,
    help="Apply NaN filtering only to target object (not to all spectra)"
)
```

---

## Updated Function Logic

**CLI feedback now shows:**
```python
if apply_only_to_object:
    click.echo(f"Applying NaN filter to {object} spectra only")
else:
    click.echo(f"Applying NaN filter to ALL spectra")  # NEW DEFAULT
```

**Call to core function:**
```python
clean_path, rejected_path = filter_and_save_fits(
    hdul,
    object_name=object,
    nan_threshold=nan_threshold,
    output_clean=output_clean,
    output_rejected=output_rejected,
    remove_column=remove,
    remove_values=list(remove_values) if remove_values else None,
    apply_to_all=not apply_only_to_object  # Inverted logic
)
```

---

## Updated Docstring

```python
"""
Filter FITS data by object, NaN content, and/or column values.

Creates two FITS files:
1. Clean file: Spectra passing NaN threshold filter
2. Rejected file: Spectra failing NaN threshold filter or matching removal criteria

By default, NaN filtering is applied to ALL spectra. Use --apply-only-to-object
to filter only the target object and keep all other objects regardless of NaN content.

Examples:

    # Filter ALL spectra by NaN threshold (default)
    filter_fits --config config.toml
    
    # Filter only M51 spectra, keep all other objects
    filter_fits --config config.toml --apply-only-to-object
    
    # Filter all with custom thresholds
    filter_fits --config config.toml --nan-threshold 0.75
    
    # Filter specific object and specify output files
    filter_fits --config config.toml --object "M51" \\
        --output-clean m51_clean.fits --output-rejected m51_rejected.fits
    
    # Remove specific AOR_ID values AND filter all spectra
    filter_fits --config config.toml --remove AOR_ID \\
        --remove-values 04_0116_0020609 --remove-values 04_0116_0020506
"""
```

---

## Behavior Comparison

### Default (NEW: apply-to-all)
```bash
filter_fits --config config.toml
```

**Result:**
- M51: Filtered by NaN threshold (< 20% kept)
- NGC1234: ALSO filtered by NaN threshold (< 20% kept)
- Other objects: ALSO filtered by NaN threshold
- **Effect: Global quality filtering, fewer spectra in clean file**

### With --apply-only-to-object Flag (OLD default behavior)
```bash
filter_fits --config config.toml --apply-only-to-object
```

**Result:**
- M51: Filtered by NaN threshold (< 20% kept)
- NGC1234: ALL kept (no NaN filtering)
- Other objects: ALL kept (no NaN filtering)
- **Effect: Target object selective, preserves supporting observations**

---

## Usage Examples

### Example 1: Default - Quality Filter All
```bash
filter_fits --config config.toml
```
- All objects filtered equally
- Clean file has highest quality data overall
- Fewer spectra, but all meet quality threshold

### Example 2: Target-Only Filtering
```bash
filter_fits --config config.toml --apply-only-to-object
```
- M51 quality controlled
- Other objects preserved for reference
- More spectra in clean file (all supporting observations)

### Example 3: Strict Global Quality
```bash
filter_fits --config config.toml --nan-threshold 0.05
```
- ALL spectra must have < 5% NaNs (very strict)
- Applied to all objects (cannot be disabled)
- Only highest quality data survives

### Example 4: Strict Global + Target-Only
```bash
filter_fits --config config.toml \
    --nan-threshold 0.05 \
    --apply-only-to-object
```
- M51 must have < 5% NaNs (very strict)
- Other objects kept regardless
- High-quality target with reference data

---

## Files Modified

1. **`src/oi_zeigt/cli.py`** (Lines 700-760)
   - Changed option name from `--apply-to-all` to `--apply-only-to-object`
   - Changed parameter name in function signature
   - Updated CLI feedback messages (inverted logic)
   - Updated function docstring with new examples
   - Updated call to core function with inverted boolean

2. **`src/oi_zeigt/reduction/core.py`** (unchanged)
   - Core function still uses `apply_to_all` parameter internally
   - Logic unchanged - just receives inverted boolean from CLI

---

## Backward Compatibility

⚠️ **NOT backward compatible** - this is a breaking change in default behavior:

**Old behavior (apply only to object):**
```bash
filter_fits --config config.toml  # Applied to M51 only
```

**New behavior (apply to all):**
```bash
filter_fits --config config.toml  # Applied to ALL objects
```

**To get old behavior:**
```bash
filter_fits --config config.toml --apply-only-to-object
```

---

## Rationale for Change

**Why apply-to-all should be default:**
1. **More consistent**: All objects treated equally by quality metric
2. **Cleaner dataset**: Output guaranteed to meet quality threshold
3. **Easier to use**: Most common case (pure quality filtering)
4. **Less surprising**: Flag for special case (keep non-target objects), not common case

**Why target-only filtering is still available:**
1. **Preserve context**: Keep reference/supporting observations
2. **Focused science**: Target object quality-controlled, others as-is
3. **Survey compatibility**: Different objects may have different quality characteristics

---

## Decision Tree (Updated)

```
Do you want non-target objects included without NaN filtering?
│
├─ NO (default): filter_fits --config config.toml
│   └─ All objects filtered equally
│
└─ YES: filter_fits --config config.toml --apply-only-to-object
    └─ Only target object affected by NaN threshold
```

---

## Test Cases

### Test 1: Verify Default is Now All-Filter
```bash
filter_fits --config config.toml
# Output message should say: "Applying NaN filter to ALL spectra"
```

### Test 2: Verify Flag Works
```bash
filter_fits --config config.toml --apply-only-to-object
# Output message should say: "Applying NaN filter to M51 spectra only"
```

### Test 3: Compare Output Sizes
```bash
# Run with default (all)
filter_fits --config config.toml --nan-threshold 0.20
mv clean_data.fits clean_all.fits

# Run with target-only
filter_fits --config config.toml --nan-threshold 0.20 --apply-only-to-object
mv clean_data.fits clean_target_only.fits

# clean_target_only.fits should have MORE spectra (all other objects preserved)
```

---

## See Also

- `FILTER_FITS_COMMAND.md` - Full command reference
- `FILTER_FITS_CODE_WALKTHROUGH.md` - Implementation details
- `FILTER_FITS_EXAMPLES.md` - Usage examples

