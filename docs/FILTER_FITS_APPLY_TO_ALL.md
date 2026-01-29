# filter_fits --apply-to-all Implementation Summary

## Changes Made

### 1. CLI Command Updated (`src/oi_zeigt/cli.py`)

**Added option:**
```python
@click.option(
    "--apply-to-all",
    is_flag=True,
    default=False,
    help="Apply NaN filtering to ALL spectra, not just target object"
)
```

**Updated function signature:**
```python
def filter_fits(config, fits, object, nan_threshold, 
                output_clean, output_rejected, remove, 
                remove_values, apply_to_all):
```

**Updated CLI feedback:**
```python
if apply_to_all:
    click.echo(f"Applying NaN filter to ALL spectra")
else:
    click.echo(f"Applying NaN filter to {object} spectra only")
```

**Updated function call:**
```python
clean_path, rejected_path = filter_and_save_fits(
    hdul,
    object_name=object,
    nan_threshold=nan_threshold,
    output_clean=output_clean,
    output_rejected=output_rejected,
    remove_column=remove,
    remove_values=list(remove_values) if remove_values else None,
    apply_to_all=apply_to_all  # NEW PARAMETER
)
```

### 2. Core Function Updated (`src/oi_zeigt/reduction/core.py`)

**Updated function signature:**
```python
def filter_and_save_fits(hdul, object_name, nan_threshold=0.20,
                        output_clean=None, output_rejected=None,
                        remove_column=None, remove_values=None,
                        apply_to_all=False):  # NEW PARAMETER
```

**Updated docstring:**
- Added `apply_to_all` parameter documentation
- Added examples for both modes
- Clarified behavior differences

**Updated filtering logic:**

**If apply_to_all=False (DEFAULT):**
```
1. Separate data by object
2. Keep ALL non-target objects (100%, no NaN filtering)
3. Filter target object by NaN threshold
4. Combine: other_data + target_clean → clean file
5. Rejected: target_rejected + removed_data → rejected file
```

**If apply_to_all=True (NEW):**
```
1. Apply NaN filtering to ALL spectra
2. Keep spectra with NaN < threshold
3. Reject spectra with NaN >= threshold
4. Output: all_clean → clean file
5. Output: all_rejected + removed_data → rejected file
```

---

## Usage Examples

### Default Behavior (Target Object Only)
```bash
filter_fits --config config.toml
# OR explicitly:
filter_fits --config config.toml --apply-to-all False
```

**Result:**
- M51 spectra: Filtered by 20% NaN threshold
- Other objects: ALL kept (no NaN filtering)

### New Behavior (All Spectra)
```bash
filter_fits --config config.toml --apply-to-all
```

**Result:**
- M51 spectra: Filtered by 20% NaN threshold
- Other objects: ALSO filtered by 20% NaN threshold
- All spectra equally treated

### Combined with Other Options
```bash
# Filter all spectra with stricter threshold and remove bad AOR_IDs
filter_fits --config config.toml --apply-to-all \
    --nan-threshold 0.10 \
    --remove AOR_ID \
    --remove-values bad_id_1 --remove-values bad_id_2
```

---

## Behavior Comparison

| Scenario | apply_to_all=False | apply_to_all=True |
|----------|-------------------|------------------|
| **M51 with < 20% NaN** | ✓ Kept in clean | ✓ Kept in clean |
| **M51 with >= 20% NaN** | ✗ Moved to rejected | ✗ Moved to rejected |
| **NGC1234 with < 20% NaN** | ✓ Kept in clean (no filter) | ✓ Kept in clean |
| **NGC1234 with >= 20% NaN** | ✓ Kept in clean (no filter) | ✗ Moved to rejected |
| **Other object with 90% NaN** | ✓ Kept in clean (no filter) | ✗ Moved to rejected |

---

## Implementation Details

### Modified Section in filter_and_save_fits()

```python
# Separate data by object (or use all data if apply_to_all)
if apply_to_all:
    # Apply NaN filtering to ALL spectra
    all_data = data
    
    # Filter all data by NaN content
    nan_fractions = []
    for spectrum in all_data['SPECTRUM']:
        _, frac = detect_nan_channels(spectrum)
        nan_fractions.append(frac)
    
    nan_fractions = np.array(nan_fractions)
    clean_mask = nan_fractions < nan_threshold
    rejected_mask = ~clean_mask
    
    clean_combined = all_data[clean_mask]
    all_rejected = all_data[rejected_mask]
    
    # Add the removed data to rejected
    if len(removed_data) > 0:
        all_rejected = np.concatenate([all_rejected, removed_data])
else:
    # Apply NaN filtering only to target object (default behavior)
    # [Original logic preserved]
    target_mask = np.array([object_name.lower() in str(obj).lower() 
                           for obj in data['OBJECT']])
    other_mask = ~target_mask
    
    other_data = data[other_mask]
    target_data = data[target_mask]
    
    # Filter target object by NaN content
    nan_fractions = []
    for spectrum in target_data['SPECTRUM']:
        _, frac = detect_nan_channels(spectrum)
        nan_fractions.append(frac)
    
    nan_fractions = np.array(nan_fractions)
    clean_target_mask = nan_fractions < nan_threshold
    rejected_target_mask = ~clean_target_mask
    
    target_clean = target_data[clean_target_mask]
    target_rejected = target_data[rejected_target_mask]
    
    # Add the removed data to rejected
    if len(removed_data) > 0:
        target_rejected = np.concatenate([target_rejected, removed_data])
    
    # Combine other objects with clean target objects
    clean_combined = np.concatenate([other_data, target_clean])
    all_rejected = target_rejected
```

---

## Backward Compatibility

✅ **Fully backward compatible**
- Default: `apply_to_all=False`
- Existing scripts work unchanged
- No breaking changes to CLI or API

---

## Testing

To verify the changes:

```bash
# Test 1: Default behavior (unchanged)
filter_fits --config config.toml
# Should behave exactly as before

# Test 2: New behavior
filter_fits --config config.toml --apply-to-all
# Should filter ALL objects

# Test 3: Verify counts are different
echo "Default mode spectra count:"
fitsheader clean_data.fits | grep NAXIS2

# Then:
filter_fits --config config.toml --apply-to-all

echo "All-filter mode spectra count (should be lower):"
fitsheader clean_data.fits | grep NAXIS2
```

---

## Files Modified

1. `src/oi_zeigt/cli.py` - Lines 663-827
   - Added `--apply-to-all` option
   - Updated function signature
   - Updated CLI feedback messages
   - Updated function call with new parameter

2. `src/oi_zeigt/reduction/core.py` - Lines 360-590
   - Updated function signature
   - Updated docstring with parameter and examples
   - Replaced filtering logic with conditional based on `apply_to_all`
   - Updated variable names (`all_rejected` instead of `target_rejected`)

