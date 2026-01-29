# filter_fits - Practical Examples & Workflows

## Quick Reference

```bash
# Basic usage (uses config)
filter_fits --config config.toml

# Custom threshold
filter_fits --config config.toml --nan-threshold 0.50

# Different object
filter_fits --config config.toml --object NGC1234

# Remove bad observations
filter_fits --config config.toml --remove AOR_ID \
    --remove-values bad_id_1 --remove-values bad_id_2

# Direct FITS file
filter_fits --fits /path/to/data.fits --object M51
```

---

## Real-World Workflows

### Workflow 1: Initial Data Quality Check

**Goal**: Separate bad spectra from good ones, review both

**Steps:**

```bash
# 1. Filter with default settings
filter_fits --config config.toml

# Output:
# clean_data.fits: 14,504 records (good spectra)
# rejected_data.fits: 2,150 records (bad spectra)
```

**Then:**
```bash
# 2. Inspect rejected spectra
plot_sample_spectra --fits rejected_data.fits \
    --object M51 --num-spectra 20

# 3. Review bad spectra visually
# Are they really bad? Or can threshold be adjusted?
```

**Decision:**
- If too many rejections: Lower threshold
- If too few rejections: Raise threshold
- If pattern visible: Use removal filtering

---

### Workflow 2: Remove Problematic Observations

**Goal**: Known bad observations in data - remove them

**Setup:**
First, identify bad observations:
```bash
# List unique AOR_ID values
fitsheader clean_data.fits | grep AOR_ID | sort | uniq
```

Output might show:
```
04_0116_0020609  <- Known bad
04_0116_0020506  <- Known bad
04_0116_0020507  <- Good
...
```

**Filter command:**
```bash
filter_fits --config config.toml \
    --remove AOR_ID \
    --remove-values 04_0116_0020609 \
    --remove-values 04_0116_0020506
```

**Result:**
```
Clean file: Clean data WITHOUT bad AOR_IDs + filtered by NaN
Rejected file: Bad AOR_IDs + high-NaN spectra
```

**Verify:**
```bash
# Check AOR_ID removed
fitsheader clean_data.fits | grep AOR_ID | grep "04_0116_0020609"
# Should return nothing
```

---

### Workflow 3: Strict Quality Filter for Publication

**Goal**: Get only the best spectra for analysis

**Commands (progressive):**

```bash
# Step 1: Stricter NaN threshold
filter_fits --config config.toml --nan-threshold 0.05
# Keep only spectra with < 5% NaNs

# Step 2: Remove known issues
filter_fits --config config.toml --nan-threshold 0.05 \
    --remove AOR_ID \
    --remove-values 04_0116_0020609 \
    --remove-values 04_0116_0020506

# Step 3: Check results
echo "Clean spectra:"
fitsheader clean_data.fits | grep NAXIS2
# Shows number of records

echo "Rejected spectra:"
fitsheader rejected_data.fits | grep NAXIS2
```

**Expected outcome:**
- ~70% of spectra kept (stricter threshold)
- High-quality dataset for publication
- Clear audit trail in rejected file

---

### Workflow 4: Compare Thresholds

**Goal**: Find optimal NaN threshold

**Script:**
```bash
#!/bin/bash

for threshold in 0.10 0.15 0.20 0.25 0.30; do
    echo "Testing threshold: $threshold"
    
    filter_fits --config config.toml \
        --nan-threshold $threshold \
        --output-clean clean_${threshold}.fits \
        --output-rejected rejected_${threshold}.fits
    
    n_clean=$(fitsheader clean_${threshold}.fits | grep NAXIS2 | awk '{print $NF}')
    n_rejected=$(fitsheader rejected_${threshold}.fits | grep NAXIS2 | awk '{print $NF}')
    
    printf "  Clean: %6d  |  Rejected: %6d  |  Ratio: %.1f%%\n" \
        $n_clean $n_rejected $(echo "scale=1; 100*$n_clean/($n_clean+$n_rejected)" | bc)
done
```

**Output:**
```
Testing threshold: 0.10
  Clean:  10234  |  Rejected:  6420  |  Ratio: 61.4%

Testing threshold: 0.15
  Clean:  11987  |  Rejected:  4667  |  Ratio: 72.0%

Testing threshold: 0.20
  Clean:  13456  |  Rejected:  3198  |  Ratio: 80.8%

Testing threshold: 0.25
  Clean:  14102  |  Rejected:  2552  |  Ratio: 84.7%

Testing threshold: 0.30
  Clean:  14504  |  Rejected:  2150  |  Ratio: 87.1%
```

**Decision process:**
1. Plot graph of rejected spectra vs threshold
2. Look for "knee" in curve
3. Choose threshold at knee
4. Use for final analysis

---

### Workflow 5: Multi-Object Filtering

**Goal**: Filter different objects with different thresholds

**Setup (in config):**
```toml
[parameters]
object = "M51"

[filters]
blank_fraction = 0.20
```

**Commands:**
```bash
# Filter M51 with default threshold
filter_fits --config config.toml --object M51 \
    --output-clean m51_clean.fits \
    --output-rejected m51_rejected.fits

# Filter NGC1234 with stricter threshold
filter_fits --config config.toml --object NGC1234 \
    --nan-threshold 0.10 \
    --output-clean ngc_clean.fits \
    --output-rejected ngc_rejected.fits

# Filter other objects - keep all
filter_fits --config config.toml --object OTHER \
    --nan-threshold 1.0 \
    --output-clean other_clean.fits \
    --output-rejected other_rejected.fits
```

**Result:**
- M51: 80% kept (20% NaN threshold)
- NGC1234: 65% kept (10% NaN threshold)
- Others: 100% kept (very permissive)

---

### Workflow 6: Iterative Refinement

**Goal**: Start rough, get progressively more selective

**Step 1: Overview**
```bash
# See what's available
filter_fits --config config.toml --nan-threshold 0.50
# Result: Keep most data, identify worst spectra
```

**Step 2: Check rejected data**
```bash
plot_sample_spectra --fits rejected_data.fits \
    --object M51 --num-spectra 30 --output rejected_overview.png

# Manual inspection: Are rejections valid?
```

**Step 3: Adjust if needed**
```bash
# If too permissive:
filter_fits --config config.toml --nan-threshold 0.25

# If too strict:
filter_fits --config config.toml --nan-threshold 0.35
```

**Step 4: Remove known issues**
```bash
# Inspect rejected_data.fits for patterns
# Identify problematic AOR_IDs

filter_fits --config config.toml \
    --nan-threshold 0.25 \
    --remove AOR_ID \
    --remove-values problem_id_1 \
    --remove-values problem_id_2
```

**Step 5: Final verification**
```bash
# Plot final clean data
plot_sample_spectra --fits clean_data.fits \
    --object M51 --num-spectra 50 --output final_review.png

# Check statistics
fitsheader clean_data.fits | grep NAXIS2
# Verify count meets expectations
```

---

## Configuration Examples

### Config 1: Conservative (High Quality)
```toml
[filters]
blank_fraction = 0.10  # Only 10% NaNs allowed

[reduction]
baseline = 3
window = [450, 500]
extract = [350, 700]
```

**Result:** ~60% of spectra kept, very high quality

### Config 2: Moderate (Balanced)
```toml
[filters]
blank_fraction = 0.20  # 20% NaNs allowed

[reduction]
baseline = 3
window = [450, 500]
extract = [350, 700]
```

**Result:** ~80% of spectra kept, good quality

### Config 3: Permissive (Include More Data)
```toml
[filters]
blank_fraction = 0.50  # 50% NaNs allowed

[reduction]
baseline = 3
window = [450, 500]
extract = [350, 700]
```

**Result:** ~95% of spectra kept, includes marginal data

---

## Common Issues & Solutions

### Issue 1: Too Many Spectra Rejected

**Problem:** Only 20% of data kept after filtering

**Solutions:**
1. **Lower NaN threshold**
   ```bash
   filter_fits --config config.toml --nan-threshold 0.30
   ```

2. **Check if target object name is correct**
   ```bash
   # See what objects are in FITS
   fitsheader data.fits | grep OBJECT | head -20
   
   # Update config or use --object flag
   ```

3. **Visualize rejected spectra**
   ```bash
   plot_sample_spectra --fits rejected_data.fits --num-spectra 20
   # Are they really bad or just different?
   ```

### Issue 2: Too Many Spectra Kept

**Problem:** 99% of data kept, no filtering happening

**Solutions:**
1. **Raise NaN threshold**
   ```bash
   filter_fits --config config.toml --nan-threshold 0.05
   ```

2. **Remove problematic observations**
   ```bash
   # Identify bad AOR_IDs first
   filter_fits --config config.toml \
       --remove AOR_ID \
       --remove-values bad_id_1 --remove-values bad_id_2
   ```

### Issue 3: Output Files Not Where Expected

**Problem:** Files saved to current directory instead of config path

**Solution:** Explicitly specify output paths
```bash
filter_fits --config config.toml \
    --output-clean /path/to/clean_data.fits \
    --output-rejected /path/to/rejected_data.fits
```

### Issue 4: "Column not found" Error

**Problem:** `--remove AOR_ID` fails with column not found

**Solutions:**
1. **Check column exists**
   ```bash
   fitsheader data.fits | grep TFIELDS  # Number of columns
   fitsheader data.fits | grep TTYPE     # Column names
   ```

2. **Use correct column name**
   ```bash
   filter_fits --config config.toml --remove OBSID
   # Instead of AOR_ID if that's the actual name
   ```

---

## Advanced Patterns

### Pattern 1: Batch Processing Multiple Files

```bash
#!/bin/bash

for file in /data/observations/*.fits; do
    filename=$(basename "$file")
    echo "Processing: $filename"
    
    filter_fits --fits "$file" \
        --object M51 \
        --nan-threshold 0.20 \
        --output-clean "clean_$filename" \
        --output-rejected "rejected_$filename"
done
```

### Pattern 2: Conditional Removal

```bash
#!/bin/bash

# Only remove if file has been flagged
if grep -q "bad_aor" bad_observations.txt; then
    bad_ids=$(grep "bad_aor" bad_observations.txt | awk '{print $2}')
    
    # Build remove-values arguments
    remove_args=""
    for id in $bad_ids; do
        remove_args="$remove_args --remove-values $id"
    done
    
    filter_fits --config config.toml \
        --remove AOR_ID $remove_args
else
    filter_fits --config config.toml
fi
```

### Pattern 3: Quality Metrics

```bash
#!/bin/bash

echo "Filtering Statistics"
echo "===================="

filter_fits --config config.toml --nan-threshold 0.20

clean_records=$(fitsheader clean_data.fits | grep NAXIS2 | awk '{print $NF}')
rejected_records=$(fitsheader rejected_data.fits | grep NAXIS2 | awk '{print $NF}')
total=$((clean_records + rejected_records))

echo "Clean records:     $clean_records"
echo "Rejected records:  $rejected_records"
echo "Total records:     $total"
echo "Keep ratio:        $(echo "scale=1; 100*$clean_records/$total" | bc)%"
echo "Rejection ratio:   $(echo "scale=1; 100*$rejected_records/$total" | bc)%"
```

---

## Integration with Other Commands

### Workflow: Filter → Reduce → Map

```bash
# 1. Filter bad spectra
filter_fits --config config.toml

# 2. Reduce clean spectra
reduce_spectra --config config.toml \
    --fits clean_data.fits \
    --baseline

# 3. Map to datacube
average_spectra --config config.toml \
    --fits reduced_data.fits
```

### Workflow: Filter → Analyze → Statistics

```bash
# 1. Get quality spectra
filter_fits --config config.toml --nan-threshold 0.10

# 2. Plot samples
plot_sample_spectra --fits clean_data.fits \
    --object M51 --num-spectra 50

# 3. Extract statistics
# (Custom script to analyze clean_data.fits)
```

---

## Tips & Best Practices

### ✓ DO:
- Start with default threshold, then adjust
- Review rejected spectra visually
- Use removal filtering for known bad observations
- Keep both clean and rejected files for audit trail
- Document your filtering choices

### ✗ DON'T:
- Use very strict threshold (< 5%) without reviewing rejected data
- Remove observations without investigation
- Assume default threshold is optimal
- Delete rejected files immediately
- Change thresholds between different object groups inconsistently

---

## See Also

- `FILTER_FITS_COMMAND.md` - Command reference
- `FILTER_FITS_CODE_WALKTHROUGH.md` - Implementation details
- `plot_sample_spectra` - Visualize filtered data
- `reduce_spectra` - Process filtered data

