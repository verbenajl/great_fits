# Original vs. Current: How Telluric Masking Works

## Direct Comparison

### Original Scripts (`pca_decompose.py` + `pca_utilities.py`)

**How telluric masking is handled:**

1. **Parameter Source**: 
   - Manual `--pca_exclude_range` command-line option
   - Must be explicitly specified by user
   - Example: `--pca_exclude_range -5 5` (velocity range in km/s)

2. **Code Location** (pca_utilities.py):
   ```python
   def apply_pca_exclude_range(spectrum, config):
       """Apply exclusion range AFTER spectrum loading"""
       if config["pca_exclude_range"]:
           channels = get_exclude_channels(config)
           # Remove those channels from spectrum
           spectrum = ma.concatenate(
               (spectrum[0:channels[0] - 1], spectrum[channels[1] - 1:]))
       return spectrum
   ```

3. **Usage** (pca_decompose.py, line 310):
   ```python
   spectrum = apply_pca_exclude_range(spectrum, config)
   ```

4. **Problem**: 
   - Telluric parameters in YAML are **NOT automatically used**
   - User must manually calculate velocity range and pass as argument
   - Tedious and error-prone
   - Different for each mission but not automatic

---

### Current Implementation (`src/oi_zeigt/pca_analysis/decompose.py`)

**How telluric masking is now handled:**

1. **Parameter Source**: 
   - Automatically from `mission_id_parameters.yml` per mission
   - Mission detected from FITS file
   - No user input needed

2. **Code Location** (decompose.py):

   **Step 1: Load during initialization** (lines 98-100)
   ```python
   line_params = config.get_line_parameters(mission_id)
   self.artifact_center = line_params.get('telluric_line_center', 50)
   self.artifact_width = line_params.get('telluric_line_width', 10)
   ```

   **Step 2: Create mask** (lines 172-174)
   ```python
   artifact_mask = create_artifact_mask(
       n_channels,
       artifact_center=self.artifact_center,
       artifact_width=self.artifact_width
   )
   ```

   **Step 3: Combine with science line mask** (lines 176-181)
   ```python
   line_mask = create_science_line_mask(spectrum, center=line_center, width=20)
   combined_mask = combine_masks(artifact_mask, line_mask)
   metadata['masked_channels'] = int(np.sum(combined_mask))
   ```

   **Step 4: Use for baseline fitting** (lines 185-196)
   ```python
   valid_mask = ~combined_mask
   baseline = fit_baseline_poly(spectrum, order=self.baseline_order, mask=valid_mask)
   ```

   **Step 5: Apply to spectrum** (lines 199-201)
   ```python
   spectrum = spectrum - baseline
   spectrum[combined_mask] = 0.0  # Zero out both telluric + science line
   ```

3. **Advantage**: 
   - ✅ Fully automatic per mission
   - ✅ No user configuration needed
   - ✅ Both telluric AND science line protected
   - ✅ Better baseline fitting (unmasked channels only)
   - ✅ Metadata tracking

---

## Side-by-Side Comparison

| Aspect | Original | Current |
|--------|----------|---------|
| **Telluric params source** | Manual argument | YAML file (mission_id_parameters.yml) |
| **When is it configured?** | Per command run | Per mission (automatic) |
| **User interaction** | Must specify `--pca_exclude_range` | None needed (automatic) |
| **Mission awareness** | No | Yes (auto-detect from FITS) |
| **Science line handling** | Not considered | Also masked + protected |
| **Baseline fitting** | May include contaminated channels | Only uses clean channels |
| **Metadata** | None recorded | Tracks masked channels |
| **Requires manual calculation?** | Yes (velocity to channels) | No (automatic from YAML) |
| **Different missions?** | Need different parameters | Handled automatically per mission |
| **Likelihood of error** | High (manual calculation) | Low (automatic) |

---

## Why the Current Approach is Better

### Problem 1: Manual Calculation Needed

**Original approach**:
```bash
# User must figure out:
# 1. What velocity range has telluric?
# 2. Convert to channel range
# 3. Pass as argument

pca_decompose.py --input data.fits --pca_exclude_range -5 5
```

**Current approach**:
```yaml
# Already in YAML file for your mission
2017-02-01_GR_F367:
  telluric_line_center: 592
  telluric_line_width: 30
```

```bash
# Just run normally, masking is automatic
pca_decompose --config config.toml
```

### Problem 2: Science Line Not Protected

**Original approach**: Only telluric masking
```python
spectrum = apply_pca_exclude_range(spectrum, config)
# Excludes telluric region
# But science line still in spectrum → baseline biased
```

**Current approach**: Telluric + science line protection
```python
# Both are masked
artifact_mask = create_artifact_mask(n, 592, 30)      # Telluric
line_mask = create_science_line_mask(spectrum, 450)   # Science
combined = combine_masks(artifact_mask, line_mask)

# Baseline fitted ONLY on unmasked regions
baseline = fit_baseline_poly(spectrum, mask=~combined)

# Both regions zeroed
spectrum[combined] = 0.0
```

### Problem 3: No Mission Awareness

**Original approach**: Same parameters for all missions
```bash
pca_decompose.py --input 2016_data.fits --pca_exclude_range -5 5
pca_decompose.py --input 2017_data.fits --pca_exclude_range -5 5  # WRONG!
# (Different missions have different telluric centers)
```

**Current approach**: Mission-specific automatically
```python
# Load mission from FITS header
mission_id = extract_mission_id_from_fits(fits_file)

# Get mission-specific params
line_params = config.get_line_parameters(mission_id)
# → automatically uses correct telluric center for this mission!
```

---

## The Telluric Parameters in mission_id_parameters.yml

From your file, each mission has:

```yaml
2016-05-12_GR_F296:
  telluric_line_center: 566      # Channel number
  telluric_line_width: 30        # Full width in channels

2017-02-01_GR_F367:
  telluric_line_center: 592      # Different for this mission!
  telluric_line_width: 30
  
2017-06-07_GR_F401:
  telluric_line_center: 564      # Different again!
  telluric_line_width: 30
```

**Key insight**: The telluric line position **varies by mission** depending on atmospheric conditions. The current approach handles this automatically.

---

## What the YAML Parameters Mean

From the code (line_detection.py):

```python
def create_artifact_mask(n_channels: int, artifact_center: int,
                        artifact_width: int) -> np.ndarray:
    """Create boolean mask for instrumental artifact region."""
    mask = np.zeros(n_channels, dtype=bool)
    
    # Center ± width/2
    ch_start = max(0, artifact_center - artifact_width // 2)
    ch_end = min(n_channels, artifact_center + artifact_width // 2 + 1)
    
    # Everything in this range is masked
    mask[ch_start:ch_end] = True
    
    return mask
```

**Example** (2017-02-01_GR_F367):
- `center = 592`, `width = 30`
- Masked region: `592 - 15` to `592 + 15` = channels 577-607
- These channels (telluric absorption) are zeroed before PCA

---

## Data Flow Comparison

### Original Implementation

```
Raw Spectrum (1264 channels)
           ↓
    apply_pca_exclude_range()
    (removes specified channels)
           ↓
    Shortened Spectrum (~1250 channels)
           ↓
    PCA Decomposition
           ↓
    Components (may still be affected by baseline bias)
```

### Current Implementation

```
Raw Spectrum (1264 channels)
           ↓
    ┌─ Create artifact_mask (telluric)
    │
    ├─ Detect science line
    │
    ├─ Create science_line_mask
    │
    └─ Combine masks (protection for both)
           ↓
    Fit baseline (ONLY on unmasked channels)
           ↓
    Subtract baseline
           ↓
    Zero both artifact + science regions
           ↓
    Normalize
           ↓
    Prepared Spectrum (1264 channels, masked regions = 0)
           ↓
    PCA Decomposition
           ↓
    Clean Components (no telluric or baseline bias)
```

---

## Key Advantages of Current Implementation

| Advantage | Benefit |
|-----------|---------|
| **Automatic mission awareness** | Correct masks for each mission without user error |
| **Science line protection** | Preserves emission line for analysis while masking for baseline |
| **Better baseline fitting** | Uses only unmasked regions → more accurate baseline |
| **Full spectrum retained** | Keeps all 1264 channels (masked = 0, not removed) |
| **Metadata tracking** | Records what was masked for quality control |
| **No manual configuration** | Just run the command, masking happens automatically |
| **Mission-specific** | Different telluric centers handled correctly per mission |

---

## Verification: Is This Actually Happening?

### Check 1: Parameters are loaded
```python
# In decompose.py, lines 98-100
line_params = config.get_line_parameters(mission_id)
self.artifact_center = line_params.get('telluric_line_center', 50)
self.artifact_width = line_params.get('telluric_line_width', 10)
# ✓ CONFIRMED: Parameters are loaded
```

### Check 2: Mask is created
```python
# In decompose.py, lines 172-174
artifact_mask = create_artifact_mask(
    n_channels,
    artifact_center=self.artifact_center,
    artifact_width=self.artifact_width
)
# ✓ CONFIRMED: Mask is created with telluric params
```

### Check 3: Mask is applied
```python
# In decompose.py, lines 176-182
combined_mask = combine_masks(artifact_mask, line_mask)

# In decompose.py, lines 185-196
baseline = fit_baseline_poly(spectrum, order=self.baseline_order, mask=~combined_mask)

# In decompose.py, lines 199-201
spectrum[combined_mask] = 0.0
# ✓ CONFIRMED: Mask is actively used
```

---

## Conclusion

**Original Scripts**: 
- Telluric parameters in YAML are **defined but not used**
- Masking is **manual via command-line argument**
- User must **calculate channel ranges themselves**
- Science line is **not protected during baseline fitting**

**Current Implementation**:
- Telluric parameters in YAML are **automatically loaded and applied**
- Masking is **per-mission, with no user configuration needed**
- Channels are **calculated automatically from mission parameters**
- Science line is **also protected during baseline fitting**
- **Result: Cleaner PCA components** ✓

The move from manual to automatic telluric masking is a significant improvement in code robustness and user experience.
