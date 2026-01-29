# Channel Trimming and PCA Parameters: Per-Mission vs. Original Code

## Summary

✅ **Yes, both approaches trim to the valid channel range**  
✅ **Yes, they use the same PCA parameters**  
✅ **This ensures consistency and will correlate properly with science spectra**

---

## 1. Channel Trimming

### Original Code (decompose.py, lines 800-804)

```python
# Find valid channels
valid_channels = ~np.all(np.isnan(spectra), axis=0)
first_valid = np.where(valid_channels)[0][0]
last_valid = np.where(valid_channels)[0][-1]
spectra = spectra[:, first_valid:last_valid+1]
```

### Per-Mission Script (pca_decompose_per_mission.py, lines 154-157)

```python
# Find valid channels (not all NaN)
valid_channels = ~np.all(np.isnan(spectra), axis=0)
first_valid = np.where(valid_channels)[0][0]
last_valid = np.where(valid_channels)[0][-1]
spectra = spectra[:, first_valid:last_valid+1]
```

### What This Does

**Identifies which channels are "valid":**
- A channel is invalid if it contains ALL NaN values across the entire spectrum set
- Creates a boolean mask: `True` = at least one non-NaN value, `False` = all NaN

**Extracts the continuous range:**
- Finds the first and last non-NaN channel index
- Trims all spectra to this range: `spectra[:, first_valid:last_valid+1]`

**Why trim?**
- FITS files have 1264 channels total, but only 1228-1229 are valid spectral data
- Channels before the first valid or after the last valid contain all NaN/zero values
- Trimming removes padding while preserving the actual spectral information

### Example from Your Data

```
Original channels:     1264 (full FITS array)
After trimming:        1229 channels (actual spectral range)
Why?:                  First 17-18 and last 16-17 channels are all NaN
```

---

## 2. NaN Filling (After Trimming)

Both implementations use **identical post-trimming logic:**

```python
# Fill remaining NaNs with per-channel mean
for ch in range(spectra.shape[1]):
    channel_data = spectra[:, ch]
    if np.any(np.isnan(channel_data)):
        mean_val = np.nanmean(channel_data)
        spectra[np.isnan(channel_data), ch] = mean_val
```

**What happens:**
- After trimming, individual channels may still have scattered NaN values
- These are filled with the mean of non-NaN values in that channel
- This preserves channel continuity while handling sparse missing data

---

## 3. PCA Parameters

### Original Code (decompose.py, lines 429-436)

```python
from sklearn.decomposition import PCA

self.pca_model = PCA(
    n_components=self.n_components,
    whiten=self.scale,
    random_state=42
)

self.pca_model.fit(spectral_matrix)
```

### Per-Mission Script (pca_decompose_per_mission.py, uses PCADecomposer class)

Same PCA initialization via the `PCADecomposer` class:

```python
decomposer = PCADecomposer(n_components=n_components, scale=False)
decomposer.fit(spectra)
```

Which internally calls (from PCADecomposer.fit() method):

```python
self.pca_model = PCA(
    n_components=self.n_components,
    whiten=self.scale,
    random_state=42
)

self.pca_model.fit(spectral_matrix)
```

### Parameter Explanation

| Parameter | Your Value | Meaning |
|-----------|-----------|---------|
| `n_components` | 5 | Extract 5 principal components |
| `whiten` (AKA `scale`) | `False` | Do NOT standardize features |
| `random_state` | 42 | Fixed seed for reproducibility |

### What These Mean

#### `whiten=False` (scale=False)
- **NOT standardizing** the spectral data
- Each component preserves the original variance scale
- Appropriate for spectral data where absolute intensity matters
- You're keeping components in their natural physical units

#### `random_state=42`
- **Reproducibility**: Same results every time you run
- Not randomness in the algorithm (PCA is deterministic)
- Only affects initialization if you had `svd_solver='randomized'` (you don't)
- Acts as a synchronization mechanism if code ever changes

---

## 4. Consistency Check: Will Science Spectra Correlate Properly?

### YES - Here's Why

When you apply the PCA components to science spectra:

#### Step 1: Load Science Spectrum
```python
science_spectrum = fits_data[spectrum_index]  # 1264 channels
```

#### Step 2: Apply Same Trimming
```python
# MUST use same first_valid and last_valid from reference decomposition!
science_spectrum_trimmed = science_spectrum[first_valid:last_valid+1]
# Now it's 1229 channels, matching reference decomposition
```

#### Step 3: Fill NaNs (Optional)
```python
# Use same per-channel mean method as reference
```

#### Step 4: Project onto Components
```python
# Components are [5, 1229]
# science_spectrum_trimmed is [1229]
# Dot product: [5] coefficients
coefficients = decomposer.pca_model.transform(science_spectrum_trimmed.reshape(1, -1))
```

### Critical Implementation Notes

**You MUST preserve channel indices:**

```python
# Load reference decomposition
with open('decomposition_2016-05-12_GR_F_20160512_components.pkl', 'rb') as f:
    result = pickle.load(f)

# Store the channel range used during decomposition
first_valid = result.metadata.get('first_valid_channel', None)
last_valid = result.metadata.get('last_valid_channel', None)
# OR: Recalculate from metadata if not stored
```

**Recommendation:** Modify `pca_decompose_per_mission.py` to **store the channel range**:

```python
metadata={
    'n_reference_spectra': len(spectra),
    'n_channels': spectra.shape[1],
    'source': 'SKYCHOPDIFF',
    'mission_id': mission_id,
    'first_valid_channel': first_valid,    # ADD THIS
    'last_valid_channel': last_valid,      # ADD THIS
}
```

---

## 5. Summary Table

| Aspect | Original Code | Per-Mission Script | Consistency |
|--------|--------------|-------------------|-------------|
| **Channel Trimming** | Lines 800-804 | Lines 154-157 | ✅ Identical |
| **Trimming Method** | Find first/last non-NaN | Find first/last non-NaN | ✅ Identical |
| **NaN Filling** | Per-channel mean | Per-channel mean | ✅ Identical |
| **PCA Algorithm** | sklearn.decomposition.PCA | Same (via PCADecomposer) | ✅ Identical |
| **whiten Parameter** | False | False | ✅ Identical |
| **random_state** | 42 | 42 | ✅ Identical |
| **Result** | 1228-1229 channels | 1228-1229 channels | ✅ Identical |

---

## 6. Recommendation for Phase 4 (Science Spectra Correction)

When implementing the correction phase for science spectra:

### Option A: Store Channel Range (Recommended)
```python
# In pca_decompose_per_mission.py, modify metadata:
metadata={
    'n_reference_spectra': len(spectra),
    'n_channels': spectra.shape[1],
    'source': 'SKYCHOPDIFF',
    'mission_id': mission_id,
    'first_valid_channel': first_valid,
    'last_valid_channel': last_valid,
    'original_n_channels': 1264,  # before trimming
}

# Then in correction script, automatically apply same trim:
first_valid = result['metadata']['first_valid_channel']
last_valid = result['metadata']['last_valid_channel']
science_spectrum_trimmed = science_spectrum[first_valid:last_valid+1]
```

### Option B: Recalculate From Reference Data
```python
# Load all SKYCHOPDIFF spectra for the mission
# Recalculate valid_channels and first_valid/last_valid
# This ensures consistency even if metadata is lost
```

### Option C: Use Full Channel Range
```python
# Load decomposition result
# Check n_channels in metadata
# Trim science spectra to this size
# Less explicit but relies on metadata consistency
```

**Recommendation:** Use **Option A** - it's explicit, documented, and foolproof.

---

## 7. Impact on Your Data

Your four missions will decompose with these characteristics:

```
2016-05-12_GR_F:
  - Reference: 2,268 SKYCHOPDIFF spectra → 1229 channels
  - Science spectra: trim [first_valid:last_valid+1] → 1229 channels
  - Correlation: ✅ Perfect alignment

2016-05-18_GR_F:
  - Reference: 6,426 SKYCHOPDIFF spectra → 1229 channels
  - Science spectra: trim [first_valid:last_valid+1] → 1229 channels
  - Correlation: ✅ Perfect alignment

2016-05-19_GR_F:
  - Reference: 3,024 SKYCHOPDIFF spectra → 1229 channels
  - Science spectra: trim [first_valid:last_valid+1] → 1229 channels
  - Correlation: ✅ Perfect alignment

2017-02-01_GR_F:
  - Reference: 13,986 SKYCHOPDIFF spectra → 1229 channels
  - Science spectra: trim [first_valid:last_valid+1] → 1229 channels
  - Correlation: ✅ Perfect alignment
```

All science spectra will have exactly the same channel indices as the reference spectra, so projection will work perfectly.

