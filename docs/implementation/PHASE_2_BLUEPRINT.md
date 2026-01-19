# Phase 2 Implementation Blueprint: Decomposition Module

**Status**: Ready to begin  
**Previous Phase**: Phase 1 (Foundation) ✅ COMPLETE  

## Overview

Phase 2 will implement the **Decomposition** phase of PCA analysis:
- Load spectral data from FITS files (SKYCHOPDIFF source)
- Prepare spectra (baseline subtraction, masking, normalization)
- Perform PCA decomposition
- Save components and statistics

## Phase 1 Foundation Available

### Configuration System ✅
```python
from src.oi_zeigt.pca_analysis import ConfigLoader, create_default_config_template

config = ConfigLoader('config.toml', 'missions.yml')
n_components = config.get('pca.decompose.number_components')
artifact_params = config.get_line_parameters('SOFIA_MISSION')
drop_filters = config.get_drop_filters('SOFIA_MISSION')
```

### FITS Data Loading ✅
```python
from src.oi_zeigt.pca_analysis import FITSIndexer, load_spectral_data

indexer = FITSIndexer('/path/to/fits/')
indexer.scan_directory()
sky_data = indexer.filter_by_source('SKYCHOPDIFF')

for data, header, metadata in sky_data.get_all_fits_data():
    # Process each spectrum
    pass
```

### Spectrum Preparation ✅
```python
from src.oi_zeigt.pca_analysis import (
    fit_baseline_poly, create_channel_mask, combine_masks,
    normalize_spectrum, detect_science_line_waterfall
)

# Create masks
artifact_mask = create_channel_mask(center=50, width=10, n_channels=256)
center, _ = detect_science_line_waterfall(spectrum_2d)
line_mask = create_channel_mask(center=center, width=20, n_channels=256)
combined_mask = combine_masks(artifact_mask, line_mask)

# Baseline fit
valid_mask = ~combined_mask
baseline = fit_baseline_poly(spectrum, mask=valid_mask)
corrected = spectrum - baseline

# Normalize
normalized, mean, std = normalize_spectrum(corrected)
```

### Error Handling ✅
```python
from src.oi_zeigt.pca_analysis import (
    PCAError, FITSError, LineDetectionError,
    BaselineError, DecompositionError
)

try:
    # Phase 2 operations
    pass
except DecompositionError as e:
    logger.error(f"Decomposition failed: {e}")
except PCAError as e:
    logger.error(f"PCA error: {e}")
```

## Phase 2 Module Structure

### Main Module: `decompose.py`

```python
src/oi_zeigt/pca_analysis/decompose.py
```

#### Classes to Implement

##### 1. DecompositionConfig
Configuration holder for decomposition parameters.

```python
@dataclass
class DecompositionConfig:
    """Configuration for PCA decomposition."""
    n_components: int
    pca_source: str  # 'SKYCHOPDIFF'
    normalize: bool = True
    scale: bool = False
    add_sky_diff: bool = False
    noise_cutoff: bool = False
    scramble: bool = False
```

##### 2. SpectrumPreparator
Prepare individual spectra for decomposition.

```python
class SpectrumPreparator:
    """Prepare spectra for PCA decomposition."""
    
    def __init__(self, config: ConfigLoader):
        """Initialize with configuration."""
        
    def prepare(self, spectrum: np.ndarray, 
                spectrum_2d: np.ndarray,
                mission_id: str) -> np.ndarray:
        """
        Prepare spectrum for decomposition.
        
        Steps:
        1. Detect science line from 2D spectrum
        2. Create artifact mask (from config)
        3. Create science line mask (from detection)
        4. Fit baseline to masked regions
        5. Subtract baseline
        6. Normalize
        7. Return prepared spectrum
        """
```

##### 3. DecompositionDataLoader
Load and organize all FITS data for decomposition.

```python
class DecompositionDataLoader:
    """Load and organize FITS data for decomposition."""
    
    def __init__(self, fits_dir: Path, config: ConfigLoader):
        """Initialize loader."""
        
    def load_decomposition_data(self) -> np.ndarray:
        """
        Load all SKYCHOPDIFF spectra.
        
        Returns:
        - 2D array: [n_spectra, n_channels]
        - Properly prepared (baseline, mask, normalized)
        """
        
    def load_reference_spectra(self) -> Tuple[np.ndarray, List[dict]]:
        """Load SKYCHOPDIFF reference spectra with metadata."""
```

##### 4. PCADecomposer
Perform PCA decomposition on prepared data.

```python
class PCADecomposer:
    """Perform PCA decomposition on spectral data."""
    
    def __init__(self, n_components: int, scale: bool = False):
        """Initialize decomposer."""
        
    def fit(self, spectral_matrix: np.ndarray) -> 'PCADecomposer':
        """
        Fit PCA to spectral matrix.
        
        Input:
        - spectral_matrix: [n_spectra, n_channels]
        
        Stores:
        - self.pca_model (sklearn.decomposition.PCA)
        - self.components (fitted components)
        - self.explained_variance_ratio
        """
        
    def get_components(self) -> np.ndarray:
        """Get PCA components: [n_components, n_channels]"""
        
    def get_explained_variance(self) -> np.ndarray:
        """Get explained variance ratio per component."""
        
    def transform(self, spectral_matrix: np.ndarray) -> np.ndarray:
        """Project data onto components: [n_spectra, n_components]"""
        
    def reconstruct(self, coefficients: np.ndarray) -> np.ndarray:
        """Reconstruct from components."""
```

##### 5. DecompositionResult
Store decomposition results for later use.

```python
@dataclass
class DecompositionResult:
    """Results from PCA decomposition."""
    components: np.ndarray  # [n_components, n_channels]
    explained_variance_ratio: np.ndarray  # [n_components]
    explained_variance: np.ndarray  # [n_components]
    pca_model: sklearn.decomposition.PCA
    config: DecompositionConfig
    metadata: dict  # mission, date, tag, etc.
    
    def save(self, filepath: Path) -> None:
        """Save decomposition results to pickle."""
        
    @classmethod
    def load(cls, filepath: Path) -> 'DecompositionResult':
        """Load decomposition results from pickle."""
```

#### Main Function

```python
def decompose_spectra(fits_directory: Path, 
                     config_file: Path,
                     mission_file: Path) -> DecompositionResult:
    """
    Perform PCA decomposition on FITS spectra.
    
    Steps:
    1. Load configuration
    2. Index FITS files
    3. Load SKYCHOPDIFF spectra
    4. Prepare each spectrum (baseline, mask, normalize)
    5. Create spectral matrix
    6. Fit PCA
    7. Return results
    
    Raises:
    - FITSError: If FITS files can't be read
    - ConfigurationError: If config is invalid
    - LineDetectionError: If science line detection fails
    - DecompositionError: If PCA fitting fails
    """
```

## Data Flow in Phase 2

```
Config File (TOML)
    │
    ▼
ConfigLoader ──────┐
                   │
FITS Directory     │
    │              │
    ▼              │
FITSIndexer ───┐   │
               │   │
               ▼   ▼
           DecompositionDataLoader
               │
               ▼
           Load & Filter:
           - Source = 'SKYCHOPDIFF'
           - Apply drop_filters
               │
               ▼
           [n_spectra, multiple_readouts, n_channels]
               │
               ├──────────────────┐
               │                  │
               ▼                  ▼
         Detect Science Line  Create Artifact Mask
         (per spectrum_2d)    (from config)
               │                  │
               └──────────┬───────┘
                          ▼
                   Combine Masks
                   (artifact + line)
                          │
                          ▼
                   SpectrumPreparator
                   (for each spectrum)
                   ├─ fit baseline
                   ├─ subtract baseline
                   ├─ apply masks
                   └─ normalize
                          │
                          ▼
                   Spectral Matrix
                   [n_spectra, n_channels]
                          │
                          ▼
                    PCADecomposer
                    ├─ fit PCA
                    ├─ extract components
                    └─ calculate variance
                          │
                          ▼
                   DecompositionResult
                   ├─ components
                   ├─ explained_variance
                   ├─ config metadata
                   └─ save to file
```

## Implementation Checklist

### Step 1: Data Loading (DecompositionDataLoader)
- [ ] Initialize with FITS directory and config
- [ ] Scan FITS files
- [ ] Filter by source (SKYCHOPDIFF)
- [ ] Apply drop filters (exclude specified missions/telescopes/scans)
- [ ] Organize into 2D array with metadata
- [ ] Handle missing/invalid data gracefully

### Step 2: Spectrum Preparation (SpectrumPreparator)
- [ ] For each spectrum:
  - [ ] Detect science line from 2D data
  - [ ] Create artifact mask
  - [ ] Create science line mask
  - [ ] Combine masks
  - [ ] Fit polynomial baseline
  - [ ] Subtract baseline
  - [ ] Apply mask (set masked regions to zero/NaN)
  - [ ] Normalize

### Step 3: PCA Decomposition (PCADecomposer)
- [ ] Create 2D spectral matrix [n_spectra, n_channels]
- [ ] Initialize sklearn PCA with n_components
- [ ] Fit PCA model
- [ ] Extract components
- [ ] Extract explained variance
- [ ] Provide transform/reconstruct methods

### Step 4: Results Management (DecompositionResult)
- [ ] Store components and variance
- [ ] Store configuration and metadata
- [ ] Implement pickle save/load
- [ ] Generate summary statistics

### Step 5: Integration & Testing
- [ ] Create unit tests for each class
- [ ] Create integration test (full decomposition pipeline)
- [ ] Add logging at key steps
- [ ] Handle error cases
- [ ] Update `__init__.py` with new exports

## Error Handling

Use these exceptions from Phase 1:
- `DecompositionError` - Raised by PCADecomposer if fitting fails
- `FITSError` - If FITS files can't be read
- `ConfigurationError` - If config is invalid
- `LineDetectionError` - If science line detection fails
- `BaselineError` - If baseline fitting fails
- `PCAError` - Catch all (parent class)

## Testing Strategy

### Unit Tests
```python
# test_decomposition.py

def test_spectrum_preparator():
    """Test single spectrum preparation."""
    
def test_decomposition_data_loader():
    """Test data loading and filtering."""
    
def test_pca_decomposer():
    """Test PCA fitting and component extraction."""
    
def test_decomposition_result():
    """Test result saving/loading."""

def test_decompose_spectra_integration():
    """Test full decomposition pipeline."""
```

### Test Data Strategy
- Create mock FITS files with known spectra
- Use synthetic data (simple Gaussian lines)
- Test with actual FITS files if available
- Verify components make sense (narrow, positive contributions)

## Expected Output

After Phase 2 completion:

```python
from src.oi_zeigt.pca_analysis import decompose_spectra

result = decompose_spectra(
    fits_directory=Path('/data/FITS'),
    config_file=Path('config.toml'),
    mission_file=Path('missions.yml')
)

print(f"Components shape: {result.components.shape}")  # [5, 256]
print(f"Explained variance: {result.explained_variance_ratio}")  # [0.45, 0.25, 0.15, ...]
print(f"Total variance explained: {result.explained_variance_ratio.sum():.1%}")  # 90%+

# Save for Phase 3
result.save(Path('decomposition_result.pkl'))
```

## Next: Phase 3 (Correction)

Phase 3 will:
1. Load decomposition results from Phase 2
2. Load M51CENTER (science) spectra
3. Project science spectra onto components
4. Calculate optimal correction (least-squares fit)
5. Apply correction to remove receiver variations
6. Preserve science line (via masking)

## Timeline

Estimated effort for Phase 2:
- Core implementation: 4-6 hours
- Testing & debugging: 2-3 hours
- Documentation: 1-2 hours
- **Total**: ~1 week of part-time work

Ready to begin Phase 2 implementation whenever you're ready!
