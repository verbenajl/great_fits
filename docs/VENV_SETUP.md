# Virtual Environments Setup

This project uses two separate virtual environments to handle conflicting dependencies.

## Environments

### 1. `~/.venv/oi` (Main environment - Python 3.12)
**Use for**: PCA decomposition and correction
- Contains: OpenCV, NumPy 2.x, and all PCA analysis tools
- Commands available: `pca_decompose`, `pca_correct`, etc.

**Activate**:
```bash
source ~/.venv/oi/bin/activate
# or
source ~/software/oi_zeigt/activate_oi.sh
```

### 2. `~/.venv/cygrid` (Gridding environment - Python 3.12)
**Use for**: Astronomical gridding with cygrid
- Contains: cygrid, NumPy 1.26, and gridding tools
- Commands available: `map_integrated` (will use cygrid instead of scipy)
- ⚠️ WARNING: `pca_correct` NOT available in this environment

**Activate**:
```bash
source ~/.venv/cygrid/bin/activate
# or
source ~/software/oi_zeigt/activate_cygrid.sh
```

## Why Two Environments?

**Dependency Conflict**: 
- OpenCV requires NumPy 2.x
- cygrid 2.0.4 wheels are pre-built for NumPy 1.x
- Having both in one environment causes binary incompatibility errors

**Solution**: 
- Use `~/.venv/oi` for PCA workflows (with OpenCV)
- Use `~/.venv/cygrid` for gridding workflows (with cygrid)
- Most users will primarily use `~/.venv/oi`

## Typical Workflow

```bash
# For PCA analysis:
source activate_oi.sh
pca_decompose --input skychopdiff.fits --config config.toml
pca_correct --input reduced_data.fits --config config.toml --plot

# For gridding (if you need cygrid specifically):
source activate_cygrid.sh
map_integrated --config config.toml --fits reduced_data.fits

# Back to PCA analysis:
source activate_oi.sh
# ... continue with PCA workflows
```

## Error Messages

If you see:
```
ERROR: OpenCV (cv2) is required for pca_correct line detection
```

**Solution**: Make sure you're in the correct environment:
```bash
source activate_oi.sh
pca_correct ...
```

If you need cygrid for gridding instead of scipy:
```bash
source activate_cygrid.sh
map_integrated --config config.toml --fits reduced_data.fits
```

## Notes

- Both environments use Python 3.12
- The `.venv/cygrid` environment is minimal - it only has gridding dependencies
- Most work should be done in `.venv/oi` since it has OpenCV required for line detection
- If cygrid is upgraded to support NumPy 2.x, these can be merged into one environment
