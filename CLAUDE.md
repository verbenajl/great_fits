# CLAUDE.md — oi_zeigt project instructions

This file is read automatically at the start of every Claude Code session.
It captures project conventions, constraints, and context so we don't
re-derive them each time.

---

## What this project is

A Python pipeline for processing SOFIA GREAT [OI] 63 µm spectral data.
Raw data is a large FITS binary table (one row = one spectrum + metadata).
The pipeline filters, reduces, PCA-corrects, and grids the spectra into maps.
See `CODEBASE.md` for a full reference.

---

## Virtualenvs — which to activate

Two separate environments exist because cygrid and OpenCV conflict:

| Task | Environment | Command |
|------|-------------|---------|
| Gridding / mapping | `~/.venv/cygrid` | `source activate_cygrid.sh` |
| PCA decompose / correct | `~/.venv/oi` | `source activate_oi.sh` |

Never assume the active environment. Check before running commands that
depend on cygrid or OpenCV.

---

## Gridding rules

- **Always use cygrid** for spatial gridding. scipy is a last-resort fallback
  only — do not write new code that depends on it.
- For post-gridding smoothing (e.g. in `collapse_cube`), use
  `astropy.convolution` (Gaussian2DKernel + convolve with
  `nan_treatment='interpolate'`). Not scipy.ndimage.
- Kernel defaults differ from CLASS: kernel FWHM = beam, pixel = beam/3
  (CLASS uses FWHM = beam/3, pixel = beam/2). See `gridding_vs_class.md`.
- Default display colormap is `rainbow`, full range `(0, 100)` percentile
  clip — matching CLASS `gt_lut_default`. Old default was `inferno` at
  `(2, 98)`; revert with `--colormap inferno --percentile-clip 2 98`.

---

## Display / colormap notes

- For `--mode peak-intensity` maps, edge pixels have high noise and dominate
  the colormap, making the core invisible. Suggest `--percentile-clip 0 98`
  or `--suppress-high <value>` when this happens.

---

## FITS data model — critical facts

### Velocity axis
```
v(i) = VELOCITY[0] + (i - (CRPIX1 - 1)) * DELTAV
```
- All units are **m/s** internally; display in km/s.
- `CRPIX1` lives in the **table header** (`hdul[1].header`), not the primary
  header and not in the data columns. This is a common source of bugs.
- `VELOCITY` and `DELTAV` are per-row columns but are constant across rows.

### Coordinates
- `RA  = CRVAL2 + CDELT2` (degrees)
- `Dec = CRVAL3 + CDELT3` (degrees)
- RA axis is inverted in display (East = left, standard sky orientation).

### Calibration rows (OBJECT column values)
These must be preserved through filtering — never apply NaN/RMS filters to them:

| OBJECT | Role |
|--------|------|
| `TSYS` | System temperature — needed for noise weighting |
| `TAU_SIG` | Atmospheric opacity — must stay in [0.001, 1.0] per channel |
| `SKYCHOPDIFF` | Sky chop-diff — PCA decomposition source |
| `SKYDIFF`, `SKY-DIFF` | Equivalent sky references |
| `TREC (SSB)` | Receiver noise only — can be excluded from science processing |

### Index columns
`TSYS_INDEX` and `TAU_SIG_INDEX` (added by `prepare_for_pca`) are row-index
pointers. They must be **remapped** whenever rows are deleted. A value of `-1`
means no valid calibration pairing. This remapping already happens in
`filter_and_save_fits` and `split_fits_by_mission`.

---

## Pipeline stage order

```
combine_fits → split_fits → filter_fits → apply_baseline → prepare_for_pca
→ pca_decompose → pca_correct → post_process_data → map/datacube
```

Config output paths for each stage:
- `[output].clean_fits` — after filter_fits
- `[output].reduced_fits` — after apply_baseline
- `[output].prepared_for_pca` — after prepare_for_pca
- `[output].pcad_fits` — after pca_correct
- `[output].post_processed_fits` — after post_process_data
- `[output].datacube` — spectral cube

---

## PCA correction — known open issues (todo.txt)

Two algorithmic issues vs the legacy pyclass behaviour are documented in
`todo.txt` and not yet fixed:

1. **Re-fit after component rejection** — our code re-fits coefficients using
   only surviving components; pyclass uses the original full-set coefficients.
   This inflates subtraction.

2. **Mean centering missing** — sklearn PCA subtracts `pca_model.mean_` at
   decompose time, so the correct projection is
   `coeff = components · (spectrum - mean_sky)`. We currently omit the
   `mean_sky` term, causing a fixed bias on every coefficient.

Do not silently fix these without discussing the approach first — both
touch the core subtraction path and need careful verification.

### Completed PCA fixes (2026-05-11)
- Double-smoothing of components removed
- Spectrum smoothing order fixed (smooth first, then mask)
- Per-component noise ratio cutoffs wired in
- `cut_coefficients` normalisation fixed (raw K units, not sigma units)

---

## Quality filtering

- `RMSRATIO` / `RMSRATIOB` ≈ 1.0 is ideal. Gridding uses Gaussian weight
  `exp(-(value - 1.0)² / (2 × 0.5²))`.
- `filter_fits --filter-below RMSRATIOB 2.0` keeps spectra with value ≤ 2.

---

## Code conventions

- Process spectra in chunks of 2000 to limit peak RAM (already established
  pattern in `filter_and_save_fits` and `reduce_spectra`).
- Do not load entire FITS files into RAM — use `memmap=True` where possible.
- Calibration rows (TSYS, TAU_SIG) must never be modified by baseline
  subtraction, smoothing, or any reduction step.
- When removing rows, always remap TSYS_INDEX and TAU_SIG_INDEX.
