# oi_zeigt — Codebase Reference

A Python package for processing and analyzing SOFIA GREAT observations of the
[OI] 63 µm fine-structure line. The name means roughly "OI shows" in German.

---

## What this package does

Raw SOFIA GREAT data arrives as FITS binary tables — one row per spectrum per
pointing, with many calibration rows mixed in. The package turns that into:

1. Clean science spectra (filtered, baseline-subtracted, quality-checked)
2. A PCA-corrected version (telluric/instrumental artefacts removed)
3. Spatial maps and 3D datacubes ready for analysis or CLASS/CASA

---

## Data model — the FITS binary table

Every FITS file handled here is a **binary table** (BinTableHDU) where each
row is one spectrum. The SPECTRUM column is a 1D array of ~1264 float values.

### Key columns

| Column | Type | Meaning |
|--------|------|---------|
| `OBJECT` | str | Sky source name or row type (see below) |
| `SPECTRUM` | float32[N] | Intensity values in K (antenna temperature) |
| `VELOCITY` | float | Reference velocity in **m/s** at CRPIX1 |
| `DELTAV` | float | Channel spacing in **m/s** (can be negative) |
| `CDELT2`, `CRVAL2` | float | RA offset and reference (deg); `RA = CRVAL2 + CDELT2` |
| `CDELT3`, `CRVAL3` | float | Dec offset and reference (deg); `Dec = CRVAL3 + CDELT3` |
| `MISSION_ID` | str | Flight identifier (e.g. `F528`) |
| `TELESCOP` | str | Receiver/pixel ID (e.g. `LFAV_PX00_S`) |
| `AOR_ID` | str | Astronomical Observation Request identifier |
| `SCAN` | int | Scan number within a flight |
| `OBSMODE` | str | Observing mode (e.g. `BSAB`, `OTFSWB`) |
| `RMSRATIO` | float | Off-source RMS / expected noise — quality metric |
| `RMSRATIOB` | float | Baseline-region version of RMSRATIO |
| `TSYS_INDEX` | int32 | Row index of the paired TSYS calibration spectrum |
| `TAU_SIG_INDEX` | int32 | Row index of the paired TAU_SIG opacity spectrum |

`TSYS_INDEX` and `TAU_SIG_INDEX` are added by `prepare_for_pca` and are
remapped correctly whenever rows are removed (e.g. by `filter_fits`,
`split_fits`).

### Row types (OBJECT column)

Science rows have the target name (e.g. `M51CENTER`). Calibration rows are
mixed in and must usually be kept:

| OBJECT value | Meaning |
|---|---|
| `TSYS` | System temperature spectrum — used for noise weighting |
| `TAU_SIG` | Atmospheric opacity spectrum — values should be in ~[0.001, 1.0] |
| `SKYCHOPDIFF` | Sky chop-difference — PCA source spectrum |
| `SKYDIFF`, `SKY-DIFF` | Equivalent sky references |
| `TREC (SSB)` | Receiver noise only — no sky signal, can be skipped |

### Velocity axis formula

```
v(i) = VELOCITY[0] + (i - (CRPIX1 - 1)) * DELTAV
```

`CRPIX1` is stored in the **table header** (1-indexed, FITS convention).
`VELOCITY` and `DELTAV` are per-row columns but are constant across all rows.
The velocity axis is always in m/s internally; the CLI plots in km/s.

---

## Pipeline overview

```
raw FITS
   │
   ├── combine_fits       combine multiple mission files into one
   │
   ├── split_fits         split by (MISSION_ID, TELESCOP, SCAN) → chunks
   │
   ├── filter_fits        NaN/zero/column filtering → clean + rejected
   │
   ├── apply_baseline     polynomial baseline subtraction
   │  (= reduce_spectra with baseline step)
   │
   ├── prepare_for_pca    pair TSYS/TAU_SIG rows, fill NaN with noise,
   │                      add TSYS_INDEX and TAU_SIG_INDEX columns
   │
   ├── pca_decompose      SVD on sky (SKYCHOPDIFF) spectra per mission
   │
   ├── pca_correct        project science spectra onto PCA components,
   │                      subtract the fitted sky contribution
   │
   ├── post_process_data  any final cleaning step after PCA
   │
   ├── map_integrated     integrated-intensity 2D map (via cygrid)
   ├── create_datacube    full 3D spectral cube (nvel × ny × nx)
   └── collapse_cube      collapse cube to 2D (peak, moment, etc.)
```

Each step reads from a FITS file and writes to a new FITS file; the paths are
controlled by `config.toml`.

---

## config.toml structure

The config file (`config.toml`) drives almost every step. All CLI tools accept
`--config` to point to it; many search the current and parent directories
automatically.

```toml
[input]
fits_file = "..."          # raw/combined input

[parameters]
object = "M51CENTER"       # science target substring for filtering

[filters]
blank_fraction = 0.20      # max NaN fraction before rejecting
remove_aor_id = "..."      # optional AOR_ID to always reject

[reduction]
baseline = 3               # polynomial order
line_window = [450, 550]   # channels to exclude from baseline fit
extract = [350, 700]       # channel range to keep after reduction
smooth = 5                 # boxcar smoothing kernel
decimate = true

[gridding]
beamsize_arcsec = 14.1
method = "cygrid"          # or "scipy" fallback
spectra_weighting = "rmsratio"
channel_weighting = true
telescop = "IRAM-30M"

[pca_prepare]
fill_with_noise = true
tau_min = 0.001
tau_max = 1.0

[pca]
n_components = 5
pca_source = "SKYCHOPDIFF"
smoothing_kernel_size = 5
variance_cutoff = 0.005
noise_ratio_cutoff = 1.0
line_kernel_size = 21
line_cutoff_std = 1.5
least_squares = true

[output]
clean_fits = "..."
rejected_fits = "..."
reduced_fits = "..."
prepared_for_pca = "..."
pcad_fits = "..."
datacube = "..."
pca_plots_dir = "..."
```

---

## Source layout

```
src/oi_zeigt/
├── __init__.py
├── basic_io.py            I/O: read/write FITS, load config, combine files
├── cli.py                 All CLI entry points (Click commands)
├── reduction/
│   └── core.py            Filtering, baseline, smoothing, extraction,
│                          split_fits_by_mission, filter_and_save_fits
├── mapping/
│   └── gridding.py        WCS map creation, cygrid gridding, datacube,
│                          collapse, integrated-intensity map
├── pca_analysis/
│   ├── config.py          Read PCA parameters from config/YAML
│   ├── core.py            sklearn PCA wrapper on spectral datacubes
│   ├── decompose.py       SVD decomposition per mission
│   ├── correct.py         Apply PCA correction to science spectra
│   ├── fits_indexing.py   Pair TSYS/TAU_SIG rows to science rows
│   ├── line_detection.py  Identify spectral line channels
│   ├── utilities.py       Shared PCA helpers
│   ├── pca_decompose_per_mission.py  CLI entry for pca_decompose
│   ├── pca_correct_fits.py           CLI entry for pca_correct
│   └── prepare_for_pca.py            TSYS/TAU pairing + noise fill
└── statistics/
    └── quality.py         RMSRATIO histograms, quality metrics
```

---

## Module summaries

### `basic_io.py`

- `get_config(path)` — load TOML config (searches up the directory tree)
- `read_fits(path)` / `read_fits_from_config()` — open FITS via astropy
- `combine_fits_files(list, single_hdu=True)` — stream-merge multiple FITS
  into one, disk-backed to avoid RAM blowup
- `combine_fits_from_list(list_txt, output)` — read paths from a text file
- `reconstruct_velocity_axis(hdu)` — build the velocity axis numpy array

### `reduction/core.py`

- `detect_nan_channels(spectrum)` → `(mask, fraction)` — NaN check
- `detect_blank_channels(spectrum)` → `(mask, fraction)` — zero/blank check
- `filter_and_save_fits(hdul, ...)` — the main filtering engine:
  - Separates science vs calibration rows
  - Applies NaN, zero, peak, tau, column-value, flight, OBSMODE filters
  - Remaps TSYS_INDEX / TAU_SIG_INDEX after rows are removed
  - Returns `(clean_path, rejected_path, stats_dict)`
- `split_fits_by_mission(fits_path, output_dir)` — splits into per-scan chunks,
  remaps index columns, writes science-only manifest
- `baseline_subtract(spectrum, order, window)` — polynomial baseline with
  iterative sigma-clipping; does not touch calibration rows
- `reduce_spectra(hdul, methods)` — orchestrates unblank → extract → baseline
  → smooth in sequence

### `mapping/gridding.py`

- `grid_to_map(ras, decs, values, ...)` — the core gridding primitive;
  uses **cygrid** (preferred) or scipy fallback; handles weighted averaging
- `create_integrated_map(hdul, ...)` — integrate spectra, grid to 2D map,
  supports per-spectrum RMSRATIO weights and per-channel tau/Tsys weights
- `create_spectral_datacube(hdul, ...)` — grid every channel → 3D cube
  `(nvel, ny, nx)`; parallelised over channels with multiprocessing
- `save_map_to_fits(...)` — write map with proper WCS header, beam keywords,
  optional spectral table and coverage/weight extensions
- Coordinate convention: `RA = CRVAL2 + CDELT2`, `Dec = CRVAL3 + CDELT3`.
  RA axis is inverted in display (increases right → left).

### `pca_analysis/`

The PCA pipeline works in two steps:

1. **`pca_decompose`** — for each mission, extract SKYCHOPDIFF spectra, run
   SVD, keep components that explain significant variance (`variance_cutoff`)
   but don't lie along the emission line (`line_detection`). Results saved as
   a `.pkl` file per mission.

2. **`pca_correct`** — for each science spectrum, project onto saved PCA
   components using least-squares, subtract the reconstructed sky contribution.
   Line channels are protected from contaminating the fit.

`prepare_for_pca` runs before both: it pairs each science row with its TSYS
and TAU_SIG rows (adding `TSYS_INDEX`/`TAU_SIG_INDEX`), optionally fills NaN
channels with noise, and clips bad tau values.

---

## CLI commands (entry points)

| Command | Does |
|---------|------|
| `print_oifits_info` | Print HDU/column/object/AOR summary |
| `plot_sample_spectra` | Plot N random spectra for a target object |
| `plot_skies` | Plot sky (SKYCHOPDIFF) spectra |
| `plot_skyobsfit` | Compare S-H_OBS vs S-H_SKY pairs |
| `analyze_blanks` | Stats on blank/zero/NaN values |
| `filter_fits` | Main NaN/quality filter; produces clean + rejected |
| `filter_missions` | Keep/drop specific missions |
| `split_fits` | Split by mission/scan |
| `combine_fits` | Merge a list of FITS files |
| `apply_baseline` | Polynomial baseline subtraction |
| `reduce_spectra` | Multi-step reduction (unblank/extract/baseline/smooth) |
| `average` | Average spectra |
| `prepare_for_pca` | Add TSYS/TAU pairing + noise fill |
| `pca_decompose` | SVD on sky spectra per mission |
| `pca_correct` | Apply PCA sky correction to science spectra |
| `post_process_data` | Post-PCA filtering |
| `spechistogram` | Histogram of a FITS column |
| `rmsratio` | RMSRATIO distribution plot |
| `map_column` | 2D map from a scalar column |
| `map_integrated` | Integrated-intensity 2D map |
| `compare_map_integrated` | Side-by-side comparison of two maps |
| `create_datacube` | Full 3D spectral cube |
| `collapse_cube` | Collapse cube to 2D |
| `examine_telluric` | Inspect telluric/atmospheric features |
| `print_pca_parameters` | Show PCA parameters from config |

Most commands accept `--config config.toml` plus stage flags like `--clean`,
`--reduced`, `--prepared`, `--pcad`, `--post` to select the input FITS from
the corresponding `[output]` path.

---

## Quality columns

| Column | Meaning | Ideal value |
|--------|---------|-------------|
| `RMSRATIO` | Measured off-source RMS / theoretical noise | ≈ 1.0 |
| `RMSRATIOB` | Same but computed from baseline channels | ≈ 1.0 |
| `TSYS` spectrum | System temperature in K (per channel) | instrument-dependent |
| `TAU_SIG` spectrum | Atmospheric opacity (dimensionless, per channel) | 0.001–1.0 |

`filter_fits --filter-below RMSRATIOB 2.0` keeps spectra with `RMSRATIOB ≤ 2`.
The gridding uses a Gaussian weight `exp(-(RMSRATIO - 1)² / (2×0.5²))` so
near-ideal spectra (RMSRATIO ≈ 1) get full weight.

---

## Important implementation notes

- **Gridding tool is cygrid** — not scipy. scipy is a fallback only.
  Activate the cygrid virtualenv with `source activate_cygrid.sh` before
  running mapping commands.

- **Index remapping** — whenever rows are deleted from a FITS file
  (`filter_and_save_fits`, `split_fits_by_mission`), `TSYS_INDEX` and
  `TAU_SIG_INDEX` are remapped to reflect new row positions. A value of `-1`
  means no valid calibration pairing.

- **Calibration rows are protected** — `filter_fits` never applies NaN
  filtering, RMSRATIO filtering, or OBSMODE exclusion to `TSYS` or `TAU_SIG`
  rows. They are always passed through to the clean file.

- **Chunked processing** — `filter_and_save_fits` and `reduce_spectra` process
  spectra in chunks of 2000 to limit peak RAM regardless of file size.

- **CRPIX1 location** — lives in the **table header** (`hdul[1].header`), not
  the primary header, and not in the data columns. This is a common gotcha.

- **RA convention** — CDELT2 can be positive or negative depending on the
  observation. The display always inverts the RA axis so East is left.

- **PCA source** — default `pca_source = "SKYCHOPDIFF"`. Each mission gets its
  own PCA decomposition because the sky spectrum shape varies flight to flight.

---

## Typical workflow

```bash
# 1. Inspect raw data
print_oifits_info --config config.toml

# 2. Filter bad spectra
filter_fits --config config.toml

# 3. Baseline subtraction
apply_baseline --config config.toml

# 4. Prepare for PCA (adds TSYS/TAU indexing)
prepare_for_pca --config config.toml

# 5. PCA decomposition (per mission)
pca_decompose --config config.toml

# 6. PCA correction
pca_correct --config config.toml

# 7. Make maps
map_integrated --config config.toml --pcad
create_datacube --config config.toml --pcad
```
