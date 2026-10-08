# oi_zeigt

A spectral reduction and mapping pipeline for the SOFIA **GREAT** archive.
It takes calibrated GREAT spectra stored as FITS binary tables and turns them
into cleaned, baseline-subtracted, PCA-corrected spectra and gridded 3D
datacubes and 2D maps. It was built for [OI] 63 µm but works for any GREAT
line (e.g. [CII] 158 µm).

The pipeline is a set of standalone command-line tools. Each one reads a
shared `config.toml` and can also be run on its own for inspection and
diagnostics.

> FITS data files are not part of this repository and are distributed
> separately.

## Features

- **Filtering**: reject spectra by object, blank (NaN) fraction, all-zero
  content, column values (τ, Tsys, RMS, …), and per-mission drop rules.
- **Reduction**: unblanking, velocity extraction, telluric-window noise
  filling, polynomial baselines, smoothing and decimation.
- **PCA correction**: builds principal components from reference spectra
  (by default `SKYCHOPDIFF`) per mission and receiver pixel, then subtracts
  the fitted systematics (standing waves, baseline ripples) from the
  science spectra. Science lines are detected and protected during fitting.
- **Quality metrics**: RMS, RMS ratio, whiteness (Allan variance) and
  histograms of any metric.
- **Gridding**: 3D datacubes with [cygrid](https://github.com/bwinkel/cygrid)
  Gaussian-kernel gridding, per-spectrum and per-channel weighting, and
  `COVERAGE` / `WEIGHT_MAP` / `WEIGHT_CUBE` extensions. The output can be
  read by CLASS.
- **Maps**: moment-0, peak and peak-range maps, velocity-field stacking
  (`--shuffle`), side-by-side comparisons and colour + contour overlays.
- **Scaling to large data**: split archives per mission/pixel/scan and run
  the whole chain in parallel with Snakemake.

## Installation

Two virtual environments are needed because of a NumPy conflict: OpenCV
(used for line detection) requires NumPy 2.x, while the cygrid wheels are
built for NumPy 1.x.

```bash
# Main environment: filtering, reduction, PCA
python3 -m venv ~/.venv/oi
source ~/.venv/oi/bin/activate
pip install -e ".[dev]" opencv-python

# Gridding environment: create_datacube and mapping
python3 -m venv ~/.venv/cygrid
source ~/.venv/cygrid/bin/activate
pip install "numpy<2" cygrid
pip install -e .
```

`activate_oi.sh` and `activate_cygrid.sh` are shortcuts for activating each
environment. See [docs/venv_setup.md](docs/venv_setup.md) for details.

## Quick start

1. Copy `config.toml` and set `[input].fits_file` and the `[output]` paths.
2. Run the pipeline step by step:

```bash
source activate_oi.sh
filter_fits      --config config.toml            # → [output].clean_fits / rejected_fits
reduce_spectra   --config config.toml            # → [output].reduced_fits
prepare_for_pca  --config config.toml            # → [output].prepared_for_pca
pca_decompose    --config config.toml --plot     # components per mission/pixel
pca_correct      --config config.toml            # → [output].pcad_fits
post_process_data --config config.toml           # optional: baseline, RMS, whiteness

source activate_cygrid.sh
create_datacube  --config config.toml            # → [output].datacube
collapse_cube    --fits datacube.fits --velocity-range 450 550
```

To process a full archive in parallel:

```bash
snakemake --cores 8 -np    # dry run
snakemake --cores 8        # split → filter → reduce → prepare → decompose → correct → merge → datacube
```

Snakemake calls each step with the right virtual environment automatically.
To use a different config file, pass `--config config_file=other.toml`.

Every command has `--help`. Most commands accept `--clean`, `--reduced`,
`--prepared`, `--pcad` or `--rejected` to read the matching `[output]` file
from the config instead of an explicit `--fits` path.

## Commands

### Inspection
| Command | Purpose |
|---|---|
| `print_fits_info` | Summary of a FITS file: columns, objects, row counts |
| `print_pca_parameters` | Unique `PCAPARAM` values in a file |
| `analyze_blanks` | Identify the value used as blank/missing marker |
| `plot_sample_spectra` | Plot a sample of spectra, optionally filtered by object |
| `plot_sample_raw` | Same as `plot_sample_spectra`, but plots the `RAW` column |
| `plot_spectra` | Averaged spectra per receiver, polarization or pixel; `--group-by` scan/telescope diagnostics |
| `plot_raw` | Same as `plot_spectra`, but averages the `RAW` column |
| `plot_shobs_window` | Scatter of S-H_OBS intensity averaged over a small window; one panel per receiver (TELESCOP), one colour per input FITS |
| `plot_skyfit_pairs` | S-H_OBS overlaid with its paired S-H_FIT (+ residual), one page per scan.subscan, multi-page PDF per input |
| `chi_sqr_skyfit` | Compare sky-fit CHI_SQR (one per scan.subscan) between files: stats, cycle-by-cycle wins, verdict; also printed by `plot_skyfit_pairs` |
| `plot_skies` | Plot sky/reference spectra (`SKYCHOPDIFF`, `SKY-DIFF`) |
| `plot_skyobsfit` | Compare observed vs fitted sky spectra |
| `examine_telluric` | Average telluric spectrum per mission/pixel |
| `cascade_plots` | Waterfall images of spectra per mission/scan/pixel |

### Filtering and reduction
| Command | Purpose |
|---|---|
| `filter_fits` | Split into clean/rejected files by object, NaN fraction, zeros, column cuts |
| `filter_missions` | Drop rows according to `mission_id_parameters.yml` rules |
| `reduce_spectra` | Unblank, extract velocity range, fill telluric window, baseline, smooth |
| `apply_baseline` | Standalone polynomial baseline subtraction |
| `average` | Average spectra (grouped by object by default) |

### PCA correction
| Command | Purpose |
|---|---|
| `prepare_for_pca` | Select PCA source + object, fill telluric regions, build Tsys/τ indices |
| `pca_decompose` | PCA of reference spectra per `MISSION_ID` / `TELESCOP` |
| `pca_correct` | Fit and subtract components from science spectra |
| `determine_pca_parameters` | Recommend `n_components` from spatial pixel coherence |
| `post_process_data` | Post-correction baseline, smoothing, RMS, whiteness check |

### Quality metrics
| Command | Purpose |
|---|---|
| `rmsratio` | RMS-ratio analysis and histogram |
| `spechistogram` | Multi-panel histograms of quality metrics |
| `map_column` | Spatial map of any column (e.g. `TSYS`, `TAU-ATM`, `RMSRATIO`) |

### Gridding and maps
| Command | Purpose |
|---|---|
| `create_datacube` | Grid spectra into a 3D (v, dec, ra) cube with WCS and weights |
| `collapse_cube` | Moment-0 / peak maps from a cube; `--shuffle`, masks, clipping |
| `compare_maps` | Collapse several cubes and show them side by side |
| `overlay_maps` | One cube as colour, another as contours |
| `map_integrated` | Integrated-intensity map directly from spectra |
| `compare_map_integrated` | Side-by-side integrated maps of pipeline stages |

### File handling
| Command | Purpose |
|---|---|
| `split_fits` | Split a large file per (`MISSION_ID`, `TELESCOP`, `SCAN`) |
| `combine_fits` | Merge FITS files (from a list file or arguments) |

## Configuration

Everything is driven by a TOML file (see `config.toml`):

| Section | Contents |
|---|---|
| `[input]` | Input FITS file |
| `[parameters]` | Object name filter (e.g. `M51CENTER`) |
| `[filters]` | Blank fraction, AOR IDs to remove |
| `[reduction]` | Baseline order, line window, extraction range, smoothing |
| `[pca_prepare]` | Telluric noise filling, τ limits |
| `[pca]` | Components, source, smoothing, noise/variance cutoffs, line detection |
| `[gridding]` | Method, beam size, pixel size, spectrum/channel weighting |
| `[output]` | Paths of all intermediate and final products |

Per-mission PCA overrides and drop rules are in
`src/oi_zeigt/pca_analysis/mission_id_parameters.yml`.

## Package layout

```
src/oi_zeigt/
├── cli.py              # click entry points for most commands
├── basic_io.py         # FITS reading/writing helpers
├── cascade_plots.py
├── reduction/          # blank detection, filtering, splitting, baselines
├── pca_analysis/       # prepare, decompose, correct, line detection, parameter search
├── mapping/            # cygrid gridding, datacubes, map collapse, CLASS colour tables
└── statistics/         # spectral quality metrics
docs/                   # reference documentation
tests/                  # tests and inspection scripts
Snakefile               # parallel pipeline workflow
```

## Documentation

More detail is in [`docs/`](docs/README.md):

- [PCA parameters](docs/pca_parameters.md) and [differences from the legacy pyclass PCA](docs/pca_differences.md)
- [Line detection](docs/line_detection.md) and [Tsys/τ pairing](docs/tsys_pairing.md)
- [Datacube weighting](docs/create_datacube_weighting.md), [gridding scales](docs/gridding_scales.md) and [comparison with CLASS](docs/gridding_vs_class.md)
- [WCS conventions](docs/wcs_coordinates.md)

## Tests

```bash
pytest tests/
```

Many tests need local FITS data, which is not included in the repository.
