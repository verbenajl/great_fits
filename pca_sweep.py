#!/usr/bin/env python3
"""
pca_sweep.py — Grid search over pca_correct parameters.

For each parameter combination the script runs:
  1. pca_correct        — apply PCA correction
  2. post_process_data  — compute RMS / refill telluric noise
  3. filter_fits        — apply RMSRATIOB, peak, and tau filters
  4. create_datacube    — grid post_1.5.fits into a temporary 3D cube
  5. collapse_cube      — moment-0 map + spectra plot (saved with parameter tag)

Plot files are saved to PLOT_DIR with filenames that encode
all parameter values.  FITS files are overwritten on every run.

Edit the GRID section below to change the parameter space.
"""

import itertools
import subprocess
import sys
import time
from pathlib import Path

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
CONFIG = "config.toml"          # path to config.toml (relative or absolute)

# Virtual environment executables
# pca_correct / post_process_data / filter_fits need OpenCV → oi venv
# create_datacube / collapse_cube need cygrid → cygrid venv
OI_BIN      = Path("~/.venv/oi/bin").expanduser()
CYGRID_BIN  = Path("~/.venv/cygrid/bin").expanduser()

# Fixed paths (must match config.toml [output] section)
POST_FITS      = "/home/diskB/data_soft/m51/processed_fits/post_processed.fits"
POST_15_FITS   = "/home/diskB/data_soft/m51/processed_fits/post_1.5.fits"
CLEAN_FITS     = "/home/diskB/data_soft/m51/processed_fits/clean_data.fits"
REDUCED_FITS   = "/home/diskB/data_soft/m51/processed_fits/reduced_data.fits"
PREP_FITS      = "/home/diskB/data_soft/m51/processed_fits/prepared_for_pca.fits"
PCAD_FITS      = "/home/diskB/data_soft/m51/processed_fits/pca_corrected.fits"
DATACUBE_TEMP  = "/home/diskB/data_soft/m51/processed_fits/datacube_temp.fits"

# Where collapse_cube plots are saved
PLOT_DIR = Path("/home/diskB/data_soft/m51/pca_sweep_plots_01")

# collapse_cube parameters
COLLAPSE_VEL_MIN  = 425
COLLAPSE_VEL_MAX  = 575
COLLAPSE_ZOOM     = 1.0
COLLAPSE_REGION_X = 25
COLLAPSE_REGION_Y = 25
COLLAPSE_REGION_R = 0.2

# ---------------------------------------------------------------------------
# GRID — edit values here
# ---------------------------------------------------------------------------
# variance_cutoff: per-component explained variance threshold.
#   A component is used only if its individual variance ratio > cutoff.
#   Lower → more components used → more aggressive correction.
#   Higher → fewer components → less aggressive correction.
#   Config default: 0.005
VARIANCE_CUTOFFS = [0.005, 0.01, 0.02]

# noise_ratio_cutoff: skip a component if its noise ratio exceeds this.
#   Lower → fewer components applied.
#   Config default: 2.5
NOISE_RATIO_CUTOFFS = [1.0, 2.5, 5.0]

# line_kernel_size: Gaussian blur kernel for OpenCV science-line detection.
#   Must be odd.  Larger → smoother, misses narrow features.
#   Config default: 21
LINE_KERNEL_SIZES = [11, 21, 31]

# smoothing_kernel: boxcar size applied to PCA components before subtraction.
#   None → disabled.  Larger → smoother correction, lower spectral resolution.
#   Config default: 5
SMOOTHING_KERNELS = [5, 7]

# cut_coefficients: skip component if |coeff| < this value.
#   0 → disabled (no coefficient cutoff).
#   Larger → fewer components applied (only those with strong projection).
#   Original pyclass default: 0.1; M51 config: disabled (0).
CUT_COEFFICIENTS = [0, 0.05, 0.1]

# global_noise_ratio_cutoff: KDE-based rejection — removed from pca_correct.
# Kept as a single None entry so the grid structure is unchanged.
GLOBAL_NOISE_RATIO_CUTOFFS = [None]

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def run(cmd: list[str], label: str, env: dict | None = None) -> bool:
    """Run a command, stream its output, return True on success."""
    print(f"\n{'='*70}")
    print(f"  {label}")
    print(f"  {' '.join(cmd)}")
    print(f"{'='*70}")
    result = subprocess.run(cmd, env=env)
    if result.returncode != 0:
        print(f"  *** FAILED (exit {result.returncode}) — skipping rest of combination ***")
        return False
    return True


def param_tag(vc, nr, lk, sk, cc, gnr) -> str:
    sk_str = str(sk) if sk is not None else "none"
    return f"vc{vc}_nr{nr}_lk{lk}_sk{sk_str}_cc{cc}"


def build_pca_correct_cmd(vc, nr, lk, sk, cc, gnr) -> list[str]:
    cmd = [
        str(OI_BIN / "pca_correct"),
        "--config", CONFIG,
        "--variance-cutoff", str(vc),
        "--noise-ratio-cutoff", str(nr),
        "--line-kernel-size", str(lk),
    ]
    if sk is not None:
        cmd += ["--smoothing-kernel", str(sk)]
    if cc:
        cmd += ["--cut-coefficients", str(cc)]
    return cmd


def build_post_process_cmd() -> list[str]:
    return [
        str(OI_BIN / "post_process_data"),
        "--config", CONFIG,
        "--pcad",
        "--refill-telluric-noise",
    ]


def build_filter_fits_cmd() -> list[str]:
    return [
        str(OI_BIN / "filter_fits"),
        "--config", CONFIG,
        "--fits", POST_FITS,
        "--output-clean", POST_15_FITS,
        "--filter-below", "RMSRATIOB", "1.5",
        "--filter-spectrum-peaks", "500",
        "--filter-tau",
    ]


def build_create_datacube_cmd() -> list[str]:
    return [
        str(CYGRID_BIN / "create_datacube"),
        "--config", CONFIG,
        "--fits", POST_15_FITS,
        "--output", DATACUBE_TEMP,
        "--weight-spectra",
        "--weight-channels",
    ]


def build_collapse_cube_cmd(plot_path: Path) -> list[str]:
    return [
        str(CYGRID_BIN / "collapse_cube"),
        "--fits", DATACUBE_TEMP,
        "--velocity-range", str(COLLAPSE_VEL_MIN), str(COLLAPSE_VEL_MAX),
        "--zoom", str(COLLAPSE_ZOOM),
        "--region-x", str(COLLAPSE_REGION_X),
        "--region-y", str(COLLAPSE_REGION_Y),
        "--region-radius", str(COLLAPSE_REGION_R),
        "--use-wcs",
        "--no-show",
        "--plot", str(plot_path),
    ]


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    PLOT_DIR.mkdir(parents=True, exist_ok=True)

    grid = list(itertools.product(
        VARIANCE_CUTOFFS,
        NOISE_RATIO_CUTOFFS,
        LINE_KERNEL_SIZES,
        SMOOTHING_KERNELS,
        CUT_COEFFICIENTS,
        GLOBAL_NOISE_RATIO_CUTOFFS,
    ))

    n_total = len(grid)
    print(f"\nPCA parameter sweep: {n_total} combinations")
    print(f"Plots → {PLOT_DIR}\n")

    results = []   # (tag, step_failed_at_or_None)

    for i, (vc, nr, lk, sk, cc, gnr) in enumerate(grid, 1):
        tag = param_tag(vc, nr, lk, sk, cc, gnr)
        plot_path = PLOT_DIR / f"compare_{tag}.png"

        print(f"\n{'#'*70}")
        print(f"  Combination {i}/{n_total}: {tag}")
        print(f"  variance_cutoff={vc}  noise_ratio_cutoff={nr}  line_kernel={lk}  "
              f"smoothing_kernel={sk}  cut_coefficients={cc}")
        print(f"{'#'*70}")

        t0 = time.time()
        failed_at = None

        if not run(build_pca_correct_cmd(vc, nr, lk, sk, cc, gnr),
                   "Step 1/5 — pca_correct"):
            failed_at = "pca_correct"
        elif not run(build_post_process_cmd(),
                     "Step 2/5 — post_process_data"):
            failed_at = "post_process_data"
        elif not run(build_filter_fits_cmd(),
                     "Step 3/5 — filter_fits"):
            failed_at = "filter_fits"
        elif not run(build_create_datacube_cmd(),
                     "Step 4/5 — create_datacube"):
            failed_at = "create_datacube"
        elif not run(build_collapse_cube_cmd(plot_path),
                     "Step 5/5 — collapse_cube"):
            failed_at = "collapse_cube"

        elapsed = time.time() - t0
        status = f"FAILED at {failed_at}" if failed_at else f"OK → {plot_path.name}"
        results.append((tag, status, elapsed))
        print(f"\n  [{i}/{n_total}] {tag}: {status}  ({elapsed:.0f}s)")

    # Summary
    print(f"\n\n{'='*70}")
    print(f"  SWEEP COMPLETE — {n_total} combinations")
    print(f"{'='*70}")
    n_ok = sum(1 for _, s, _ in results if s.startswith("OK"))
    n_fail = n_total - n_ok
    print(f"  Succeeded: {n_ok}   Failed: {n_fail}")
    print()
    for tag, status, elapsed in results:
        mark = "✓" if status.startswith("OK") else "✗"
        print(f"  {mark}  {tag:<45}  {status}")

    print(f"\nPlots saved to: {PLOT_DIR}")
    return 0 if n_fail == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
