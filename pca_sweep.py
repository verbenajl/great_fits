#!/usr/bin/env python3
"""
pca_sweep.py — Grid search over pca_correct parameters.

For each parameter combination the script runs:
  1. pca_correct        — apply PCA correction
  2. post_process_data  — compute RMS / refill telluric noise
  3. filter_fits        — apply RMSRATIOB, peak, and tau filters
  4. compare_map_integrated — produce a comparison plot

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
# compare_map_integrated needs cygrid → cygrid venv
OI_BIN      = Path("~/.venv/oi/bin").expanduser()
CYGRID_BIN  = Path("~/.venv/cygrid/bin").expanduser()

# Fixed paths (must match config.toml [output] section)
POST_FITS    = "/home/diskB/data_soft/m51/processed_fits/post_processed.fits"
POST_15_FITS = "/home/diskB/data_soft/m51/processed_fits/post_1.5.fits"
CLEAN_FITS   = "/home/diskB/data_soft/m51/processed_fits/clean_data.fits"
REDUCED_FITS = "/home/diskB/data_soft/m51/processed_fits/reduced_data.fits"
PREP_FITS    = "/home/diskB/data_soft/m51/processed_fits/prepared_for_pca.fits"
PCAD_FITS    = "/home/diskB/data_soft/m51/processed_fits/pca_corrected.fits"

# Where comparison plots are saved
PLOT_DIR = Path("/home/diskB/data_soft/m51/pca_sweep_plots_00")

# Velocity range for compare_map_integrated (km/s)
VEL_MIN, VEL_MAX = 500, 550

# ---------------------------------------------------------------------------
# GRID — edit values here
# ---------------------------------------------------------------------------
# variance_cutoff: per-component explained variance threshold.
#   A component is used only if its individual variance ratio > cutoff.
#   Lower → more components used → more aggressive correction.
#   Higher → fewer components → less aggressive correction.
#   Config default: 0.005
VARIANCE_CUTOFFS = [0.005, 0.008, 0.01, 0.02]

# noise_ratio_cutoff: skip a component if its noise ratio exceeds this.
#   Lower → fewer components applied.
#   Config default: 2.5
NOISE_RATIO_CUTOFFS = [1.0, 2.5, 5.0]

# line_kernel_size: Gaussian blur kernel for OpenCV science-line detection.
#   Must be odd.  Larger → smoother, misses narrow features.
#   Config default: 21
LINE_KERNEL_SIZES = [11, 21, 51]

# smoothing_kernel: boxcar size applied to PCA components before subtraction.
#   None → disabled.  Larger → smoother correction, lower spectral resolution.
#   Config default: 5
SMOOTHING_KERNELS = [None, 3, 5, 7, 11]

# cut_coefficients: skip component if |coeff| < this value.
#   0 → disabled (no coefficient cutoff).
#   Larger → fewer components applied (only those with strong projection).
#   Original pyclass default: 0.1; M51 config: disabled (0).
CUT_COEFFICIENTS = [0, 0.05, 0.1]

# global_noise_ratio_cutoff: enables KDE-based per-group component rejection.
#   None → disabled (falls back to simple noise_ratio_cutoff threshold).
#   When set, fits a KDE to the noise ratio distribution across all spectra
#   per group per component; rejects component globally if KDE peak > this value.
GLOBAL_NOISE_RATIO_CUTOFFS = [None, 3.0, 5.0]

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
    gnr_str = str(gnr) if gnr is not None else "none"
    return f"vc{vc}_nr{nr}_lk{lk}_sk{sk_str}_cc{cc}_gnr{gnr_str}"


def build_pca_correct_cmd(vc, nr, lk, sk, cc, gnr, plot_dir: Path) -> list[str]:
    cmd = [
        str(OI_BIN / "pca_correct"),
        "--config", CONFIG,
        "--variance-cutoff", str(vc),
        "--noise-ratio-cutoff", str(nr),
        "--line-kernel-size", str(lk),
        "--plot",
        "--plot-dir", str(plot_dir),
    ]
    if sk is not None:
        cmd += ["--smoothing-kernel", str(sk)]
    if cc:
        cmd += ["--cut-coefficients", str(cc)]
    if gnr is not None:
        cmd += ["--global-noise-ratio-cutoff", str(gnr)]
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


def build_compare_cmd(plot_path: Path) -> list[str]:
    return [
        str(CYGRID_BIN / "compare_map_integrated"),
        "--config", CONFIG,
        "--fits", CLEAN_FITS,
        "--fits", REDUCED_FITS,
        "--fits", PREP_FITS,
        "--fits", PCAD_FITS,
        "--fits", POST_FITS,
        "--fits", POST_15_FITS,
        "--velocity-range", str(VEL_MIN), str(VEL_MAX),
        "--weight-spectra",
        "--weight-channels",
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
        pca_plot_dir = PLOT_DIR / tag
        plot_path = PLOT_DIR / f"compare_{tag}.png"

        print(f"\n{'#'*70}")
        print(f"  Combination {i}/{n_total}: {tag}")
        print(f"  variance_cutoff={vc}  noise_ratio_cutoff={nr}  line_kernel={lk}  "
              f"smoothing_kernel={sk}  cut_coefficients={cc}  global_noise_ratio_cutoff={gnr}")
        print(f"{'#'*70}")

        t0 = time.time()
        failed_at = None

        if not run(build_pca_correct_cmd(vc, nr, lk, sk, cc, gnr, pca_plot_dir),
                   "Step 1/4 — pca_correct"):
            failed_at = "pca_correct"
        elif not run(build_post_process_cmd(),
                     "Step 2/4 — post_process_data"):
            failed_at = "post_process_data"
        elif not run(build_filter_fits_cmd(),
                     "Step 3/4 — filter_fits"):
            failed_at = "filter_fits"
        elif not run(build_compare_cmd(plot_path),
                     "Step 4/4 — compare_map_integrated"):
            failed_at = "compare_map_integrated"

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
