#!/usr/bin/env python3
"""
pca_sweep_complete.py — Grid search over pca_correct parameters for the M51 complete dataset.

For each parameter combination the script runs:
  1. pca_correct        (~/.venv/oi)      — apply PCA correction
  2. post_process_data  (~/.venv/oi)      — refill telluric noise, compute RMS
  3. filter_fits        (~/.venv/oi)      — apply RMSRATIOB, peak, and tau filters
  4. create_datacube    (~/.venv/cygrid)  — grid post_1.5.fits into a temporary cube
  5. collapse_cube      (~/.venv/cygrid)  — moment-0 map saved with parameter tag

Plot files are saved to PLOT_DIR with filenames encoding all parameter values.
FITS intermediate files are overwritten on every iteration.

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
CONFIG = str(Path("~/.venv/oi").expanduser().parent.parent /
             "software/oi_zeigt/config_complete.toml")
# Use absolute path so the script can be run from any directory
CONFIG = str(Path("/home/verbena/software/oi_zeigt/config_complete.toml"))

# Virtual environment executables
OI_BIN      = Path("~/.venv/oi/bin").expanduser()
CYGRID_BIN  = Path("~/.venv/cygrid/bin").expanduser()

# Fixed paths (must match config_complete.toml [output] section)
POST_FITS      = "/home/diskB/data_soft/m51/processed_complete_fits/post_processed.fits"
POST_15_FITS   = "/home/diskB/data_soft/m51/processed_complete_fits/post_1.5.fits"
DATACUBE_TEMP  = "/home/diskB/data_soft/m51/processed_complete_fits/datacube_sweep_temp.fits"

# Where collapse_cube plots are saved
PLOT_DIR = Path("/home/diskB/data_soft/m51/pca_sweep_complete")

# collapse_cube parameters (from procedure_complete.txt second block)
COLLAPSE_VEL_MIN  = 400
COLLAPSE_VEL_MAX  = 570
COLLAPSE_ZOOM     = 1
COLLAPSE_ZOOM_RA  = "13:29:53"
COLLAPSE_ZOOM_DEC = "47:11:30"

# Fixed pca_correct parameters (not swept)
LINE_CUTOFF_STD = 1.0

# ---------------------------------------------------------------------------
# GRID — edit values here
# ---------------------------------------------------------------------------
# variance_cutoff: component used only if its variance ratio > cutoff.
#   Lower → more components (more aggressive).  Config default: 0.005.
VARIANCE_CUTOFFS = [0.005, 0.01]

# noise_ratio_cutoff: skip a component if spectrum_std / component_std > this.
#   Lower → fewer components applied.  Procedure range: 1.0 – 5.0.
NOISE_RATIO_CUTOFFS = [1.0, 3.0, 5.0, 10.0]

# line_kernel_size: Gaussian blur kernel for OpenCV science-line detection.
#   Must be odd.  Procedure uses 31.
LINE_KERNEL_SIZES = [11, 21, 31]

# smoothing_kernel: boxcar applied to PCA components before subtraction.
#   Procedure uses 5.
SMOOTHING_KERNELS = [5, 11, 21]

# cut_coefficients: skip component if |coeff| < this value.  0 = disabled.
CUT_COEFFICIENTS = [0, 0.05, 0.1]

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def run(cmd: list[str], label: str) -> bool:
    """Run a command, stream its output, return True on success."""
    print(f"\n{'='*70}")
    print(f"  {label}")
    print(f"  {' '.join(cmd)}")
    print(f"{'='*70}")
    result = subprocess.run(cmd)
    if result.returncode != 0:
        print(f"  *** FAILED (exit {result.returncode}) — skipping rest of combination ***")
        return False
    return True


def param_tag(vc, nr, lk, sk, cc) -> str:
    sk_str = str(sk) if sk is not None else "none"
    return f"vc{vc}_nr{nr}_lk{lk}_sk{sk_str}_cc{cc}"


AOR_ID_FILTER = "0507,0508,0408,0407,506,406,606,607,608"


def build_pca_correct_cmd(vc, nr, lk, sk, cc) -> list[str]:
    cmd = [
        str(OI_BIN / "pca_correct"),
        "--config", CONFIG,
        "--variance-cutoff", str(vc),
        "--noise-ratio-cutoff", str(nr),
        "--line-kernel-size", str(lk),
        "--line-cutoff-std", str(LINE_CUTOFF_STD),
        "--aor-id", AOR_ID_FILTER,
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
        "--zoom-ra", COLLAPSE_ZOOM_RA,
        "--zoom-dec", COLLAPSE_ZOOM_DEC,
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
    ))

    n_total = len(grid)
    print(f"\nPCA parameter sweep (complete dataset): {n_total} combinations")
    print(f"Config:  {CONFIG}")
    print(f"Plots → {PLOT_DIR}\n")

    results = []   # (tag, status, elapsed)

    for i, (vc, nr, lk, sk, cc) in enumerate(grid, 1):
        tag = param_tag(vc, nr, lk, sk, cc)
        plot_path = PLOT_DIR / f"collapse_{tag}.png"

        print(f"\n{'#'*70}")
        print(f"  Combination {i}/{n_total}: {tag}")
        print(f"  variance_cutoff={vc}  noise_ratio_cutoff={nr}  line_kernel={lk}  "
              f"smoothing_kernel={sk}  cut_coefficients={cc}")
        print(f"{'#'*70}")

        t0 = time.time()
        failed_at = None

        if not run(build_pca_correct_cmd(vc, nr, lk, sk, cc),
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
        print(f"  {mark}  {tag:<50}  ({elapsed:.0f}s)  {status}")

    print(f"\nPlots saved to: {PLOT_DIR}")
    return 0 if n_fail == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
