#!/usr/bin/env python3
"""determine_pca_parameters — recommend per-flight PCA n_components from reduced data.

Given a directory of reduced_data_*.fits, this ALWAYS re-runs the per-pixel PCA
decomposition (via `pca_decompose`) for every flight, then chooses each flight's
n_components by PIXEL-TO-PIXEL COHERENCE of the resulting components.

Method
------
The decomposition runs INDEPENDENTLY per detector pixel (GREX/LFA: up to 14
pixels, LFAH/LFAV x PX00..PX06), so each flight yields one component pickle per
pixel. A PCA component that reflects a real, array-wide common mode (atmosphere /
receiver drift — what the correction should subtract) has the SAME spectral shape
in every pixel; a component that is just a pixel's own noise does not.

So for component index i we take the i-th component from every pixel, align them
on their shared velocity grid (velocity_axis_kms — pixels can be edge-trimmed to
slightly different channel counts), correlate every pair, and take the MEDIAN
|Pearson r|. High == coherent (keep); low == noise (drop). |r| handles PCA's
arbitrary per-pixel sign/ordering.

n_components = number of LEADING components whose coherence stays >= the
threshold (default 0.40); once one drops below, the trailing ones are noise.

Output
------
Prints a per-flight coherence table and a ready-to-paste yaml block (n_components
from coherence + the standard "correct-time" gate block). With --write, merges
those values into the mission_id_pca_parameters.yml in place, preserving the
file's comments and any flights not covered by this run.

Examples
--------
  determine_pca_parameters --dir /home/diskB/data_soft/m51/clean_all_missions
  determine_pca_parameters --dir /path/to/reduced --config config.toml --write
"""
import argparse
import glob
import os
import pickle
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

# ── Standard "correct-time" gate block (shared by every tuned flight) ──────────
# The coherence method only sets n_components; these round out each yaml block.
DEFAULT_GATES = dict(
    global_noise_ratio_cutoff=40,
    noise_ratio_cutoff=20,
    cutoff_variance=0.002,
    cut_coefficients=0.0001,
    line_kernel_size=31,
    line_cutoff_std=4,
    least_squares=False,
)
GATE_ORDER = list(DEFAULT_GATES.keys())

PKL_RE = re.compile(r"decomposition_(?P<mid>.+)_(?P<tel>LFA[HV]_PX\d\d_S)_components\.pkl$")
FLIGHT_KEY_RE = re.compile(r"^\d{4}-\d{2}-\d{2}_\S*F\d+:\s*$")


# ── Coherence scoring ──────────────────────────────────────────────────────────
def load_pixel_components(components_dir, mission_id):
    """Return (comps, evrs, vaxes, n_ref) for every pixel pickle of a mission."""
    comps, evrs, vaxes, n_ref = [], [], [], 0
    pat = os.path.join(components_dir, f"decomposition_{mission_id}_*_components.pkl")
    for p in sorted(glob.glob(pat)):
        try:
            with open(p, "rb") as fh:
                d = pickle.load(fh)
        except Exception as e:  # noqa: BLE001
            print(f"    ! could not read {os.path.basename(p)}: {e}")
            continue
        c = np.asarray(d.get("components"))
        if c is None or c.ndim != 2 or c.shape[0] == 0:
            continue
        comps.append(c)
        evrs.append(np.asarray(d.get("explained_variance_ratio")))
        va = d.get("metadata", {}).get("velocity_axis_kms")
        vaxes.append(np.asarray(va) if va is not None else None)
        n_ref += int(d.get("metadata", {}).get("n_reference_spectra", 0) or 0)
    return comps, evrs, vaxes, n_ref


def component_coherence(comps, vaxes, i):
    """Median |Pearson r| of component i across pixels, aligned on velocity.

    Pixels can be edge-trimmed to slightly different channel counts, so vectors
    are aligned on the velocities shared by ALL contributing pixels (via
    velocity_axis_kms) before correlating.
    """
    pairs = [(c[i], va) for c, va in zip(comps, vaxes)
             if c.shape[0] > i and va is not None and va.shape == c[i].shape]
    pairs = [(v, va) for v, va in pairs if np.isfinite(v).all() and np.std(v) > 0]
    if len(pairs) < 2:
        return np.nan, len(pairs)
    common = set.intersection(*[set(np.round(va, 2)) for _, va in pairs])
    if len(common) < 8:
        return np.nan, len(pairs)
    aligned = []
    for v, va in pairs:
        r = np.round(va, 2)
        mask = np.isin(r, list(common))
        order = np.argsort(r[mask])
        aligned.append(v[mask][order])
    L = min(len(a) for a in aligned)
    aligned = [a[:L] for a in aligned]
    R = np.corrcoef(np.vstack(aligned))
    iu = np.triu_indices(len(aligned), k=1)
    return float(np.median(np.abs(R[iu]))), len(pairs)


def choose_n_components(coh, threshold):
    """Leading contiguous components with coherence >= threshold (>=1)."""
    n = 0
    for c in coh:
        if np.isfinite(c) and c >= threshold:
            n += 1
        else:
            break
    return max(1, n)


def score_mission(components_dir, mission_id, threshold, min_pixels):
    """Return (n_components, coh_list, npix, n_ref) for a mission, printing a table."""
    comps, evrs, vaxes, n_ref = load_pixel_components(components_dir, mission_id)
    npix = len(comps)
    if npix == 0:
        print(f"{mission_id}: NO pickles found — skipped\n")
        return None
    max_nc = max(c.shape[0] for c in comps)
    print(f"{mission_id}   ({npix} pixels, {n_ref} ref spectra, up to {max_nc} comps)")
    coh = []
    for i in range(max_nc):
        c, npair = component_coherence(comps, vaxes, i)
        coh.append(c)
        evs = [e[i] for e in evrs if e is not None and len(e) > i]
        ev = float(np.mean(evs)) if evs else float("nan")
        flag = "keep" if (np.isfinite(c) and c >= threshold) else "noise"
        print(f"    C{i+1}: coh|r|={c:5.2f}  meanEVR={ev:6.3f}  ({npair} px)  {flag}")
    nc = choose_n_components(coh, threshold)
    if npix < min_pixels:
        print(f"    ⚠ only {npix} pixels (< {min_pixels}) — coherence unreliable, review by hand")
    print(f"    => n_components = {nc}\n")
    return nc, coh, npix, n_ref


# ── Decomposition driver (always re-run) ───────────────────────────────────────
def _keep(name, only, exclude):
    """Substring include/exclude filter (matches flight id or mission_id)."""
    if only and not any(s in name for s in only):
        return False
    if exclude and any(s in name for s in exclude):
        return False
    return True


def run_decompositions(reduced_dir, config, n_components, verbose, only, exclude,
                       file_glob="reduced_data_*.fits"):
    """Re-run pca_decompose for every input FITS in reduced_dir.

    file_glob selects which files to feed (default reduced_data_*.fits); pass e.g.
    "post_cleaned_nopca_*.fits" to score already post-processed/filtered inputs.
    Returns the wall-clock start time so the caller can restrict scoring to the
    pickles this run (re)wrote.
    """
    files = sorted(glob.glob(os.path.join(reduced_dir, file_glob)))
    files = [f for f in files if _keep(os.path.basename(f), only, exclude)]
    if not files:
        sys.exit(f"No {file_glob} found in {reduced_dir}"
                 + (" matching --only/--exclude" if (only or exclude) else ""))
    exe = shutil.which("pca_decompose")
    base = [exe] if exe else [sys.executable, "-m",
                              "oi_zeigt.pca_analysis.pca_decompose_per_mission"]
    t0 = time.time()
    print(f"Re-running decomposition for {len(files)} flight(s) "
          f"(n_components={n_components} ceiling)...\n")
    for i, f in enumerate(files, 1):
        print(f"  [{i}/{len(files)}] {os.path.basename(f)} ...", flush=True)
        cmd = base + ["--config", config, "--fits", f,
                      "--n-components", str(n_components)]
        res = subprocess.run(cmd, text=True,
                             capture_output=not verbose)
        if res.returncode != 0:
            if not verbose:
                sys.stderr.write((res.stdout or "") + (res.stderr or ""))
            sys.exit(f"pca_decompose failed on {f} (exit {res.returncode})")
    print()
    return t0


def discover_fresh_missions(components_dir, t0):
    """mission_ids whose pickles were written by this run (mtime >= t0)."""
    mids = set()
    for p in glob.glob(os.path.join(components_dir, "decomposition_*_components.pkl")):
        if os.path.getmtime(p) >= t0 - 1.0:
            m = PKL_RE.search(os.path.basename(p))
            if m:
                mids.add(m.group("mid"))
    return sorted(mids)


# ── yaml emit / merge ──────────────────────────────────────────────────────────
def _fmt(v):
    return str(v).lower() if isinstance(v, bool) else str(v)


def emit_yaml(results, gates):
    """Print a ready-to-paste yaml block for every scored mission."""
    print("# ---- yaml block (n_components from coherence; standard gates) ----")
    for mid, (nc, comment) in results.items():
        print(f"{mid}:")
        print(f"  n_components: {nc}" + (f"           # {comment}" if comment else ""))
        for k in GATE_ORDER:
            print(f"  {k}: {_fmt(gates[k])}")
        print()


def _flight_blocks(lines):
    """Yield (name, body_start, body_end) for each top-level flight entry.

    body excludes trailing blank/comment lines (those belong to the next entry
    or are section headers) so their comments are preserved untouched.
    """
    idxs = [i for i, l in enumerate(lines) if FLIGHT_KEY_RE.match(l)]
    for j, i in enumerate(idxs):
        name = lines[i].rstrip().rstrip(":")
        end = idxs[j + 1] if j + 1 < len(idxs) else len(lines)
        k = end - 1
        while k > i and (lines[k].strip() == "" or lines[k].lstrip().startswith("#")):
            k -= 1
        yield name, i + 1, k + 1


def _set_key(body, key, value, comment=None):
    """Set `  key: value` in a block body.

    If a new `comment` is given it is used; otherwise any inline comment already
    on the existing key line is preserved (so unchanged gate keys keep their
    rationale notes). Inserts the key if absent.
    """
    def _compose(text_comment):
        line = f"  {key}: {_fmt(value)}"
        if text_comment:
            line = line + " " * max(1, 26 - len(line)) + f"# {text_comment}"
        return line

    for idx, l in enumerate(body):
        if re.match(rf"^\s+{re.escape(key)}\s*:", l):
            existing = None
            if "#" in l:
                existing = l[l.index("#") + 1:].strip()
            body[idx] = _compose(comment if comment is not None else existing)
            return
    line = _compose(comment)
    if key == "n_components":
        body.insert(0, line)
    else:
        body.append(line)


def merge_into_yaml(path, results, gates, n_components_only=False):
    """In-place update of n_components + gate keys for scored flights.

    Preserves all comments and any flights not in `results`. Existing flights are
    updated key-by-key; flights not present are appended as a new block. When
    `n_components_only` is True, only the n_components line is touched and every
    existing gate key (and its hand-tuned value/comment) is left untouched — use
    this to refresh coherence-derived n_components without disturbing gates that
    were tuned by hand or that the gate-override flags can't express (e.g. a
    float line_cutoff_std).
    """
    text = Path(path).read_text()
    lines = text.split("\n")
    blocks = list(_flight_blocks(lines))
    present = {name for name, _, _ in blocks}
    # Update existing blocks in reverse so earlier slices stay valid.
    for name, bs, be in reversed(blocks):
        if name not in results:
            continue
        nc, comment = results[name]
        body = lines[bs:be]
        _set_key(body, "n_components", nc, comment)
        if not n_components_only:
            for k in GATE_ORDER:
                _set_key(body, k, gates[k])
        lines[bs:be] = body
    # Append flights not already in the file.
    for name, (nc, comment) in results.items():
        if name in present:
            continue
        block = [f"{name}:", f"  n_components: {_fmt(nc)}"
                 + (f"           # {comment}" if comment else "")]
        if not n_components_only:
            block += [f"  {k}: {_fmt(gates[k])}" for k in GATE_ORDER]
        if lines and lines[-1].strip() != "":
            lines.append("")
        lines.extend(block)
    Path(path).write_text("\n".join(lines))
    updated = sum(1 for n in results if n in present)
    appended = len(results) - updated
    print(f"✓ wrote {path}: {updated} entr{'y' if updated == 1 else 'ies'} updated, "
          f"{appended} appended")


def main_cli():
    parser = argparse.ArgumentParser(
        description="Recommend per-flight PCA n_components from a directory of "
                    "reduced_data_*.fits, by re-running the per-pixel decomposition "
                    "and scoring pixel-to-pixel coherence.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("--dir", required=True,
                        help="Directory containing the input FITS (see --file-glob)")
    parser.add_argument("--file-glob", default="reduced_data_*.fits", dest="file_glob",
                        help="Filename glob selecting inputs inside --dir. Default: "
                             "reduced_data_*.fits. Use e.g. 'post_cleaned_nopca_*.fits' "
                             "to score already post-processed/filtered data.")
    parser.add_argument("--config", default="config.toml",
                        help="config.toml (for pca_source, smoothing, components_dir, "
                             "mission_pca_parameters). Default: config.toml")
    parser.add_argument("--n-components", type=int, default=5, dest="n_components",
                        help="Ceiling of components to EXTRACT during decomposition "
                             "(coherence then picks how many to keep). Default: 5")
    parser.add_argument("--coherence-threshold", type=float, default=0.40,
                        dest="threshold",
                        help="Keep leading components with median |r| >= this. Default: 0.40")
    parser.add_argument("--min-pixels", type=int, default=3, dest="min_pixels",
                        help="Warn if a flight has fewer pixels than this. Default: 3")
    parser.add_argument("--only", nargs="+", default=None,
                        help="Restrict to flights whose reduced filename / mission_id "
                             "contains any of these substrings (e.g. --only 352 372).")
    parser.add_argument("--exclude", nargs="+", default=None,
                        help="Skip flights whose reduced filename / mission_id contains "
                             "any of these substrings (e.g. --exclude 296 407).")
    parser.add_argument("--components-dir", default=None, dest="components_dir",
                        help="Where decomposition pickles are read from "
                             "(default: config [output][components_dir])")
    parser.add_argument("--skip-decompose", action="store_true",
                        help="Do NOT re-run decomposition; score existing pickles "
                             "in components-dir. (Default is to always re-decompose.)")
    parser.add_argument("--write", action="store_true",
                        help="Merge n_components + standard gates into the params yaml "
                             "(config [input][mission_pca_parameters]), preserving comments.")
    parser.add_argument("--params-file", default=None, dest="params_file",
                        help="Params yaml to update with --write "
                             "(default: config [input][mission_pca_parameters])")
    parser.add_argument("--n-components-only", action="store_true",
                        dest="n_components_only",
                        help="With --write, update ONLY the n_components line for each "
                             "scored flight and leave every existing gate key untouched "
                             "(preserves hand-tuned gates/comments).")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="Stream pca_decompose output instead of capturing it")
    # Gate overrides (defaults = the standard tuned block)
    for key, val in DEFAULT_GATES.items():
        flag = "--" + key.replace("_", "-")
        if isinstance(val, bool):
            parser.add_argument(flag, dest=f"gate_{key}", default=None,
                                choices=["true", "false"],
                                help=f"Gate override (default: {str(val).lower()})")
        else:
            typ = type(val)
            parser.add_argument(flag, dest=f"gate_{key}", type=typ, default=None,
                                help=f"Gate override (default: {val})")
    args = parser.parse_args()

    reduced_dir = os.path.abspath(args.dir)
    if not os.path.isdir(reduced_dir):
        sys.exit(f"--dir not a directory: {reduced_dir}")

    # Resolve config-derived paths
    from oi_zeigt.basic_io import get_config
    config = get_config(args.config)
    out_cfg = config.get("output", {})
    components_dir = args.components_dir or out_cfg.get("components_dir", "output/pca_components")
    components_dir = os.path.abspath(components_dir)

    # Assemble gate values (defaults + overrides)
    gates = dict(DEFAULT_GATES)
    for key in DEFAULT_GATES:
        ov = getattr(args, f"gate_{key}")
        if ov is not None:
            gates[key] = (ov == "true") if isinstance(DEFAULT_GATES[key], bool) else ov

    print(f"reduced data dir : {reduced_dir}")
    print(f"components dir    : {components_dir}")
    print(f"coherence thresh : |r| >= {args.threshold}\n")

    # 1) (Re)decompose, then figure out which missions were produced
    if args.skip_decompose:
        t0 = 0.0
    else:
        t0 = run_decompositions(reduced_dir, args.config, args.n_components,
                                args.verbose, args.only, args.exclude,
                                file_glob=args.file_glob)
    missions = discover_fresh_missions(components_dir, t0)
    missions = [m for m in missions if _keep(m, args.only, args.exclude)]
    if not missions:
        sys.exit(f"No decomposition pickles found in {components_dir}"
                 + (" matching --only/--exclude" if (args.only or args.exclude) else ""))

    # 2) Score each mission by coherence
    results = {}  # mission_id -> (n_components, comment)
    for mid in missions:
        scored = score_mission(components_dir, mid, args.threshold, args.min_pixels)
        if scored is None:
            continue
        nc, coh, npix, n_ref = scored
        parts = "; ".join(f"C{i+1}={c:.2f}" for i, c in enumerate(coh) if np.isfinite(c))
        comment = f"coherence {parts}"
        results[mid] = (nc, comment)

    # 3) Emit / write
    print()
    emit_yaml(results, gates)
    if args.write:
        params_file = args.params_file or config.get("input", {}).get("mission_pca_parameters")
        if not params_file:
            sys.exit("--write: no params file (set --params-file or config "
                     "[input][mission_pca_parameters])")
        merge_into_yaml(params_file, results, gates, n_components_only=args.n_components_only)


if __name__ == "__main__":
    main_cli()
