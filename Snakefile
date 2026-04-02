# oi_zeigt pipeline — Snakemake workflow
#
# Usage (from the project directory, with ~/.venv/oi active):
#
#   snakemake --cores 8          # run full pipeline, 8 parallel jobs
#   snakemake --cores 8 -np      # dry-run: show what would execute
#   snakemake --cores 8 merge_chunks   # stop before create_datacube
#
# Two environments are used automatically:
#   ~/.venv/oi      — steps 0–6 (split, filter, reduce, prepare, decompose, correct, merge)
#   ~/.venv/cygrid  — step 7   (create_datacube / gridding)
#
# The split step produces two manifests:
#   manifest.txt         — all 1106 chunks (including TREC-only scans)
#   manifest_science.txt — 518 science chunks only (skip TREC-only in pipeline)
#
# All paths come from config.toml.  Override with:
#   snakemake --config config_file=other_config.toml

import tomllib
from pathlib import Path

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
CONFIG_FILE = config.get("config_file", "config.toml")

# Full paths to executables — correct environment used per step regardless
# of which venv is active when snakemake is launched.
OI     = lambda cmd: str(Path("~/.venv/oi/bin").expanduser()     / cmd)
CYGRID = lambda cmd: str(Path("~/.venv/cygrid/bin").expanduser() / cmd)

with open(CONFIG_FILE, "rb") as _f:
    CFG = tomllib.load(_f)

FITS_INPUT   = CFG["input"]["fits_file"]
WORK_DIR     = Path("pipeline_work")
CHUNKS_DIR   = WORK_DIR / "chunks"
CLEAN_DIR    = WORK_DIR / "clean"
REDUCED_DIR  = WORK_DIR / "reduced"
PREPARED_DIR = WORK_DIR / "prepared"
DECOMP_DIR   = WORK_DIR / "decompositions"
PCAD_DIR     = WORK_DIR / "pcad"

DATACUBE_OUT = CFG["output"]["datacube"]

# ---------------------------------------------------------------------------
# Helpers: read chunk IDs from the manifests after the split checkpoint
# ---------------------------------------------------------------------------
def _read_manifest(path):
    with open(path) as fh:
        return [l.strip() for l in fh if l.strip() and not l.startswith("#")]

def _science_chunk_ids(wildcards):
    """Chunk stems for science scans only (excludes TREC-only)."""
    checkpoints.split_fits.get(**wildcards)
    return [Path(p).stem for p in _read_manifest(CHUNKS_DIR / "manifest_science.txt")]

def expand_science(wildcards, subdir, suffix):
    return expand(str(subdir / "{chunk}") + suffix, chunk=_science_chunk_ids(wildcards))

# ---------------------------------------------------------------------------
# Default target
# ---------------------------------------------------------------------------
rule all:
    input:
        DATACUBE_OUT,

# ---------------------------------------------------------------------------
# Step 0 — split the big file into per-(mission, telescope, scan) chunks
#           Two manifests: all chunks and science-only chunks
# ---------------------------------------------------------------------------
checkpoint split_fits:
    input:
        fits = FITS_INPUT,
    output:
        manifest         = str(CHUNKS_DIR / "manifest.txt"),
        manifest_science = str(CHUNKS_DIR / "manifest_science.txt"),
    params:
        outdir = str(CHUNKS_DIR),
        exe    = OI("split_fits"),
    shell:
        "{params.exe} --fits {input.fits} --output-dir {params.outdir}"
        " --manifest {output.manifest}"
        " --manifest-science {output.manifest_science}"

# ---------------------------------------------------------------------------
# Step 1 — filter each science chunk  (~/.venv/oi)
# ---------------------------------------------------------------------------
rule filter_chunk:
    input:
        fits   = str(CHUNKS_DIR / "{chunk}.fits"),
        config = CONFIG_FILE,
    output:
        clean  = str(CLEAN_DIR / "{chunk}_clean.fits"),
    params:
        exe = OI("filter_fits"),
    shell:
        "{params.exe} --fits {input.fits} --config {input.config} --output-clean {output.clean}"

# ---------------------------------------------------------------------------
# Step 2 — reduce spectra  (~/.venv/oi)
# ---------------------------------------------------------------------------
rule reduce_chunk:
    input:
        fits   = str(CLEAN_DIR / "{chunk}_clean.fits"),
        config = CONFIG_FILE,
    output:
        reduced = str(REDUCED_DIR / "{chunk}_reduced.fits"),
    params:
        exe = OI("reduce_spectra"),
    shell:
        "{params.exe} --fits {input.fits} --config {input.config} --output {output.reduced}"

# ---------------------------------------------------------------------------
# Step 3 — prepare for PCA  (~/.venv/oi)
# ---------------------------------------------------------------------------
rule prepare_pca_chunk:
    input:
        fits   = str(REDUCED_DIR / "{chunk}_reduced.fits"),
        config = CONFIG_FILE,
    output:
        prepared = str(PREPARED_DIR / "{chunk}_prepared.fits"),
    params:
        exe = OI("prepare_for_pca"),
    shell:
        "{params.exe} --fits {input.fits} --config {input.config} --output {output.prepared}"

# ---------------------------------------------------------------------------
# Step 4 — PCA decomposition  (~/.venv/oi)
# ---------------------------------------------------------------------------
rule pca_decompose_chunk:
    input:
        fits   = str(PREPARED_DIR / "{chunk}_prepared.fits"),
        config = CONFIG_FILE,
    output:
        decomp_dir = directory(str(DECOMP_DIR / "{chunk}")),
    params:
        exe = OI("pca_decompose"),
    shell:
        "{params.exe} --config {input.config} --fits {input.fits} --output-dir {output.decomp_dir}"

# ---------------------------------------------------------------------------
# Step 5 — PCA correction  (~/.venv/oi)
# ---------------------------------------------------------------------------
rule pca_correct_chunk:
    input:
        fits       = str(PREPARED_DIR / "{chunk}_prepared.fits"),
        decomp_dir = str(DECOMP_DIR   / "{chunk}"),
        config     = CONFIG_FILE,
    output:
        pcad = str(PCAD_DIR / "{chunk}_pcad.fits"),
    params:
        exe = OI("pca_correct"),
    shell:
        "{params.exe} --config {input.config} --input {input.fits}"
        " --decomposition {input.decomp_dir}"
        " --output {output.pcad}"

# ---------------------------------------------------------------------------
# Step 6 — merge all corrected science chunks  (~/.venv/oi)
# ---------------------------------------------------------------------------
rule merge_chunks:
    input:
        pcad_files = lambda wc: expand_science(wc, PCAD_DIR, "_pcad.fits"),
    output:
        manifest = str(WORK_DIR / "pcad_manifest.txt"),
        merged   = str(WORK_DIR / "merged_pcad.fits"),
    params:
        combine = OI("combine_fits"),
    run:
        with open(output.manifest, "w") as fh:
            for p in input.pcad_files:
                fh.write(str(p) + "\n")
        shell("{params.combine} --input {output.manifest} --output {output.merged} --single-hdu")

# ---------------------------------------------------------------------------
# Step 7 — build the data cube  (~/.venv/cygrid)
# ---------------------------------------------------------------------------
rule create_datacube:
    input:
        fits   = str(WORK_DIR / "merged_pcad.fits"),
        config = CONFIG_FILE,
    output:
        cube = DATACUBE_OUT,
    params:
        exe = CYGRID("create_datacube"),
    shell:
        "{params.exe} --fits {input.fits} --config {input.config} --output {output.cube}"
