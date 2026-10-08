"""
Observed vs fitted sky: S-H_OBS spectra overlaid with their S-H_FIT counterparts.

Each S-H_OBS row is paired with the S-H_FIT row of the same TELESCOP, SCAN
and SUBSCAN.  One page per scan.subscan, one panel per receiver (TELESCOP),
each panel showing OBS and FIT overlaid above an OBS - FIT residual strip.

A *cycle* is one calibration subscan, named SCAN.SUBSCAN (e.g. 39658.1).
In a GREAT scan, calibration subscans (1, 5, 9, 13: sky, hot and cold loads,
giving S-H_OBS / S-H_FIT) alternate with science subscans (2, 6, 10, 14);
each calibration subscan presumably calibrates the science subscan that
follows it.  kalibrate fits the atmosphere once per cycle, jointly for all
pixels, so MH2O and CHI_SQR have one value per cycle, copied onto every row.
"""

from pathlib import Path
from typing import Callable, Dict, List, Optional, Sequence, Tuple, Union

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from astropy.io import fits

from .spectra_plots import _column_units, _find_spectrum_hdu, _to_str

PairKey = Tuple[str, int, int]  # (TELESCOP, SCAN, SUBSCAN)
CycleKey = Tuple[int, int]  # (SCAN, SUBSCAN): one calibration cycle, see module docstring


def parse_cycle(text: str) -> CycleKey:
    """'39658.1' -> (39658, 1).  Parsed as text, so '39646.10' is subscan 10, not 1."""
    try:
        scan, sub = text.strip().split('.')
        return int(scan), int(sub)
    except ValueError:
        raise ValueError(f"cycle must be SCAN.SUBSCAN, e.g. 39658.1; got '{text}'")


def _select_cycles(keys, scans, subscans, cycles) -> np.ndarray:
    """Boolean mask over (scan, subscan) keys for the --scan/--subscan/--cycle selection."""
    keys = list(keys)
    keep = np.ones(len(keys), dtype=bool)
    for i, (scan, sub) in enumerate(keys):
        if scans and scan not in scans:
            keep[i] = False
        elif subscans and sub not in subscans:
            keep[i] = False
        elif cycles and (scan, sub) not in cycles:
            keep[i] = False
    return keep


def pair_rows(data, obs_object: str = 'S-H_OBS',
              fit_object: str = 'S-H_FIT') -> Tuple[Dict[PairKey, Tuple[int, int]], List[PairKey], List[PairKey]]:
    """
    Pair ``obs_object`` rows with ``fit_object`` rows on (TELESCOP, SCAN, SUBSCAN).

    Returns ({key: (obs_row, fit_row)}, obs keys without a fit, fit keys without an obs).
    A key that occurs more than once for the same OBJECT raises ValueError,
    since the pairing would then be ambiguous.
    """
    objects = np.array([_to_str(o).upper() for o in data['OBJECT']])
    rows_by_object = {}
    for name in (obs_object, fit_object):
        rows = {}
        for i in np.flatnonzero(objects == name.upper()):
            key = (_to_str(data['TELESCOP'][i]).upper(), int(data['SCAN'][i]), int(data['SUBSCAN'][i]))
            if key in rows:
                raise ValueError(f"{name}: two rows ({rows[key]}, {i}) for TELESCOP/SCAN/SUBSCAN {key}")
            rows[key] = i
        if not rows:
            raise ValueError(f"No spectra with OBJECT '{name}'")
        rows_by_object[name] = rows
    obs, fit = rows_by_object[obs_object], rows_by_object[fit_object]
    pairs = {k: (obs[k], fit[k]) for k in sorted(obs.keys() & fit.keys())}
    return pairs, sorted(obs.keys() - fit.keys()), sorted(fit.keys() - obs.keys())


def _velocity(hdu, row: int, nchan: int) -> np.ndarray:
    """Velocity axis (km/s) of one row, from its own VELOCITY/DELTAV and CRPIX1."""
    crpix1 = float(hdu.header.get('CRPIX1', 1.0))
    ch = np.arange(nchan, dtype=float) - (crpix1 - 1.0)
    return (float(hdu.data['VELOCITY'][row]) + ch * float(hdu.data['DELTAV'][row])) / 1000.0


def _page(hdu, title: str, telescops: List[str], pairs: Dict[PairKey, Tuple[int, int]],
          scan: int, subscan: int, column: str, ncols: int,
          velocity_range: Optional[Tuple[float, float]], labels: Tuple[str, str]) -> plt.Figure:
    """One figure: a panel per TELESCOP (OBS/FIT overlay + residual strip) for one scan.subscan."""
    ylabel, unit = _column_units(column)
    ncols = max(1, min(ncols, len(telescops)))
    nrows = -(-len(telescops) // ncols)
    fig = plt.figure(figsize=(4.5 * ncols, 3.6 * nrows + 1.0), layout='constrained')
    gs = fig.add_gridspec(2 * nrows, ncols, height_ratios=[3, 1] * nrows)
    has_chi = 'CHI_SQR' in hdu.data.dtype.names

    for k, tel in enumerate(telescops):
        r, c = divmod(k, ncols)
        ax = fig.add_subplot(gs[2 * r, c])
        ax_res = fig.add_subplot(gs[2 * r + 1, c], sharex=ax)
        ax.tick_params(labelbottom=False, labelsize=7)
        ax_res.tick_params(labelsize=7)
        ax.grid(True, alpha=0.3)
        ax_res.grid(True, alpha=0.3)
        if (tel, scan, subscan) not in pairs:
            ax.set_title(f"{tel}: no pair", fontsize=8, loc='left')
            continue
        i_obs, i_fit = pairs[(tel, scan, subscan)]
        obs = np.asarray(hdu.data[column][i_obs], dtype=float)
        fit = np.asarray(hdu.data[column][i_fit], dtype=float)
        vel = _velocity(hdu, i_obs, len(obs))
        res = obs - fit
        shown = np.ones(len(vel), dtype=bool) if velocity_range is None else \
            (vel >= min(velocity_range)) & (vel <= max(velocity_range))
        finite = np.isfinite(res[shown])
        rms = float(np.std(res[shown][finite])) if finite.any() else np.nan

        # Rasterized: 16k-channel curves x 14 panels x tens of pages would bloat a vector PDF.
        ax.plot(vel, obs, linewidth=0.5, color='tab:blue', label=labels[0], rasterized=True)
        ax.plot(vel, fit, linewidth=0.8, color='tab:red', label=labels[1], rasterized=True)
        ax_res.plot(vel, res, linewidth=0.4, color='k', rasterized=True)
        ax_res.axhline(0, color='tab:red', linewidth=0.6)
        chi = f"  χ²={float(hdu.data['CHI_SQR'][i_obs]):.3g}" if has_chi else ""
        ax.set_title(f"{tel}{chi}", fontsize=8, loc='left')
        ax_res.text(0.99, 0.95, f"rms {rms:.3g}{unit}" if finite.any() else "no residual",
                    transform=ax_res.transAxes, ha='right', va='top', fontsize=7)
        for name, spec in zip(labels, (obs, fit)):
            if not np.isfinite(spec).any():
                ax.text(0.5, 0.5, f"{name} all NaN", transform=ax.transAxes, ha='center',
                        va='center', fontsize=9, color='tab:red')
        if velocity_range is not None:
            ax.set_xlim(min(velocity_range), max(velocity_range))
            for a, y in ((ax, np.r_[obs[shown], fit[shown]]), (ax_res, res[shown])):
                y = y[np.isfinite(y)]
                if y.size:
                    pad = 0.05 * (y.max() - y.min() or 1.0)
                    a.set_ylim(y.min() - pad, y.max() + pad)
        if c == 0:
            ax.set_ylabel(ylabel, fontsize=8)
            ax_res.set_ylabel('OBS−FIT', fontsize=7)
        if k + ncols >= len(telescops):  # lowest panel in its column
            ax_res.set_xlabel('Velocity (km/s)', fontsize=8)

    handles = [plt.Line2D([], [], color='tab:blue', label=labels[0]),
               plt.Line2D([], [], color='tab:red', label=labels[1]),
               plt.Line2D([], [], color='k', label=f"{labels[0]} − {labels[1]}")]
    fig.legend(handles=handles, loc='outside upper center', ncol=3, fontsize=9,
               frameon=False, title=title, title_fontsize=12)
    return fig


def plot_sky_fit_pairs(input_fits: Union[str, Path],
                       output: Optional[Union[str, Path]] = None,
                       obs_object: str = 'S-H_OBS',
                       fit_object: str = 'S-H_FIT',
                       column: str = 'SPECTRUM',
                       scans: Optional[Sequence[int]] = None,
                       subscans: Optional[Sequence[int]] = None,
                       telescops: Optional[Sequence[str]] = None,
                       cycles: Optional[Sequence[CycleKey]] = None,
                       velocity_range: Optional[Tuple[float, float]] = None,
                       ncols: int = 5,
                       max_screen_pages: int = 5,
                       echo: Callable[[str], None] = print) -> int:
    """
    Plot every observed sky spectrum with its fitted sky, one page per scan.subscan.

    ``obs_object`` rows are paired with ``fit_object`` rows on TELESCOP, SCAN
    and SUBSCAN (see :func:`pair_rows`).  Each page has one panel per
    TELESCOP (``ncols`` per row): OBS and FIT overlaid, and below it the
    OBS - FIT residual with its rms (over ``velocity_range`` when given).

    ``scans``, ``subscans``, ``cycles`` ((scan, subscan) tuples, see the
    module docstring) and ``telescops`` restrict what is plotted; requested
    cycles that are not in the file are reported.
    With ``output`` (a .pdf), all pages go into one multi-page PDF; without
    it, pages are shown on screen, which is refused for more than
    ``max_screen_pages`` pages.  Returns the number of pages.
    """
    with fits.open(input_fits) as hdul:
        hdu = _find_spectrum_hdu(hdul, column)
        for col in ('OBJECT', 'TELESCOP', 'SCAN', 'SUBSCAN', 'VELOCITY', 'DELTAV'):
            if col not in hdu.data.dtype.names:
                raise ValueError(f"{input_fits}: no {col} column")
        pairs, obs_only, fit_only = pair_rows(hdu.data, obs_object, fit_object)
        echo(f"{input_fits}: {len(pairs)} {obs_object}/{fit_object} pairs")
        for name, missing, other in ((obs_object, obs_only, fit_object), (fit_object, fit_only, obs_object)):
            if missing:
                echo(f"  {len(missing)} {name} spectra have no {other}, e.g. {missing[:3]}")

        if telescops:
            wanted = {t.upper() for t in telescops}
            pairs = {k: v for k, v in pairs.items() if k[0] in wanted}
        if cycles:
            missing = sorted(set(cycles) - {(scan, sub) for _, scan, sub in pairs})
            if missing:
                echo(f"  Cycle(s) not in this file: {', '.join(f'{s}.{ss}' for s, ss in missing)}")
        keys = list(pairs)
        keep = _select_cycles([(scan, sub) for _, scan, sub in keys], scans, subscans, cycles)
        pairs = {k: pairs[k] for k, kp in zip(keys, keep) if kp}
        if not pairs:
            raise ValueError(f"{input_fits}: no pairs left after the --scan/--subscan/--cycle/--telescop selection")

        pages = sorted({(scan, sub) for _, scan, sub in pairs})
        panel_tels = sorted({tel for tel, _, _ in pairs})
        if output is None and len(pages) > max_screen_pages:
            raise ValueError(f"{len(pages)} pages is too many to show on screen; give --output-dir "
                             f"or narrow with --cycle/--scan/--subscan")
        echo(f"  {len(pages)} page(s), one per cycle (scan.subscan), {len(panel_tels)} receivers per page")

        labels = (obs_object, fit_object)
        name = Path(input_fits).name
        if output is None:
            for scan, sub in pages:
                _page(hdu, f"{name}: cycle {scan}.{sub}", panel_tels, pairs, scan, sub,
                      column, ncols, velocity_range, labels)
            return len(pages)

        output = Path(output)
        output.parent.mkdir(parents=True, exist_ok=True)
        with PdfPages(output) as pdf:
            for scan, sub in pages:
                fig = _page(hdu, f"{name}: cycle {scan}.{sub}", panel_tels, pairs, scan, sub,
                            column, ncols, velocity_range, labels)
                pdf.savefig(fig, dpi=110)
                plt.close(fig)
        echo(f"  Saved {output}")
    return len(pages)


def chi_sqr_by_cycle(input_fits: Union[str, Path], obs_object: str = 'S-H_OBS',
                     scans: Optional[Sequence[int]] = None,
                     subscans: Optional[Sequence[int]] = None,
                     cycles: Optional[Sequence[CycleKey]] = None) -> Tuple[Dict[CycleKey, float], List[CycleKey]]:
    """
    CHI_SQR of the sky fit per calibration cycle (scan.subscan), from the ``obs_object`` rows.

    kalibrate writes one CHI_SQR per scan.subscan onto every row of that
    cycle.  Returns ({(scan, subscan): chi_sqr}, cycles whose rows disagree);
    a disagreeing cycle gets the mean of its values.
    """
    with fits.open(input_fits) as hdul:
        hdu = _find_spectrum_hdu(hdul, 'CHI_SQR')
        data = hdu.data
        objects = np.array([_to_str(o).upper() for o in data['OBJECT']])
        rows = np.flatnonzero(objects == obs_object.upper())
        if len(rows) == 0:
            raise ValueError(f"{input_fits}: no spectra with OBJECT '{obs_object}'")
        scan = np.asarray(data['SCAN'][rows], dtype=int)
        sub = np.asarray(data['SUBSCAN'][rows], dtype=int)
        chi = np.asarray(data['CHI_SQR'][rows], dtype=float)

    keep = _select_cycles(zip(scan.tolist(), sub.tolist()), scans, subscans, cycles)
    values: Dict[CycleKey, List[float]] = {}
    for s, ss, c in zip(scan[keep], sub[keep], chi[keep]):
        values.setdefault((int(s), int(ss)), []).append(float(c))
    result = {k: float(np.mean(v)) for k, v in sorted(values.items())}
    mixed = [k for k, v in values.items() if np.ptp(v) > 1e-6 * max(1.0, abs(np.mean(v)))]
    return result, sorted(mixed)


def chi_sqr_report(input_fits: Sequence[Union[str, Path]], obs_object: str = 'S-H_OBS',
                   scans: Optional[Sequence[int]] = None,
                   subscans: Optional[Sequence[int]] = None,
                   cycles: Optional[Sequence[CycleKey]] = None,
                   per_cycle: bool = False) -> str:
    """
    Text report comparing the sky-fit CHI_SQR of several files.

    Per file: number of cycles and median / mean / min / max CHI_SQR.  With two
    or more files, a paired comparison over the cycles they share: how often
    each file has the lowest CHI_SQR and the median difference, and which file
    fits best.  ``per_cycle`` appends a table of every cycle's value per file.
    """
    names = [Path(p).name for p in input_fits]
    tags = [chr(ord('A') + i) if i < 26 else f"F{i}" for i in range(len(names))]
    per_file = []
    for path in input_fits:
        chis, mixed = chi_sqr_by_cycle(path, obs_object, scans, subscans, cycles)
        if not chis:
            raise ValueError(f"{path}: no {obs_object} rows left after the scan/subscan/cycle selection")
        per_file.append((chis, mixed))

    lines = [f"CHI_SQR of the sky fit — one value per calibration cycle (scan.subscan), "
             f"from the {obs_object} rows. Lower is better.",
             "A cycle is one calibration subscan (e.g. 39658.1: scan 39658, subscan 1). Calibration",
             "subscans 1, 5, 9, 13 alternate with science subscans 2, 6, 10, 14; kalibrate fits the",
             "atmosphere once per cycle for all pixels together, so CHI_SQR is per cycle, not per pixel.",
             ""]
    w = max(len(n) for n in names)
    lines.append(f"    {'file':<{w}}  cycles   median     mean      min      max")
    medians = []
    for tag, name, (chis, mixed) in zip(tags, names, per_file):
        v = np.array(list(chis.values()))
        medians.append(float(np.median(v)))
        lines.append(f"{tag}:  {name:<{w}}  {len(v):>6}  {np.median(v):>7.3f}  {v.mean():>7.3f}  "
                     f"{v.min():>7.3f}  {v.max():>7.3f}")
        if mixed:
            lines.append(f"    note: {len(mixed)} cycle(s) have differing CHI_SQR across rows "
                         f"(mean used), e.g. {mixed[:3]}")

    if len(per_file) > 1:
        common = sorted(set.intersection(*(set(c) for c, _ in per_file)))
        lines.append("")
        if not common:
            lines.append("No scan.subscan cycles in common: cannot compare cycle by cycle.")
        else:
            table = np.array([[chis[k] for chis, _ in per_file] for k in common])
            best = table.argmin(axis=1)
            ties = (table == table.min(axis=1, keepdims=True)).sum(axis=1) > 1
            lines.append(f"Cycle by cycle over the {len(common)} common cycles "
                         f"(lowest CHI_SQR per cycle):")
            for i, tag in enumerate(tags):
                lines.append(f"    {tag}: lowest in {int(((best == i) & ~ties).sum())}/{len(common)}")
            if ties.any():
                lines.append(f"    ties: {int(ties.sum())}")
            if len(per_file) == 2:
                diff = table[:, 0] - table[:, 1]
                lines.append(f"    median difference A − B: {np.median(diff):+.3f} "
                             f"(mean {diff.mean():+.3f})")
            wins = np.bincount(best[~ties], minlength=len(tags))
            by_median, by_wins = int(np.argmin(medians)), int(np.argmax(wins))
            lines.append("")
            if by_median == by_wins:
                lines.append(f"Best fits: {tags[by_median]} ({names[by_median]}) — lowest median "
                             f"CHI_SQR and lowest in most cycles.")
            else:
                lines.append(f"Mixed result: {tags[by_median]} has the lowest median CHI_SQR, "
                             f"but {tags[by_wins]} is lowest in most cycles.")

    if per_cycle:
        cycles = sorted(set().union(*(set(c) for c, _ in per_file)))
        lines += ["", "Per cycle:", "    scan.subscan  " + "  ".join(f"{t:>8}" for t in tags)
                  + ("    lowest" if len(tags) > 1 else "")]
        for k in cycles:
            vals = [chis.get(k) for chis, _ in per_file]
            cells = "  ".join(f"{v:>8.3f}" if v is not None else f"{'—':>8}" for v in vals)
            present = [(v, t) for v, t in zip(vals, tags) if v is not None]
            low = ""
            if len(tags) > 1 and len(present) > 1:
                lowest = min(v for v, _ in present)
                winners = [t for v, t in present if v == lowest]
                low = winners[0] if len(winners) == 1 else "tie"
            lines.append(f"    {f'{k[0]}.{k[1]}':<12}  {cells}    {low}".rstrip())
    return "\n".join(lines)
