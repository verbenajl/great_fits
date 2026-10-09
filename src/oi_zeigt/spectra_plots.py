"""
Averaged-spectrum plots per receiver, polarization or pixel, and a
per-scan / subscan / pixel diagnostic (``plot_grouped_spectra``).

GREAT TELESCOP values encode receiver, polarization and pixel, e.g.
``LFAH_PX03_S`` = receiver LFA, polarization H, pixel 03.  Spectra are grouped
on those fields and averaged channel by channel (NaN-aware), either as a plain
mean or weighted by 1/RMS² per spectrum.
"""

import re
from pathlib import Path
from typing import Callable, Dict, List, Optional, Tuple, Union

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

from .basic_io import reconstruct_velocity_axis, row_crpix1


# <receiver><pol>_PX<nn>[_<suffix>], e.g. LFAH_PX00_S, HFAV_PX06_S.
# The receiver is the name minus its last letter (LFAH -> LFA, HFAV -> HFA).
_TELESCOP_RE = re.compile(r'^(?P<receiver>[A-Z0-9]+?)(?P<pol>[HVI])_PX(?P<pixel>\d+)(?:_\w+)?$')

# OBJECT values (prefixes) of housekeeping rows, never treated as science.
# HOT-COLD rows carry SOBSMODE = 'ON', so they must be excluded by name.
_CALIBRATION_OBJECTS = ('TSYS', 'TAU_SIG', 'SKY-DIFF', 'SKYCHOPDIFF', 'TREC', 'S-H', 'HOT-COLD')

# Columns tried, in order, for the per-spectrum noise used by weighted averaging.
WEIGHT_RMS_COLUMNS = ('RMS', 'RMS_BASELINE')


def _to_str(v):
    return v.decode().strip() if isinstance(v, bytes) else str(v).strip()


def parse_telescop(name: str) -> Optional[Tuple[str, str, int]]:
    """
    Split a GREAT TELESCOP value into (receiver, polarization, pixel).

    >>> parse_telescop('LFAH_PX03_S')
    ('LFA', 'H', 3)

    Returns None if the value does not follow the ``<rx><pol>_PX<nn>`` pattern.
    """
    m = _TELESCOP_RE.match(_to_str(name).upper())
    if m is None:
        return None
    return m.group('receiver'), m.group('pol'), int(m.group('pixel'))


def select_science_rows(data, object_filter: Optional[str] = None) -> np.ndarray:
    """
    Boolean mask of science spectra, optionally restricted to an OBJECT substring.

    Science rows are SOBSMODE == 'ON' when that column exists, and in every case
    rows whose OBJECT is a known calibration/housekeeping name are dropped.
    """
    names = data.dtype.names
    objects = np.array([_to_str(o).upper() for o in data['OBJECT']]) if 'OBJECT' in names else None

    if 'SOBSMODE' in names:
        mask = np.array([_to_str(s).upper() == 'ON' for s in data['SOBSMODE']])
    else:
        mask = np.ones(len(data), dtype=bool)
    if objects is not None:
        mask &= ~np.array([any(o.startswith(c) for c in _CALIBRATION_OBJECTS) for o in objects])

    if object_filter and objects is not None:
        mask &= np.array([object_filter.upper() in o for o in objects])
    return mask


def average_group(spectra: np.ndarray, weights: Optional[np.ndarray] = None,
                  noise_mask: Optional[np.ndarray] = None) -> dict:
    """
    NaN-aware average of a (n_spectra, n_channels) block.

    With ``weights`` (one per spectrum), each channel is Σ(w·T)/Σ(w) over the
    spectra that are finite in that channel.  The reported rms is the scatter
    of the average over ``noise_mask`` channels (all channels if None).
    """
    finite = np.isfinite(spectra)
    if weights is None:
        avg = np.nanmean(spectra, axis=0)
    else:
        w = np.where(finite, weights[:, None], 0.0)
        wsum = w.sum(axis=0)
        with np.errstate(invalid='ignore', divide='ignore'):
            avg = np.where(wsum > 0, (w * np.where(finite, spectra, 0.0)).sum(axis=0) / wsum, np.nan)
    return {
        'avg_spectrum': avg,
        'count': len(spectra),
        'rms': float(np.nanstd(avg if noise_mask is None else avg[noise_mask])),
    }


def average_by(spectra: np.ndarray, keys: List[Optional[tuple]],
               weights: Optional[np.ndarray] = None,
               noise_mask: Optional[np.ndarray] = None) -> Dict[tuple, dict]:
    """Average ``spectra`` per distinct key; rows with key None are skipped."""
    keys_arr = np.empty(len(keys), dtype=object)
    keys_arr[:] = keys
    result = {}
    for key in sorted({k for k in keys if k is not None}):
        sel = np.array([k == key for k in keys_arr])
        result[key] = average_group(spectra[sel], None if weights is None else weights[sel], noise_mask)
    return result


def _rms_weights(data) -> Tuple[np.ndarray, str]:
    """1/RMS² per row from the first available RMS column; invalid RMS -> NaN."""
    for col in WEIGHT_RMS_COLUMNS:
        if col in data.dtype.names:
            rms = np.asarray(data[col], dtype=float)
            with np.errstate(divide='ignore', invalid='ignore'):
                w = np.where(np.isfinite(rms) & (rms > 0), 1.0 / rms**2, np.nan)
            return w, col
    raise ValueError(f"--weighted needs one of the columns {WEIGHT_RMS_COLUMNS}; none found")


def _find_spectrum_hdu(hdul: fits.HDUList, column: str = 'SPECTRUM'):
    for hdu in hdul:
        if isinstance(hdu, fits.BinTableHDU) and column in hdu.columns.names:
            return hdu
    raise ValueError(f"No HDU with a {column} column found")


def _column_units(column: str) -> Tuple[str, str]:
    """(axis label, rms unit suffix) for a spectral column: SPECTRUM is T_A* in K, others (RAW) unitless."""
    return ('T$_A^*$ (K)', ' K') if column == 'SPECTRUM' else (column, '')


def _output_path(output: Path, label: Optional[str]) -> Path:
    return output if label is None else output.with_name(f"{output.stem}_{label}{output.suffix}")


def _plot_on(ax, velocity, group: dict, label: str, vrange, color=None, unit: str = ' K'):
    # Histogram style (CLASS convention): each channel is a flat bin centred on its velocity.
    ax.step(velocity, group['avg_spectrum'], where='mid', linewidth=1.0, color=color,
            label=f"{label} (N={group['count']}, rms={group['rms']:.3f}{unit})")
    ax.axhline(0, color='k', linewidth=0.5, linestyle='--')
    ax.grid(True, alpha=0.3)
    if vrange is not None:
        ax.set_xlim(*vrange)


def _autoscale_y(axes, velocity, groups: List[dict], vrange):
    """Fit the shared y-range to the data inside the plotted velocity range."""
    sel = np.ones_like(velocity, dtype=bool) if vrange is None else \
        (velocity >= min(vrange)) & (velocity <= max(vrange))
    vals = np.concatenate([g['avg_spectrum'][sel] for g in groups])
    vals = vals[np.isfinite(vals)]
    if vals.size:
        lo, hi = vals.min(), vals.max()
        pad = 0.05 * (hi - lo or 1.0)
        for ax in axes:
            ax.set_ylim(lo - pad, hi + pad)


# Hexagonal 7-pixel array as (row, start half-column) in a 3 x 6 half-column grid:
# PX00 in the centre, PX01 top-left, then counterclockwise
# (PX02 left, PX03 bottom-left, PX04 bottom-right, PX05 right, PX06 top-right).
HEX_PIXEL_SLOTS = {1: (0, 1), 6: (0, 3),
                   2: (1, 0), 0: (1, 2), 5: (1, 4),
                   3: (2, 1), 4: (2, 3)}


def _pixel_axes(pixels: List[int], ylabel: str = 'T$_A^*$ (K)') -> Tuple[plt.Figure, Dict[int, plt.Axes]]:
    """
    One axes per pixel, sharing x and y.  Pixels 0-6 are laid out as the
    hexagonal array (HEX_PIXEL_SLOTS); any other set falls back to a 4-column grid.
    """
    hexagonal = set(pixels) <= set(HEX_PIXEL_SLOTS)
    if hexagonal:
        fig = plt.figure(figsize=(15, 9))
        gs = fig.add_gridspec(3, 6)
        slots = {px: gs[HEX_PIXEL_SLOTS[px][0], HEX_PIXEL_SLOTS[px][1]:HEX_PIXEL_SLOTS[px][1] + 2]
                 for px in pixels}
    else:
        ncols = min(4, len(pixels))
        nrows = int(np.ceil(len(pixels) / ncols))
        fig = plt.figure(figsize=(4 * ncols, 3 * nrows))
        gs = fig.add_gridspec(nrows, ncols)
        slots = {px: gs[i // ncols, i % ncols] for i, px in enumerate(pixels)}

    axes, first = {}, None
    for px, spec in slots.items():
        ax = fig.add_subplot(spec, sharex=first, sharey=first)
        first = first or ax
        axes[px] = ax
    # Axis labels on the outer panels: lowest row gets x, leftmost panel of each row gets y.
    rows = {}
    for px, ax in axes.items():
        ss = ax.get_subplotspec()
        rows.setdefault(ss.rowspan.start, []).append((ss.colspan.start, ax))
    bottom = max(rows)
    for row, cells in rows.items():
        min(cells, key=lambda c: c[0])[1].set_ylabel(ylabel)
        if row == bottom:
            for _, ax in cells:
                ax.set_xlabel('Velocity (km/s)')
    return fig, axes


# --group-by fields, in the order they appear in a group key / label.
GROUP_BY_FIELDS = ('scan', 'subscan', 'telescope')


def _group_keys(data, mask: np.ndarray, group_by: List[str]) -> Tuple[List[tuple], List[str]]:
    """
    One key per selected row, built from the requested fields in GROUP_BY_FIELDS order.

    SUBSCAN numbers restart in every scan, so 'subscan' always brings SCAN
    into the key with it.
    """
    fields = [f for f in GROUP_BY_FIELDS if f in group_by or (f == 'scan' and 'subscan' in group_by)]
    columns = {'scan': 'SCAN', 'subscan': 'SUBSCAN', 'telescope': 'TELESCOP'}
    values = []
    for f in fields:
        col = columns[f]
        if col not in data.dtype.names:
            raise ValueError(f"--group-by {f} needs a {col} column; none found")
        v = data[col][mask]
        values.append([_to_str(x) for x in v] if f == 'telescope' else [int(x) for x in v])
    return [tuple(k) for k in zip(*values)], fields


def _group_label(key: tuple, fields: List[str]) -> str:
    parts = dict(zip(fields, key))
    label = []
    if 'scan' in parts:
        label.append(f"{parts['scan']}" + (f".{parts['subscan']}" if 'subscan' in parts else ''))
    if 'telescope' in parts:
        label.append(parts['telescope'])
    return ' '.join(label)


def _flag_report_yaml(flagged: List[tuple], fields: List[str], header_lines: List[str],
                      unit: str = ' K') -> str:
    """
    Turn flagged groups into recommended drop rules, as a commented YAML snippet.

    ``flagged`` holds (key, group, excess, missions) per group, worst first.
    """
    from .drop_rules import report_yaml

    recs: Dict[str, dict] = {}
    for key, g, exc, mids in flagged:
        parts = dict(zip(fields, key))
        why = f"excess {exc:.2f}, N={g['count']}, rms {g['rms']:.3f}{unit}"
        for mid in mids:
            rec = recs.setdefault(mid, {'telescope': {}, 'scans_complete': {},
                                        'scans_telescope': {}, 'notes': []})
            if 'subscan' in parts:
                where = f"scan {parts['scan']} subscan {parts['subscan']}" + \
                        (f" {parts['telescope']}" if 'telescope' in parts else '')
                rec['notes'].append(f"{where}: {why} — no subscan drop rule; drop the whole "
                                    f"scan{' for that telescope' if 'telescope' in parts else ''} if it matters")
            elif 'scan' in parts and 'telescope' in parts:
                rec['scans_telescope'].setdefault(parts['telescope'], {})[parts['scan']] = why
            elif 'scan' in parts:
                rec['scans_complete'][parts['scan']] = why
            else:
                rec['telescope'][parts['telescope']] = why
    return report_yaml(recs, header_lines)


def plot_grouped_spectra(input_fits: Union[str, Path],
                         group_by: List[str],
                         object_filter: Optional[str] = None,
                         weighted: bool = False,
                         velocity_range: Optional[Tuple[float, float]] = None,
                         line_window: Optional[Tuple[float, float]] = None,
                         flag_sigma: float = 5.0,
                         output: Optional[Union[str, Path]] = None,
                         flag_report: Optional[Union[str, Path]] = None,
                         column: str = 'SPECTRUM',
                         echo: Callable[[str], None] = print) -> List[Path]:
    """
    Diagnostic plot: the averaged spectrum of every scan / subscan / pixel group.

    Science rows are grouped on any combination of ``group_by`` fields
    ('scan', 'subscan', 'telescope') and averaged per group.  One figure:

    * left  – heatmap, one row per group, of the group's averaged spectrum
              (colour = T_A*), so a scan or pixel with its own baseline shape
              shows up as a stripe;
    * right – the group's *excess*: rms of its average over the line-free
              channels, divided by the rms expected from radiometric noise
              alone (the per-spectrum rms combined over its N spectra).
              White noise gives ~1; structure that does not average down
              (standing waves, a bad scan) gives more.  Groups whose excess
              lies more than ``flag_sigma`` robust sigmas (1.4826·MAD) above
              the median are flagged in red and listed.

    The excess, not the plain rms, is what is compared: groups hold very
    different numbers of spectra, and a bigger group's average is quieter
    whether or not its data are any good.

    The flagged groups are also turned into recommended ``drop:`` rules for
    the per-mission YAML (see :mod:`oi_zeigt.drop_rules`), one entry per
    MISSION_ID the group's spectra belong to, and printed; with
    ``flag_report`` the same YAML snippet is written to that file.
    Subscan groups have no drop rule and are listed as comments.

    Parameters are as for :func:`plot_averaged_spectra` (including ``column``), plus ``group_by``,
    ``flag_sigma`` and ``flag_report``.  Returns the saved path(s).
    """
    unknown = set(group_by) - set(GROUP_BY_FIELDS)
    if not group_by or unknown:
        raise ValueError(f"--group-by takes any of {', '.join(GROUP_BY_FIELDS)}; got {sorted(unknown) or 'nothing'}")

    with fits.open(input_fits) as hdul:
        hdu = _find_spectrum_hdu(hdul, column)
        data = hdu.data
        velocity = reconstruct_velocity_axis(hdu) / 1000.0  # m/s -> km/s

        mask = select_science_rows(data, object_filter)
        if not mask.any():
            raise ValueError(f"No science spectra left after OBJECT filter '{object_filter}'")
        keys, fields = _group_keys(data, mask, group_by)
        missions = (np.array([_to_str(m) for m in data['MISSION_ID'][mask]])
                    if 'MISSION_ID' in data.dtype.names else np.full(int(mask.sum()), 'UNKNOWN_MISSION'))

        weights = None
        if weighted:
            weights, rms_col = _rms_weights(data[mask])
            bad = ~np.isfinite(weights)
            if bad.any():
                echo(f"Warning: {int(bad.sum())} spectra have no valid {rms_col}; excluded from the weighted average")
                keys = [None if b else k for k, b in zip(keys, bad)]
            echo(f"Weighting each spectrum by 1/{rms_col}²")

        spectra = np.asarray(data[column][mask], dtype=float)

    ylabel, unit = _column_units(column)
    if 'subscan' in group_by and 'scan' not in group_by:
        echo("Note: SUBSCAN numbers restart in every scan, so groups are per scan.subscan")
    echo(f"{input_fits}: {int(mask.sum())} science spectra"
         + (f" matching OBJECT '{object_filter}'" if object_filter else "")
         + f", grouped by {' + '.join(fields)}")

    noise_mask = np.isfinite(velocity)
    if velocity_range is not None:
        noise_mask &= (velocity >= min(velocity_range)) & (velocity <= max(velocity_range))
    if line_window is not None:
        noise_mask &= ~((velocity >= min(line_window)) & (velocity <= max(line_window)))

    # Radiometric noise of each single spectrum, measured on the same line-free channels.
    with np.errstate(invalid='ignore'):
        sigma = np.nanstd(spectra[:, noise_mask], axis=1)

    groups = average_by(spectra, keys, weights, noise_mask)
    keys_arr = np.empty(len(keys), dtype=object)
    keys_arr[:] = keys
    for key, g in groups.items():
        sel = np.array([k == key for k in keys_arr]) & np.isfinite(sigma)
        s = sigma[sel]
        if weights is None:
            expected = np.sqrt(np.sum(s ** 2)) / len(s) if len(s) else np.nan
        else:
            w = weights[sel]
            expected = np.sqrt(np.sum(w ** 2 * s ** 2)) / np.sum(w) if len(s) else np.nan
        g['expected_rms'] = expected
        g['excess'] = g['rms'] / expected if expected > 0 else np.nan

    order = list(groups)  # sorted by key in average_by
    labels = [_group_label(k, fields) for k in order]
    excess = np.array([groups[k]['excess'] for k in order])

    finite = np.isfinite(excess)
    median = float(np.median(excess[finite])) if finite.any() else np.nan
    mad = 1.4826 * float(np.median(np.abs(excess[finite] - median))) if finite.any() else np.nan
    threshold = median + flag_sigma * mad
    flagged = finite & (excess > threshold)

    echo(f"  {len(order)} groups; excess median {median:.2f}, robust sigma {mad:.2f}, "
         f"flag threshold {threshold:.2f} ({flag_sigma:g} sigma)")
    if flagged.any():
        echo(f"  Flagged groups ({int(flagged.sum())}), worst first:")
        for i in sorted(np.flatnonzero(flagged), key=lambda i: -excess[i]):
            g = groups[order[i]]
            echo(f"    {labels[i]:<24} N={g['count']:<6} rms={g['rms']:.4f}{unit}  "
                 f"expected={g['expected_rms']:.4f}{unit}  excess={excess[i]:.2f}")
    else:
        echo("  No group stands out")

    report = _flag_report_yaml(
        [(order[i], groups[order[i]], excess[i],
          sorted(set(missions[np.array([k == order[i] for k in keys_arr])])))
         for i in sorted(np.flatnonzero(flagged), key=lambda i: -excess[i])],
        fields,
        header_lines=[
            f"Recommended drop rules — {'plot_spectra' if column == 'SPECTRUM' else 'plot_raw'} --group-by {','.join(group_by)} --flag-sigma {flag_sigma:g}",
            f"input: {input_fits}" + (f"  (OBJECT '{object_filter}')" if object_filter else ""),
            f"{len(order)} groups; excess = line-free rms / radiometric rms; median {median:.2f}, "
            f"robust sigma {mad:.2f}, flagged above {threshold:.2f}",
            "",
            "Review the plot before adopting any of these. To adopt one, MERGE it into the",
            "mission's existing entry in mission_id_pca_parameters.yml — do not paste a second",
            "key with the same name (YAML silently keeps only the last). For a rule meant for",
            "every mission, put the drop: block under a top-level 'all_missions:' key.",
            "Applied by prepare_for_pca / post_process_data --filter-missions and filter_missions.",
        ],
        unit=unit)
    if flagged.any():
        echo("\n" + report)
    if flag_report is not None:
        flag_report = Path(flag_report)
        flag_report.parent.mkdir(parents=True, exist_ok=True)
        flag_report.write_text(report)
        echo(f"  Flag report saved: {flag_report}")

    # ── Figure ────────────────────────────────────────────────────────────────
    vsel = np.isfinite(velocity)
    if velocity_range is not None:
        vsel &= (velocity >= min(velocity_range)) & (velocity <= max(velocity_range))
    vel = velocity[vsel]
    image = np.array([groups[k]['avg_spectrum'][vsel] for k in order])
    # Velocity may run either way; imshow wants it ascending left to right.
    if vel.size > 1 and vel[0] > vel[-1]:
        vel, image = vel[::-1], image[:, ::-1]
    dv = abs(vel[1] - vel[0]) if vel.size > 1 else 1.0

    finite_img = image[np.isfinite(image)]
    vmax = float(np.percentile(np.abs(finite_img), 99)) if finite_img.size else 1.0

    n = len(order)
    height = min(max(4.0, 0.16 * n + 1.5), 60.0)
    fig, (ax_img, ax_bar) = plt.subplots(1, 2, figsize=(14, height), sharey=True,
                                         gridspec_kw={'width_ratios': [4, 1]})
    im = ax_img.imshow(image, aspect='auto', cmap='RdBu_r', vmin=-vmax, vmax=vmax,
                       interpolation='nearest', origin='upper',
                       extent=(vel[0] - dv / 2, vel[-1] + dv / 2, n - 0.5, -0.5))
    if line_window is not None:
        for v in line_window:
            ax_img.axvline(v, color='k', linestyle='--', linewidth=0.8)
    ax_img.set_yticks(range(n))
    ax_img.set_yticklabels(labels, fontsize=6 if n > 40 else 8)
    for tick, f in zip(ax_img.get_yticklabels(), flagged):
        if f:
            tick.set_color('red')
    ax_img.set_xlabel('Velocity (km/s)')
    ax_img.set_title(f"Averaged {'' if column == 'SPECTRUM' else column + ' '}spectrum per {' + '.join(fields)}"
                     + (f" — {object_filter}" if object_filter else "")
                     + (" (1/RMS² weighted)" if weighted else ""))
    fig.colorbar(im, ax=ax_img, label=ylabel, pad=0.01, fraction=0.03)

    ax_bar.barh(range(n), np.nan_to_num(excess), color=np.where(flagged, 'tab:red', 'steelblue'))
    ax_bar.axvline(1.0, color='k', linewidth=0.5)
    ax_bar.axvline(threshold, color='tab:red', linestyle='--', linewidth=0.8,
                   label=f"flag > {threshold:.2f}")
    ax_bar.set_xlabel('rms / radiometric rms')
    ax_bar.set_title('Excess (line-free)', fontsize=10)
    ax_bar.legend(fontsize=8, loc='lower right')
    ax_bar.grid(True, axis='x', alpha=0.3)
    fig.tight_layout()

    saved = []
    if output is not None:
        output = Path(output)
        output.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output, dpi=150)
        saved.append(output)
        echo(f"  Plot saved: {output}")
        plt.close(fig)
    else:
        plt.show()
    return saved


def plot_averaged_spectra(input_fits: Union[str, Path],
                          mode: str = 'receiver',
                          object_filter: Optional[str] = None,
                          weighted: bool = False,
                          velocity_range: Optional[Tuple[float, float]] = None,
                          line_window: Optional[Tuple[float, float]] = None,
                          output: Optional[Union[str, Path]] = None,
                          column: str = 'SPECTRUM',
                          echo: Callable[[str], None] = print) -> List[Path]:
    """
    Plot averaged spectra from a FITS file, grouped per receiver, polarization or pixel.

    Parameters
    ----------
    input_fits : str or Path
        FITS file with a ``column`` and TELESCOP column.
    mode : {'receiver', 'polarization', 'pixel'}
        'receiver'     – one figure, one averaged spectrum per receiver (all pixels, both pols).
        'polarization' – one figure per receiver, one panel per polarization.
        'pixel'        – one figure per backend (e.g. LFAH, LFAV), one panel per pixel.
    object_filter : str, optional
        Only use science rows whose OBJECT contains this substring.
    weighted : bool
        Weight each spectrum by 1/RMS² (RMS, else RMS_BASELINE column).
        Rows without a valid RMS are excluded from the average.
    velocity_range : (vmin, vmax), optional
        Velocity range in km/s shown on the x-axis (y-axis scales to it).
    line_window : (vmin, vmax), optional
        Line window in km/s left out of the reported rms (and, with
        velocity_range, only line-free channels inside that range are used).
    output : str or Path, optional
        Output image. With several figures, the group label is appended to
        the file stem (out.png -> out_LFAH.png). If None, figures are shown.
    column : str
        Spectral column to average: 'SPECTRUM' (T_A*, K) or e.g. 'RAW' (unitless).
    echo : callable
        Message sink (``click.echo`` from the CLI).

    Returns
    -------
    list of Path
        Saved figure paths (empty when shown on screen).
    """
    if mode not in ('receiver', 'polarization', 'pixel'):
        raise ValueError(f"Unknown mode '{mode}'")

    with fits.open(input_fits) as hdul:
        hdu = _find_spectrum_hdu(hdul, column)
        data = hdu.data
        if 'TELESCOP' not in data.dtype.names:
            raise ValueError("No TELESCOP column: cannot tell receivers/pixels apart")
        velocity = reconstruct_velocity_axis(hdu) / 1000.0  # m/s -> km/s

        mask = select_science_rows(data, object_filter)
        if not mask.any():
            raise ValueError(f"No science spectra left after OBJECT filter '{object_filter}'")

        parsed = [parse_telescop(t) for t in data['TELESCOP'][mask]]
        unparsed = sorted({_to_str(t) for t, p in zip(data['TELESCOP'][mask], parsed) if p is None})
        if unparsed:
            echo(f"Warning: skipping rows with unrecognised TELESCOP {unparsed}")

        weights = None
        if weighted:
            weights, rms_col = _rms_weights(data[mask])
            bad = ~np.isfinite(weights)
            if bad.any():
                echo(f"Warning: {int(bad.sum())} spectra have no valid {rms_col}; excluded from the weighted average")
                parsed = [None if b else p for p, b in zip(parsed, bad)]
            echo(f"Weighting each spectrum by 1/{rms_col}²")

        spectra = np.asarray(data[column][mask], dtype=float)

    echo(f"{input_fits}: {int(mask.sum())} science spectra"
         + (f" matching OBJECT '{object_filter}'" if object_filter else ""))

    noise_mask = np.isfinite(velocity)
    if velocity_range is not None:
        noise_mask &= (velocity >= min(velocity_range)) & (velocity <= max(velocity_range))
    if line_window is not None:
        noise_mask &= ~((velocity >= min(line_window)) & (velocity <= max(line_window)))

    ylabel, unit = _column_units(column)
    what = 'spectrum' if column == 'SPECTRUM' else f"{column} spectrum"
    title_suffix = (f" — {object_filter}" if object_filter else "") + (" (1/RMS² weighted)" if weighted else "")
    figures: List[Tuple[Optional[str], plt.Figure]] = []

    if mode == 'receiver':
        groups = average_by(spectra, [p and (p[0],) for p in parsed], weights, noise_mask)
        fig, ax = plt.subplots(figsize=(12, 5))
        for (rx,), g in groups.items():
            _plot_on(ax, velocity, g, rx, velocity_range, unit=unit)
            echo(f"  {rx}: N={g['count']}, rms={g['rms']:.4f}{unit}")
        _autoscale_y([ax], velocity, list(groups.values()), velocity_range)
        ax.legend(fontsize=9)
        ax.set_xlabel('Velocity (km/s)')
        ax.set_ylabel(ylabel)
        ax.set_title(f"Averaged {what} per receiver" + title_suffix)
        fig.tight_layout()
        figures.append((None, fig))

    elif mode == 'polarization':
        groups = average_by(spectra, [p and (p[0], p[1]) for p in parsed], weights, noise_mask)
        receivers = sorted({rx for rx, _ in groups})
        for rx in receivers:
            pols = sorted(pol for r, pol in groups if r == rx)
            if len(pols) == 1:
                echo(f"Receiver {rx} has only one polarization ({pols[0]}); "
                     f"no more than one polarization exists, plotting it alone")
            fig, axes = plt.subplots(len(pols), 1, figsize=(12, 3.5 * len(pols)),
                                     sharex=True, squeeze=False)
            axes = axes[:, 0]
            for ax, pol in zip(axes, pols):
                g = groups[(rx, pol)]
                _plot_on(ax, velocity, g, f"{rx}{pol}", velocity_range, unit=unit)
                ax.set_ylabel(ylabel)
                ax.legend(fontsize=9, loc='upper right')
                echo(f"  {rx}{pol}: N={g['count']}, rms={g['rms']:.4f}{unit}")
            _autoscale_y(axes, velocity, [groups[(rx, p)] for p in pols], velocity_range)
            axes[-1].set_xlabel('Velocity (km/s)')
            axes[0].set_title(f"{rx}: averaged {what} per polarization" + title_suffix)
            fig.tight_layout()
            figures.append((rx if len(receivers) > 1 else None, fig))

    else:  # pixel
        groups = average_by(spectra, [p and (p[0] + p[1], p[2]) for p in parsed], weights, noise_mask)
        backends = sorted({be for be, _ in groups})
        for be in backends:
            pixels = sorted(px for b, px in groups if b == be)
            fig, axes_by_px = _pixel_axes(pixels, ylabel)
            for px, ax in axes_by_px.items():
                g = groups[(be, px)]
                _plot_on(ax, velocity, g, f"PX{px:02d}", velocity_range, color='steelblue', unit=unit)
                ax.set_title(f"{be}_PX{px:02d}  N={g['count']}  rms={g['rms']:.3f}{unit}", fontsize=9)
                echo(f"  {be}_PX{px:02d}: N={g['count']}, rms={g['rms']:.4f}{unit}")
            _autoscale_y(list(axes_by_px.values()), velocity, [groups[(be, p)] for p in pixels], velocity_range)
            fig.suptitle(f"{be}: averaged {what} per pixel" + title_suffix)
            fig.tight_layout()
            figures.append((be, fig))

    saved = []
    if output is not None:
        output = Path(output)
        output.parent.mkdir(parents=True, exist_ok=True)
        for label, fig in figures:
            path = _output_path(output, label if len(figures) > 1 else None)
            fig.savefig(path, dpi=150)
            saved.append(path)
            echo(f"  Plot saved: {path}")
            plt.close(fig)
    else:
        plt.show()
    return saved


def _parse_date(s: str) -> np.datetime64:
    """'YYYY-MM-DD[...]' -> datetime64[D], zero-padding short fields ('2017-02-1'); NaT if unparseable."""
    try:
        y, m, d = s[:10].split('T')[0].split('-')[:3]
        return np.datetime64(f"{int(y):04d}-{int(m):02d}-{int(d):02d}", 'D')
    except ValueError:
        return np.datetime64('NaT', 'D')


def _row_times(data) -> Optional[np.ndarray]:
    """datetime64 per row from DATE-OBS (date part) + UT (seconds of day); None if unavailable."""
    names = data.dtype.names
    if 'DATE-OBS' not in names or 'UT' not in names:
        return None
    strings = [_to_str(d) for d in data['DATE-OBS']]
    parsed = {s: _parse_date(s) for s in set(strings)}
    dates = np.array([parsed[s] for s in strings], dtype='datetime64[D]')
    if np.isnat(dates).all():
        return None
    ut_ms = np.round(np.asarray(data['UT'], dtype=float) * 1000.0).astype('int64')
    return dates.astype('datetime64[ms]') + ut_ms.astype('timedelta64[ms]')


def window_mean(hdu, rows: np.ndarray, column: str = 'SPECTRUM',
                window: Optional[Tuple[float, float]] = None,
                channel_window: Optional[Tuple[int, int]] = None) -> Tuple[np.ndarray, str]:
    """
    NaN-aware mean of ``column`` over a small window, one value per selected row.

    ``window`` is a velocity range in km/s, evaluated per row from its own
    VELOCITY/DELTAV (rows need not share an axis).  ``channel_window`` is an
    inclusive 0-based channel range.  With neither, the 11 central channels
    are used.  Returns (values, description of the window).
    """
    spectra = np.asarray(hdu.data[column][rows], dtype=float)
    nchan = spectra.shape[1]
    if window is not None:
        vmin, vmax = sorted(window)
        crpix1 = row_crpix1(hdu.data[rows], hdu.header)[:, None]
        ch = np.arange(nchan, dtype=float)[None, :] - (crpix1 - 1.0)
        vel = (np.asarray(hdu.data['VELOCITY'][rows], dtype=float)[:, None]
               + ch * np.asarray(hdu.data['DELTAV'][rows], dtype=float)[:, None]) / 1000.0
        in_win = (vel >= vmin) & (vel <= vmax)
        desc = f"{vmin:g} to {vmax:g} km/s"
    else:
        if channel_window is None:
            channel_window = (nchan // 2 - 5, nchan // 2 + 5)
        cmin, cmax = sorted(channel_window)
        if cmin < 0 or cmax >= nchan:
            raise ValueError(f"Channel window {cmin}-{cmax} outside 0-{nchan - 1}")
        in_win = np.zeros_like(spectra, dtype=bool)
        in_win[:, cmin:cmax + 1] = True
        desc = f"channels {cmin}-{cmax}"
    if not in_win.any():
        raise ValueError(f"Window {desc} contains no channels")
    use = in_win & np.isfinite(spectra)
    n = use.sum(axis=1)
    with np.errstate(invalid='ignore', divide='ignore'):
        values = np.where(n > 0, np.where(use, spectra, 0.0).sum(axis=1) / n, np.nan)
    return values, desc


def _window_series(input_fits: Union[str, Path], object_name: str, column: str,
                   window: Optional[Tuple[float, float]], channel_window: Optional[Tuple[int, int]],
                   x_axis: str) -> Tuple[Dict[str, Tuple[np.ndarray, np.ndarray]], str, bool]:
    """
    Window-averaged intensity of every ``object_name`` spectrum in one file, per TELESCOP.

    Returns ({telescop: (x, values)}, window description, x-is-time).  x is the
    spectrum's sequence number within its TELESCOP (file order), or its
    DATE-OBS + UT time when ``x_axis`` is 'time' and those columns are usable.
    """
    with fits.open(input_fits) as hdul:
        hdu = _find_spectrum_hdu(hdul, column)
        data = hdu.data
        if 'TELESCOP' not in data.dtype.names:
            raise ValueError(f"{input_fits}: no TELESCOP column")
        objects = np.array([_to_str(o).upper() for o in data['OBJECT']])
        rows = np.flatnonzero(objects == object_name.upper())
        if len(rows) == 0:
            raise ValueError(f"{input_fits}: no spectra with OBJECT '{object_name}'")
        telescops = np.array([_to_str(t).upper() for t in data['TELESCOP'][rows]])
        values, desc = window_mean(hdu, rows, column, window, channel_window)
        times = _row_times(data[rows]) if x_axis == 'time' else None

    series = {}
    for tel in sorted(set(telescops)):
        sel = telescops == tel
        x = times[sel] if times is not None else np.arange(int(sel.sum()))
        series[tel] = (x, values[sel])
    return series, desc, times is not None


def plot_window_intensity(input_fits: List[Union[str, Path]],
                          object_name: str = 'S-H_OBS',
                          column: str = 'SPECTRUM',
                          window: Optional[Tuple[float, float]] = None,
                          channel_window: Optional[Tuple[int, int]] = None,
                          x_axis: str = 'index',
                          ncols: int = 5,
                          output: Optional[Union[str, Path]] = None,
                          echo: Callable[[str], None] = print) -> Optional[plt.Figure]:
    """
    Scatter plot of the window-averaged intensity of every ``object_name`` spectrum.

    Each spectrum (OBJECT == ``object_name``, exact match) is reduced to its
    mean over a small window (see :func:`window_mean`).  One figure with one
    panel per receiver (TELESCOP, e.g. LFAH_PX03_S, in name order, ``ncols``
    panels per row); each panel overlays that receiver's spectra from all input files,
    colored per file (legend = file names).  x is ``x_axis``: 'index' (the
    spectrum's number within that receiver, in file order) or 'time'
    (DATE-OBS + UT; falls back to 'index' for a file without them).

    Saved to ``output`` if given (returns None), else the figure is returned
    for the caller to show.
    """
    if x_axis not in ('time', 'index'):
        raise ValueError(f"x_axis must be 'time' or 'index', got '{x_axis}'")

    per_file = []
    for path in input_fits:
        series, desc, is_time = _window_series(path, object_name, column, window, channel_window, x_axis)
        if x_axis == 'time' and not is_time:
            echo(f"  {path}: no usable DATE-OBS/UT; plotting against spectrum index")
        echo(f"{path}: {sum(len(v) for _, v in series.values())} '{object_name}' spectra, "
             f"{len(series)} receivers, window {desc}")
        per_file.append((Path(path).name, series, is_time))
    all_time = x_axis == 'time' and all(t for _, _, t in per_file)

    ylabel, unit = _column_units(column)
    receivers = sorted({tel for _, series, _ in per_file for tel in series})
    ncols = max(1, min(ncols, len(receivers)))
    nrows = -(-len(receivers) // ncols)
    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.5 * ncols, 3.0 * nrows + 1.0),
                             sharex=True, squeeze=False, layout='constrained')
    for ax in axes.flat:
        ax.set_visible(False)

    for k, tel in enumerate(receivers):
        r, c = divmod(k, ncols)
        ax = axes[r, c]
        ax.set_visible(True)
        for i, (name, series, _) in enumerate(per_file):
            if tel not in series:
                echo(f"  {tel}: not in {name}")
                continue
            x, v = series[tel]
            ax.scatter(x, v, s=8, alpha=0.7, linewidths=0, color=colors[i % len(colors)], label=name)
            echo(f"  {tel} {name}: N={len(v)}, median={np.nanmedian(v):.4g}{unit}, "
                 f"min={np.nanmin(v):.4g}, max={np.nanmax(v):.4g}")
        ax.set_title(tel, fontsize=9, loc='left')
        ax.grid(True, alpha=0.3)
        ax.tick_params(labelsize=8)
        if c == 0:
            ax.set_ylabel(ylabel, fontsize=8)
        ax.xaxis.set_tick_params(labelbottom=True)

    # x label under the lowest panel of each column (the last row may be short).
    for c in range(ncols):
        ax = axes[(len(receivers) - 1 - c) // ncols, c]
        ax.set_xlabel('Time (UT)' if all_time else f"{object_name} spectrum index (file order)", fontsize=8)
    if all_time:
        for ax in axes.flat:
            plt.setp(ax.get_xticklabels(), rotation=30, ha='right')
    handles = [plt.Line2D([], [], marker='o', linestyle='', color=colors[i % len(colors)], label=name)
               for i, (name, _, _) in enumerate(per_file)]
    # The title goes in the legend: a suptitle collides with an outside legend.
    fig.legend(handles=handles, loc='outside upper center', ncol=min(len(handles), 2),
               fontsize=9, frameon=False, title=f"{object_name} {column} averaged over {desc}",
               title_fontsize=12)

    if output is None:
        return fig
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=150)
    plt.close(fig)
    echo(f"  Plot saved: {output}")
    return None
