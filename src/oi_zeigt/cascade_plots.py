"""
Cascade (waterfall) plots of spectra, grouped by MISSION_ID/SCAN/TELESCOP.
"""

from pathlib import Path
from typing import List, Optional, Tuple, Union

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

from .basic_io import reconstruct_velocity_axis


def _to_str(v):
    return v.decode().strip() if isinstance(v, bytes) else str(v).strip()


def _draw_waterfall(ax, group_spectra: np.ndarray, velocity_axis: np.ndarray,
                    cmap_name: str = 'viridis'):
    """
    Draw a waterfall-style image: one row per spectrum, color = amplitude.

    Matches the color-coding convention used by pca_correct's diagnostic
    plots (imshow + colorbar) — unlike an offset line stack, rows never
    overlap regardless of how many spectra are in the group.

    Returns the image artist (`AxesImage`), for the caller to attach a colorbar to.
    """
    extent = [velocity_axis[0], velocity_axis[-1], len(group_spectra), 0]
    im = ax.imshow(group_spectra, aspect='auto', extent=extent,
                   interpolation='nearest', cmap=cmap_name)
    return im


def draw_cascade_plots(input_fits: Union[str, Path],
                       output_folder: Optional[Union[str, Path]] = None,
                       hdu_index: Optional[int] = None,
                       max_figsize_height: float = 30.0,
                       panel_cols: int = 4,
                       object_filter: Optional[str] = None,
                       velocity_range: Optional[Tuple[float, float]] = None) -> List[Path]:
    """
    Draw cascade (waterfall) plots of spectra, one per MISSION_ID/SCAN/TELESCOP
    group found in the input FITS file.

    Each group is rendered as an imshow image — one row per spectrum, color
    encoding amplitude, with a colorbar — the same color-coding convention
    used throughout pca_correct's diagnostic plots. This scales cleanly to
    large groups (no line overlap, unlike a stacked-line cascade).

    Parameters
    ----------
    input_fits : str or Path
        Path to the FITS file containing spectra.
    output_folder : str or Path, optional
        Directory to save the plots into. One subdirectory per MISSION_ID is
        created inside it, matching the pca_correct diagnostic-plot convention.
        If None (default), nothing is written to disk — instead, all groups
        are drawn as panels in a single interactive figure and shown on screen.
    hdu_index : int, optional
        HDU index to read. If None, the first HDU with a SPECTRUM column is used.
    max_figsize_height : float, optional
        Cap on the figure height in inches for the saved-to-file mode,
        regardless of how many spectra are in a group (default: 30).
    panel_cols : int, optional
        Number of columns in the panel grid when showing on screen (default: 4).
    object_filter : str, optional
        Only include rows whose OBJECT contains this substring (case-insensitive).
        If None (default), all rows are included regardless of OBJECT.
    velocity_range : (vmin, vmax), optional
        Only draw the channels inside this velocity range (km/s). The colour
        scale then follows the data inside the range. Ignored (with a warning)
        when the file has no velocity axis.

    Returns
    -------
    list of Path
        Paths to the saved plot files, one per group. Empty when
        `output_folder` is None (screen mode writes nothing to disk).
    """
    input_fits = Path(input_fits)

    with fits.open(input_fits) as hdul:
        hdu = None
        if hdu_index is not None:
            hdu = hdul[hdu_index]
        else:
            for h in hdul:
                if hasattr(h, 'data') and h.data is not None and h.data.dtype.names is not None \
                        and 'SPECTRUM' in h.data.dtype.names:
                    hdu = h
                    break
        if hdu is None:
            raise ValueError("No HDU with a SPECTRUM column found")

        data = hdu.data
        n = len(data)
        if n == 0:
            raise ValueError("No spectra found in the input file")

        for col in ('MISSION_ID', 'SCAN', 'TELESCOP'):
            if col not in data.dtype.names:
                raise ValueError(f"{col} column not found in FITS data — required for cascade_plots")

        spectra    = np.asarray(data['SPECTRUM'], dtype=np.float64)
        missions   = np.array([_to_str(v) for v in data['MISSION_ID']])
        scans      = np.array(data['SCAN'])
        telescopes = np.array([_to_str(v) for v in data['TELESCOP']])

        try:
            velocity_axis = reconstruct_velocity_axis(hdu) / 1000.0  # m/s -> km/s
            x_label = 'Velocity (km/s)'
        except (KeyError, IndexError):
            velocity_axis = np.arange(spectra.shape[1])
            x_label = 'Channel'
            if velocity_range is not None:
                print("Warning: no velocity axis in this file; --velocity-range ignored")
                velocity_range = None

        if velocity_range is not None:
            vsel = (velocity_axis >= min(velocity_range)) & (velocity_axis <= max(velocity_range))
            if not np.any(vsel):
                raise ValueError(f"No channels inside velocity range {velocity_range} km/s "
                                 f"(file covers {velocity_axis.min():.1f} to {velocity_axis.max():.1f})")
            spectra, velocity_axis = spectra[:, vsel], velocity_axis[vsel]

        if object_filter:
            if 'OBJECT' not in data.dtype.names:
                raise ValueError("OBJECT column not found in FITS data — cannot use --object")
            obj_strs = np.array([_to_str(v).lower() for v in data['OBJECT']])
            obj_mask = np.char.find(obj_strs, object_filter.lower()) >= 0
            n_before = n
            spectra, missions, scans, telescopes = (
                spectra[obj_mask], missions[obj_mask], scans[obj_mask], telescopes[obj_mask])
            n = int(np.sum(obj_mask))
            print(f"Object filter '{object_filter}': {n} of {n_before} spectra match")
            if n == 0:
                raise ValueError(f"No spectra match OBJECT containing '{object_filter}'")

    groups = sorted(set(zip(missions, scans, telescopes)))
    print(f"Found {len(groups)} (MISSION_ID, SCAN, TELESCOP) groups across {n} spectra")

    def _group_spectra(mission_id, scan, telescope):
        group_mask = (missions == mission_id) & (scans == scan) & (telescopes == telescope)
        return spectra[group_mask]

    if output_folder is not None:
        output_folder = Path(output_folder)
        output_folder.mkdir(parents=True, exist_ok=True)
        saved_paths = []

        for mission_id, scan, telescope in groups:
            group_spectra = _group_spectra(mission_id, scan, telescope)
            n_spec = len(group_spectra)
            if n_spec == 0:
                continue

            fig_height = float(np.clip(0.1 * n_spec, 4.0, max_figsize_height))
            fig, ax = plt.subplots(figsize=(14, fig_height))

            im = _draw_waterfall(ax, group_spectra, velocity_axis)

            ax.set_xlabel(x_label, fontsize=13)
            ax.set_ylabel('Spectrum #', fontsize=13)
            ax.set_title(f'{mission_id} | Scan {scan} | Telescope {telescope} | {n_spec} spectra',
                        fontsize=15, fontweight='bold')
            ax.tick_params(labelsize=11)
            cbar = fig.colorbar(im, ax=ax, pad=0.02)
            cbar.set_label('T$_A^*$ (K)', fontsize=12)
            cbar.ax.tick_params(labelsize=10)

            mission_subdir = output_folder / str(mission_id)
            mission_subdir.mkdir(parents=True, exist_ok=True)
            safe_telescope = str(telescope).replace('/', '_').replace(' ', '_')
            plot_path = mission_subdir / f'cascade_{mission_id}_scan{scan}_{safe_telescope}.png'

            fig.tight_layout()
            fig.savefig(plot_path, dpi=150, bbox_inches='tight')
            plt.close(fig)

            saved_paths.append(plot_path)
            print(f"  Saved: {plot_path}  ({n_spec} spectra)")

        return saved_paths

    # ---- No output_folder: show all groups as panels in one figure ----
    n_groups = len(groups)
    if n_groups == 0:
        raise ValueError("No (MISSION_ID, SCAN, TELESCOP) groups found")

    n_cols = min(panel_cols, n_groups)
    n_rows = int(np.ceil(n_groups / n_cols))
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(5 * n_cols, 3 * n_rows), squeeze=False)

    for idx, (mission_id, scan, telescope) in enumerate(groups):
        ax = axes[idx // n_cols][idx % n_cols]
        group_spectra = _group_spectra(mission_id, scan, telescope)
        n_spec = len(group_spectra)
        if n_spec == 0:
            ax.set_visible(False)
            continue

        im = _draw_waterfall(ax, group_spectra, velocity_axis)

        ax.set_title(f'{mission_id} | Scan {scan} | {telescope}\n({n_spec} spectra)', fontsize=8)
        ax.set_xlabel(x_label, fontsize=7)
        ax.set_ylabel('Spectrum #', fontsize=7)
        ax.tick_params(labelsize=6)
        cbar = fig.colorbar(im, ax=ax, pad=0.02)
        cbar.ax.tick_params(labelsize=6)

    # Hide unused axes in the last row
    for idx in range(n_groups, n_rows * n_cols):
        axes[idx // n_cols][idx % n_cols].set_visible(False)

    fig.suptitle(f'Cascade plots — {n_groups} groups', fontsize=12, fontweight='bold')
    fig.tight_layout()
    plt.show()

    return []
