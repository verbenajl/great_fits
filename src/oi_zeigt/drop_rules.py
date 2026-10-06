"""
Row drop rules from the per-mission YAML files (``drop:`` blocks).

The same schema is accepted in ``mission_id_parameters.yml`` (telluric / drop
rules) and in ``mission_id_pca_parameters.yml`` (per-mission PCA parameters)::

    all_missions:                 # optional: applies to EVERY mission in the data
      drop:
        telescope: [LFAV_PX03_S]

    2022-04-09_GR_F853:
      drop:
        telescope:                # every row of these telescopes in this mission
          - LFAV_PX03_S
        scans:
          complete:               # every row of these scans, all telescopes
            - 42205
          telescope:              # these scans, only for this telescope
            LFAV_PX01_S:
              - 42199

    2016-05-24_GR_F301:
      drop: flight                # the whole mission

Telescope names match exactly or as a prefix up to an underscore, so ``LFAV_3``
matches ``LFAV_3_S``.  A mission's rules are the union of its own block and the
``all_missions`` block, and the union over every YAML file given.  Rows of every
OBJECT are dropped (science, PCA reference, TSYS, TAU_SIG): tools that link
science rows to calibration rows must rebuild those links afterwards.
"""

from pathlib import Path
from typing import Callable, Dict, Iterable, List, Optional, Union

import numpy as np
import yaml

# Top-level YAML key whose drop block applies to every mission.
ALL_MISSIONS_KEY = 'all_missions'


def _to_str(v):
    return v.decode().strip() if isinstance(v, bytes) else str(v).strip()


def _empty_rules() -> dict:
    return {'flight': False, 'telescope': set(), 'scans_complete': set(), 'scans_telescope': {}}


def _merge(rules: dict, drop_val, source: str) -> None:
    """Merge one YAML ``drop:`` value into ``rules`` (in place)."""
    if not drop_val:
        return
    if drop_val == 'flight':
        rules['flight'] = True
        return
    if not isinstance(drop_val, dict):
        raise ValueError(f"{source}: 'drop' must be 'flight' or a mapping, got {drop_val!r}")
    if drop_val.get('flight'):
        rules['flight'] = True
    rules['telescope'].update(_to_str(t) for t in (drop_val.get('telescope') or []))
    scans = drop_val.get('scans') or {}
    rules['scans_complete'].update(int(s) for s in (scans.get('complete') or []))
    for tele, scan_list in (scans.get('telescope') or {}).items():
        rules['scans_telescope'].setdefault(_to_str(tele), set()).update(int(s) for s in (scan_list or []))


def load_drop_rules(yaml_files: Iterable[Optional[Union[str, Path]]],
                    echo: Callable[[str], None] = print) -> Dict[str, dict]:
    """
    Collect the ``drop:`` blocks of every given YAML file.

    Returns {yaml key: rules}, where the key is a MISSION_ID or ALL_MISSIONS_KEY
    and rules = {'flight', 'telescope', 'scans_complete', 'scans_telescope'}.
    Missing files are skipped with a warning; None entries are ignored.
    """
    by_key: Dict[str, dict] = {}
    for path in yaml_files:
        if not path:
            continue
        path = Path(path)
        if not path.exists():
            echo(f"Warning: drop-rule YAML not found, skipping: {path}")
            continue
        with open(path) as f:
            content = yaml.safe_load(f) or {}
        n = 0
        for key, params in content.items():
            if isinstance(params, dict) and params.get('drop'):
                _merge(by_key.setdefault(str(key), _empty_rules()), params['drop'], f"{path}:{key}")
                n += 1
        echo(f"Drop rules: {n} entr{'y' if n == 1 else 'ies'} in {path}")
    return by_key


def rules_for_mission(by_key: Dict[str, dict], mission_id: str) -> dict:
    """Effective rules for one mission: its own block(s) plus ALL_MISSIONS_KEY."""
    rules = _empty_rules()
    for key, r in by_key.items():
        if key == ALL_MISSIONS_KEY or key == mission_id:
            rules['flight'] |= r['flight']
            rules['telescope'] |= r['telescope']
            rules['scans_complete'] |= r['scans_complete']
            for tele, scans in r['scans_telescope'].items():
                rules['scans_telescope'].setdefault(tele, set()).update(scans)
    return rules


def drop_mask(data, by_key: Dict[str, dict], within: Optional[np.ndarray] = None,
              echo: Callable[[str], None] = print) -> np.ndarray:
    """
    Boolean mask of the rows of ``data`` to DROP under ``by_key``.

    ``data`` needs MISSION_ID; TELESCOP and SCAN are needed only by the rules
    that use them (a rule whose column is missing is reported and skipped).
    With ``within``, only those rows can be dropped (and counted) — pass the
    rows a tool is still keeping so the reported counts are what it removes.
    Every applied rule is reported with its row count; a rule that matches
    nothing is reported too, since that usually means a typo.
    """
    names = data.dtype.names
    drop = np.zeros(len(data), dtype=bool)
    within = np.ones(len(data), dtype=bool) if within is None else np.asarray(within, dtype=bool)
    if not by_key or len(data) == 0:
        return drop
    if 'MISSION_ID' not in names:
        raise ValueError("Drop rules need a MISSION_ID column")

    missions = np.array([_to_str(m) for m in data['MISSION_ID']])
    telescop = np.array([_to_str(t) for t in data['TELESCOP']]) if 'TELESCOP' in names else None
    scan = np.asarray(data['SCAN']).astype(int) if 'SCAN' in names else None

    def tele_match(tele):
        return (telescop == tele) | np.char.startswith(telescop, tele + '_')

    def apply(sel, label):
        sel = sel & within
        n = int(np.sum(sel & ~drop))
        echo(f"  drop {label}: {n} rows" if n else f"  Warning: drop {label}: 0 rows matched")
        drop[:] |= sel

    for mid in sorted(set(missions[within])):
        rules = rules_for_mission(by_key, mid)
        in_mid = missions == mid
        if rules['flight']:
            apply(in_mid, f"{mid} / entire flight")
            continue
        if (rules['telescope'] or rules['scans_telescope']) and telescop is None:
            echo(f"  Warning: no TELESCOP column; telescope rules for {mid} skipped")
        if (rules['scans_complete'] or rules['scans_telescope']) and scan is None:
            echo(f"  Warning: no SCAN column; scan rules for {mid} skipped")
        if telescop is not None:
            for tele in sorted(rules['telescope']):
                apply(in_mid & tele_match(tele), f"{mid} / telescope {tele}")
        if scan is not None:
            for s in sorted(rules['scans_complete']):
                apply(in_mid & (scan == s), f"{mid} / scan {s} (all telescopes)")
            if telescop is not None:
                for tele, scans in sorted(rules['scans_telescope'].items()):
                    for s in sorted(scans):
                        apply(in_mid & tele_match(tele) & (scan == s), f"{mid} / {tele} / scan {s}")
    return drop


def report_yaml(recommendations: Dict[str, dict], header_lines: List[str]) -> str:
    """
    Render recommended drop rules as a YAML snippet with comments.

    ``recommendations`` maps mission -> {'telescope': {tele: comment},
    'scans_complete': {scan: comment}, 'scans_telescope': {tele: {scan: comment}},
    'notes': [comment, ...]}.
    """
    out = [f"# {h}" if h else "#" for h in header_lines]
    if not recommendations:
        out.append("# No groups flagged: nothing to recommend.")
        return "\n".join(out) + "\n"
    for mid in sorted(recommendations):
        rec = recommendations[mid]
        out.append(f"{mid}:")
        if not (rec['telescope'] or rec['scans_complete'] or rec['scans_telescope']):
            out.append("  # (no drop rule expresses these; see notes)")
        else:
            out.append("  drop:")
            if rec['telescope']:
                out.append("    telescope:")
                for tele, c in sorted(rec['telescope'].items()):
                    out.append(f"      - {tele}    # {c}")
            if rec['scans_complete'] or rec['scans_telescope']:
                out.append("    scans:")
                if rec['scans_complete']:
                    out.append("      complete:")
                    for s, c in sorted(rec['scans_complete'].items()):
                        out.append(f"        - {s}    # {c}")
                if rec['scans_telescope']:
                    out.append("      telescope:")
                    for tele, scans in sorted(rec['scans_telescope'].items()):
                        out.append(f"        {tele}:")
                        for s, c in sorted(scans.items()):
                            out.append(f"          - {s}    # {c}")
        for note in rec['notes']:
            out.append(f"  # {note}")
    return "\n".join(out) + "\n"
