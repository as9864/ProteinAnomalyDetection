"""Loader for LIT-PCBA.

Expected layout (matching the LIT-PCBA download as distributed):

    data/raw/litpcba/<target>/actives.smi
    data/raw/litpcba/<target>/inactives.smi

Each ``.smi`` file is whitespace-delimited ``SMILES id ...`` per line.
LIT-PCBA was explicitly constructed to reduce AVE/property/scaffold bias
relative to DUD-E, so it is a useful low-bias reference point when
interpreting scores from this tool.
"""
from __future__ import annotations

from pathlib import Path

from ..features import MoleculeSet
from ._smi_utils import parse_smi_file


def list_litpcba_targets(root: Path) -> list[str]:
    root = Path(root)
    if not root.exists():
        return []
    return sorted(p.name for p in root.iterdir() if p.is_dir())


def _find(target_dir: Path, stem: str) -> Path | None:
    for suffix in (".smi", ".smi.gz"):
        candidate = target_dir / f"{stem}{suffix}"
        if candidate.exists():
            return candidate
    return None


def load_litpcba_target(root: Path, target: str) -> tuple[MoleculeSet, MoleculeSet]:
    """Load (actives, inactives) MoleculeSets for one LIT-PCBA target."""
    target_dir = Path(root) / target
    actives_path = _find(target_dir, "actives")
    inactives_path = _find(target_dir, "inactives")
    if actives_path is None or inactives_path is None:
        raise FileNotFoundError(
            f"Expected actives/inactives .smi files under {target_dir}"
        )
    actives = MoleculeSet.from_records(parse_smi_file(actives_path))
    decoys = MoleculeSet.from_records(parse_smi_file(inactives_path))
    return actives, decoys
