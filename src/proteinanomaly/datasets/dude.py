"""Loader for DUD-E (Directory of Useful Decoys, Enhanced).

Expected layout (matching the DUD-E download as distributed):

    data/raw/dude/<target>/actives_final.ism
    data/raw/dude/<target>/decoys_final.ism

Each ``.ism`` file is whitespace-delimited ``SMILES id ...`` per line.
Gzip-compressed variants (``.ism.gz``) are also accepted.
"""
from __future__ import annotations

from pathlib import Path

from ..features import MoleculeSet
from ._smi_utils import parse_smi_file


def list_dude_targets(root: Path) -> list[str]:
    root = Path(root)
    if not root.exists():
        return []
    return sorted(p.name for p in root.iterdir() if p.is_dir())


def _find(target_dir: Path, stem: str) -> Path | None:
    for suffix in (".ism", ".ism.gz", ".smi", ".smi.gz"):
        candidate = target_dir / f"{stem}{suffix}"
        if candidate.exists():
            return candidate
    return None


def load_dude_target(root: Path, target: str) -> tuple[MoleculeSet, MoleculeSet]:
    """Load (actives, decoys) MoleculeSets for one DUD-E target."""
    target_dir = Path(root) / target
    actives_path = _find(target_dir, "actives_final")
    decoys_path = _find(target_dir, "decoys_final")
    if actives_path is None or decoys_path is None:
        raise FileNotFoundError(
            f"Expected actives_final/decoys_final .ism files under {target_dir}"
        )
    actives = MoleculeSet.from_records(parse_smi_file(actives_path))
    decoys = MoleculeSet.from_records(parse_smi_file(decoys_path))
    return actives, decoys
