"""Shared parsing for whitespace-delimited ``SMILES [id ...]`` files, the
format used by both DUD-E's ``.ism`` files and LIT-PCBA's ``.smi`` files.
"""
from __future__ import annotations

import gzip
from pathlib import Path


def _open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return open(path, "r")


def parse_smi_file(path: Path) -> list[tuple[str, str]]:
    """Return a list of (id, smiles) pairs from a SMILES-per-line file.

    Each non-empty, non-comment line is ``SMILES [id] [... ignored]``. If no
    id column is present, one is synthesized from the line number.
    """
    records: list[tuple[str, str]] = []
    with _open_text(path) as fh:
        for line_no, line in enumerate(fh):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            smiles = parts[0]
            mol_id = parts[1] if len(parts) > 1 else f"{path.stem}_{line_no}"
            records.append((mol_id, smiles))
    return records
