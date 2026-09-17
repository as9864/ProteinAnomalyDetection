"""Loader for PDBbind, used here as the "normal" reference distribution
(real, experimentally-resolved protein-ligand complexes with measured
affinities) that DUD-E/LIT-PCBA actives/decoys are compared against.

Expected layout:

    data/raw/pdbbind/index/INDEX_general_PL_data.<year>
    data/raw/pdbbind/<pdbid>/<pdbid>_ligand.sdf   (or .mol2)

The index file format (unchanged across PDBbind releases) is whitespace
columns with a ``//`` separator before a trailing ``(ligand_name)``:

    3zzf  2.20  2012  2.72  Ki=1.9mM      // 3zzf.pdf (MLY)

Columns: pdb_code, resolution, release_year, neg_log_affinity, affinity_str,
[ligand_name].
"""
from __future__ import annotations

import re
from pathlib import Path

import pandas as pd
from rdkit import Chem

_LIGAND_NAME_RE = re.compile(r"\(([^)]+)\)\s*$")


def load_pdbbind_index(index_path: Path) -> pd.DataFrame:
    """Parse a PDBbind ``INDEX_general_PL_data.*`` file into a DataFrame with
    columns: pdb_code, resolution, release_year, neg_log_affinity,
    affinity_str, ligand_name.
    """
    rows = []
    with open(index_path, "r") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            ligand_match = _LIGAND_NAME_RE.search(line)
            ligand_name = ligand_match.group(1) if ligand_match else None
            head = line.split("//")[0].split()
            if len(head) < 5:
                continue
            pdb_code, resolution, release_year, neg_log_affinity, affinity_str = head[:5]
            try:
                resolution_val = float(resolution)
            except ValueError:
                resolution_val = float("nan")  # e.g. "NMR"
            try:
                rows.append(
                    {
                        "pdb_code": pdb_code,
                        "resolution": resolution_val,
                        "release_year": int(release_year),
                        "neg_log_affinity": float(neg_log_affinity),
                        "affinity_str": affinity_str,
                        "ligand_name": ligand_name,
                    }
                )
            except ValueError:
                continue
    return pd.DataFrame(rows)


def load_pdbbind_ligand_smiles(root: Path, pdb_code: str) -> str | None:
    """Best-effort load of a ligand's SMILES from the per-complex directory.

    Tries ``<pdb_code>_ligand.sdf`` first, then ``.mol2``. Returns None if
    neither is present or RDKit fails to parse the file (e.g. because it
    lacks bond-order/valence info, which is common for crystallographic
    mol2/sdf ligand records).
    """
    complex_dir = Path(root) / pdb_code
    sdf_path = complex_dir / f"{pdb_code}_ligand.sdf"
    if sdf_path.exists():
        supplier = Chem.SDMolSupplier(str(sdf_path))
        for mol in supplier:
            if mol is not None:
                return Chem.MolToSmiles(mol)

    mol2_path = complex_dir / f"{pdb_code}_ligand.mol2"
    if mol2_path.exists():
        mol = Chem.MolFromMol2File(str(mol2_path))
        if mol is not None:
            return Chem.MolToSmiles(mol)

    return None
