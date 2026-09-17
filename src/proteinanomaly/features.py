"""Molecule-level feature extraction: physicochemical descriptors, Bemis-Murcko
scaffolds, and Morgan fingerprints.

These are the building blocks the three bias metrics (property, scaffold, AVE)
are computed from. Kept dependency-free beyond RDKit/numpy so they can be reused
independently of any particular dataset loader.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Iterable, Optional

import numpy as np
from rdkit import Chem, DataStructs
from rdkit.Chem import Descriptors, rdFingerprintGenerator, rdMolDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold

# Standard property-matching descriptor set used by DUD-E / LIT-PCBA to argue
# actives and decoys are physicochemically indistinguishable. Kept as the
# default set for property_bias metrics.
DESCRIPTOR_NAMES = (
    "MolWt",
    "MolLogP",
    "TPSA",
    "NumHDonors",
    "NumHAcceptors",
    "NumRotatableBonds",
    "NetCharge",
    "NumAromaticRings",
)

_MORGAN_GEN = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)


def smiles_to_mol(smiles: str) -> Optional[Chem.Mol]:
    """Parse a SMILES string, returning None (never raising) on failure."""
    if not smiles:
        return None
    mol = Chem.MolFromSmiles(smiles)
    return mol


def compute_descriptors(mol: Chem.Mol) -> dict[str, float]:
    """Compute the standard property-matching descriptor set for one molecule."""
    return {
        "MolWt": Descriptors.MolWt(mol),
        "MolLogP": Descriptors.MolLogP(mol),
        "TPSA": rdMolDescriptors.CalcTPSA(mol),
        "NumHDonors": float(rdMolDescriptors.CalcNumHBD(mol)),
        "NumHAcceptors": float(rdMolDescriptors.CalcNumHBA(mol)),
        "NumRotatableBonds": float(rdMolDescriptors.CalcNumRotatableBonds(mol)),
        "NetCharge": float(Chem.GetFormalCharge(mol)),
        "NumAromaticRings": float(rdMolDescriptors.CalcNumAromaticRings(mol)),
    }


def murcko_scaffold_smiles(mol: Chem.Mol, generic: bool = False) -> str:
    """Bemis-Murcko scaffold of a molecule, as a canonical SMILES string.

    generic=True collapses atom/bond types (Bemis-Murcko "generic" scaffold),
    which is coarser and useful when the raw-scaffold set is too sparse to
    estimate concentration/overlap statistics reliably.
    """
    scaffold = MurckoScaffold.GetScaffoldForMol(mol)
    if generic:
        scaffold = MurckoScaffold.MakeScaffoldGeneric(scaffold)
    if scaffold is None or scaffold.GetNumAtoms() == 0:
        # Acyclic molecules have an empty Murcko scaffold; treat as its own bucket
        # rather than colliding every acyclic molecule into "".
        return "<acyclic>"
    return Chem.MolToSmiles(scaffold)


def morgan_fingerprint(mol: Chem.Mol):
    """ECFP4-equivalent (radius=2, 2048 bits) Morgan fingerprint."""
    return _MORGAN_GEN.GetFingerprint(mol)


def tanimoto_matrix(fps_a: list, fps_b: list) -> np.ndarray:
    """Dense pairwise Tanimoto similarity matrix, shape (len(fps_a), len(fps_b))."""
    out = np.empty((len(fps_a), len(fps_b)), dtype=np.float64)
    for i, fp in enumerate(fps_a):
        out[i, :] = DataStructs.BulkTanimotoSimilarity(fp, fps_b)
    return out


@dataclass
class MoleculeSet:
    """A labeled collection of molecules (e.g. all actives, or all decoys, for
    one target) with lazily-computed, cached features."""

    ids: list[str]
    smiles: list[str]
    mols: list[Chem.Mol] = field(default_factory=list, repr=False)

    def __post_init__(self):
        if not self.mols:
            parsed_ids, parsed_smiles, parsed_mols = [], [], []
            for i, s in zip(self.ids, self.smiles):
                m = smiles_to_mol(s)
                if m is not None:
                    parsed_ids.append(i)
                    parsed_smiles.append(s)
                    parsed_mols.append(m)
            self.ids, self.smiles, self.mols = parsed_ids, parsed_smiles, parsed_mols

    def __len__(self) -> int:
        return len(self.mols)

    def descriptor_frame(self):
        import pandas as pd

        rows = [compute_descriptors(m) for m in self.mols]
        return pd.DataFrame(rows, index=self.ids)

    def scaffolds(self, generic: bool = False) -> list[str]:
        return [murcko_scaffold_smiles(m, generic=generic) for m in self.mols]

    def fingerprints(self) -> list:
        return [morgan_fingerprint(m) for m in self.mols]

    @classmethod
    def from_records(cls, records: Iterable[tuple[str, str]]) -> "MoleculeSet":
        """Build from an iterable of (id, smiles) pairs."""
        ids, smiles = zip(*records) if records else ((), ())
        return cls(ids=list(ids), smiles=list(smiles))
