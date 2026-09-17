"""Shared synthetic fixtures for metric tests.

Two scenarios per test:
  - "matched": actives/decoys share scaffold family and physchem range
    (benzene-ring + 2-carbon-linker compounds differing only in the
    terminal substituent and the linker's terminal heteroatom).
  - "biased": decoys are long acyclic alkanes — a different scaffold
    ("<acyclic>" vs. the benzene ring), different physchem (no TPSA/HBD/HBA,
    much higher logP), and far away in fingerprint space.

These aren't meant to resemble real DUD-E/LIT-PCBA data quantitatively; they
exist to give each metric a case that should score low and a case that
should score unambiguously high, so a metric that's inverted, degenerate,
or ignores its input gets caught.
"""
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from proteinanomaly.features import MoleculeSet  # noqa: E402

_ACTIVES_SMILES = [
    "c1ccccc1CCN",
    "c1ccc(Cl)cc1CCN",
    "c1ccc(F)cc1CCN",
    "c1ccc(C)cc1CCN",
    "c1ccc(O)cc1CCN",
    "c1ccc(N)cc1CCN",
    "c1ccc(Br)cc1CCN",
    "c1ccc(cc1)CCN",
]

_DECOYS_MATCHED_SMILES = [
    "c1ccccc1CCO",
    "c1ccc(Cl)cc1CCO",
    "c1ccc(F)cc1CCO",
    "c1ccc(C)cc1CCO",
    "c1ccc(O)cc1CCO",
    "c1ccc(N)cc1CCO",
    "c1ccc(Br)cc1CCO",
    "c1ccc(cc1)CCO",
]

_DECOYS_BIASED_SMILES = [
    "CCCCCCCCCC",
    "CCCCCCCCCCC",
    "CCCCCCCCCCCC",
    "CCCCCCCCCCCCC",
    "CCCCCCCCCCCCCC",
    "CCCCCCCCCCCCCCC",
    "CCCCCCCCCCCCCCCC",
    "CCCCCCCCCCCCCCCCC",
]


def _to_records(smiles_list):
    return [(f"m{i}", s) for i, s in enumerate(smiles_list)]


@pytest.fixture
def actives():
    return MoleculeSet.from_records(_to_records(_ACTIVES_SMILES))


@pytest.fixture
def decoys_matched():
    return MoleculeSet.from_records(_to_records(_DECOYS_MATCHED_SMILES))


@pytest.fixture
def decoys_biased():
    return MoleculeSet.from_records(_to_records(_DECOYS_BIASED_SMILES))
