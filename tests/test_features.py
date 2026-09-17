from proteinanomaly.features import (
    compute_descriptors,
    morgan_fingerprint,
    murcko_scaffold_smiles,
    smiles_to_mol,
    tanimoto_matrix,
)


def test_smiles_to_mol_valid():
    assert smiles_to_mol("c1ccccc1") is not None


def test_smiles_to_mol_invalid_returns_none():
    assert smiles_to_mol("not a smiles ###") is None
    assert smiles_to_mol("") is None


def test_descriptors_reasonable_for_benzene():
    mol = smiles_to_mol("c1ccccc1")
    desc = compute_descriptors(mol)
    assert desc["NumAromaticRings"] == 1
    assert desc["NumHDonors"] == 0
    assert 75 < desc["MolWt"] < 80


def test_murcko_scaffold_strips_side_chain():
    mol = smiles_to_mol("c1ccccc1CCN")
    assert murcko_scaffold_smiles(mol) == "c1ccccc1"


def test_murcko_scaffold_acyclic_molecule():
    mol = smiles_to_mol("CCCCCC")
    assert murcko_scaffold_smiles(mol) == "<acyclic>"


def test_tanimoto_self_similarity_is_one():
    mol = smiles_to_mol("c1ccccc1CCN")
    fp = morgan_fingerprint(mol)
    sim = tanimoto_matrix([fp], [fp])
    assert sim.shape == (1, 1)
    assert sim[0, 0] == 1.0


def test_moleculeset_from_records_and_features(actives):
    assert len(actives) == 8
    desc_frame = actives.descriptor_frame()
    assert len(desc_frame) == 8
    assert "MolWt" in desc_frame.columns
    scaffolds = actives.scaffolds()
    assert len(scaffolds) == 8
    assert all(s == "c1ccccc1" for s in scaffolds)
    fps = actives.fingerprints()
    assert len(fps) == 8
