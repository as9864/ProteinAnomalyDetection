# Raw data layout

This directory is git-ignored (`data/raw/*`) — benchmark downloads are large
and/or license-gated, so fetch them yourself and place them here.

## DUD-E

Download per-target archives from https://dude.docking.org/ and lay them out as:

```
data/raw/dude/<target>/actives_final.ism
data/raw/dude/<target>/decoys_final.ism
```

## LIT-PCBA

Download from https://drugdesign.unistra.fr/LIT-PCBA/ and lay out as:

```
data/raw/litpcba/<target>/actives.smi
data/raw/litpcba/<target>/inactives.smi
```

LIT-PCBA ships both a "full" and an "AVE-unbiased" split per target — pick
one consistently; comparing this tool's `ave_bias` score across both is
itself a useful sanity check (the unbiased split should score near 0).

## PDBbind

Requires a (free) *registered account* at http://www.pdbbind.org.cn/ —
confirmed in practice that the download is not reachable unattended without
logging in. Used here as the "normal" reference distribution — a real,
experimentally resolved protein-ligand complex with a measured affinity, as
opposed to a computationally chosen decoy. Lay out as:

```
data/raw/pdbbind/index/INDEX_general_PL_data.<year>
data/raw/pdbbind/<pdbid>/<pdbid>_ligand.sdf
```

### Substitute actually used: RCSB binding-affinity API

Until/unless a registered PDBbind download is available, `datasets/rcsb_affinity.py`
fetches the same kind of data (real complexes with curated experimental
Ki/Kd/IC50, sourced by RCSB from BindingDB / Binding MOAD / PDBbind) from the
public, login-free RCSB PDB REST/GraphQL API:

```bash
python scripts/fetch_rcsb_binding_reference.py --limit 300 --out data/raw/rcsb_affinity/reference.csv
```

Both loaders return the same `MoleculeSet` type, so swapping in a real
PDBbind download later is a one-line change in any script that consumes it.
