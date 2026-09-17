# ProteinAnomalyDetection

Diagnoses **artificial enrichment** in virtual-screening benchmarks (DUD-E,
LIT-PCBA) — the well-documented failure mode where decoys differ from actives
in trivial, non-binding-related ways (physicochemistry, scaffold identity,
overall chemical-space location), so a classifier can score well without
learning anything about actual binding.

Three complementary metrics, computed per target:

| Metric | Question it answers | Module |
|---|---|---|
| Property bias | Are actives/decoys separable on MW, logP, TPSA, HBD/HBA, rotatable bonds, charge, aromatic ring count alone? | [`metrics/property_bias.py`](src/proteinanomaly/metrics/property_bias.py) |
| Scaffold bias | Do actives/decoys share a Bemis-Murcko scaffold vocabulary, or are they built from disjoint chemical series? | [`metrics/scaffold_bias.py`](src/proteinanomaly/metrics/scaffold_bias.py) |
| AVE bias | Would a nearest-neighbor memorizer (no real model) score above chance under random train/test splits? (Wallach & Heifets, 2018) | [`metrics/ave_bias.py`](src/proteinanomaly/metrics/ave_bias.py) |

All three are in `[0, 1]` except AVE bias, which is unbounded but centered at
0 (0 = not trivially memorizable; positive = memorizable; in principle it can
go slightly negative if decoys are, by chance, *closer* to the wrong class).

PDBbind (real, experimentally resolved complexes with measured affinity)
serves as the "normal" reference distribution against which DUD-E/LIT-PCBA
actives/decoys can be compared. PDBbind itself requires a registered account
to download, so [`datasets/rcsb_affinity.py`](src/proteinanomaly/datasets/rcsb_affinity.py)
fetches the same kind of curated-affinity real complexes from the public,
login-free RCSB PDB API as a drop-in substitute (`fetch_binding_affinity_reference()`
returns the same `MoleculeSet`-compatible data `datasets/pdbbind.py` would).
See [`data/README.md`](data/README.md) for dataset download/layout
instructions — none of the datasets are bundled here (large, and PDBbind is
license-gated).

**Validated against real data**: [`doc/제안서.md`](doc/제안서.md) §5.3 has
results from an actual run — DUD-E's AMPC target vs LIT-PCBA's TP53
target, plus both compared against the RCSB affinity reference set. All three
metrics agree LIT-PCBA is less biased than DUD-E on this target pair, matching
the literature's claim about LIT-PCBA's design goal.

## Setup

```bash
python -m venv .venv
source .venv/Scripts/activate   # Windows Git Bash; use .venv\Scripts\Activate.ps1 in PowerShell
pip install -e .
```

## Usage

```bash
python scripts/run_diagnosis.py --dataset dude --root data/raw/dude --out results/dude_bias.csv
python scripts/run_diagnosis.py --dataset litpcba --root data/raw/litpcba --out results/litpcba_bias.csv
```

Writes a per-target summary CSV (`target, n_actives, n_decoys, property_bias,
scaffold_bias, ave_bias`) plus a `*.detail.json` with the full per-descriptor
and per-scaffold breakdown for each target.

To compare a target's actives against the RCSB (PDBbind-substitute) real-binder
reference set instead of its own decoys:

```bash
python scripts/fetch_rcsb_binding_reference.py --limit 300 --out data/raw/rcsb_affinity/reference.csv
python scripts/compare_actives_to_reference.py --dataset dude --root data/raw/dude --target ampc \
    --reference data/raw/rcsb_affinity/reference.csv
```

Programmatic use:

```python
from proteinanomaly.datasets import load_dude_target
from proteinanomaly.pipeline import diagnose_target

actives, decoys = load_dude_target("data/raw/dude", "egfr")
report = diagnose_target("egfr", actives, decoys)
print(report.to_summary_row())
```

## Tests

```bash
pytest
```

Tests run against small synthetic molecule sets (no benchmark download
needed) with a deliberately "matched" case and a deliberately "biased" case
for each metric, plus a class-size-imbalance regression case for the
scaffold metric, so a metric that's inverted, degenerate, or ignoring its
input gets caught immediately.
