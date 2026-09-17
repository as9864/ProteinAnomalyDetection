#!/usr/bin/env python
"""Compare a benchmark target's *actives* against the real-binder reference
set (fetch_rcsb_binding_reference.py output) using the same property/scaffold
metrics used for actives-vs-decoys diagnosis.

This answers a different question than run_diagnosis.py: not "are this
target's actives and decoys separable," but "does this target's active set
itself look like the general population of real, experimentally confirmed
binders, or is it a skewed sample?" (proposal.md, section 5.1).

Usage:
    python scripts/compare_actives_to_reference.py --dataset dude --root data/raw/dude --target ampc \
        --reference data/raw/rcsb_affinity/reference.csv
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from proteinanomaly.datasets import load_dude_target, load_litpcba_target  # noqa: E402
from proteinanomaly.datasets.rcsb_affinity import load_rcsb_affinity_reference  # noqa: E402
from proteinanomaly.metrics import property_bias, scaffold_bias  # noqa: E402

LOADERS = {"dude": load_dude_target, "litpcba": load_litpcba_target}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset", choices=sorted(LOADERS), required=True)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--target", required=True)
    parser.add_argument("--reference", type=Path, required=True)
    args = parser.parse_args()

    actives, _decoys = LOADERS[args.dataset](args.root, args.target)
    reference = load_rcsb_affinity_reference(args.reference)

    prop = property_bias(actives.descriptor_frame(), reference.descriptor_frame())
    scaf = scaffold_bias(actives.scaffolds(), reference.scaffolds())

    print(f"target: {args.target} ({args.dataset})")
    print(f"  n_actives={len(actives)}  n_reference={len(reference)}")
    print(f"  property_bias vs reference  = {prop.score:.4f}")
    print(f"  scaffold_bias vs reference  = {scaf.score:.4f}  (active_scaffold_coverage={scaf.active_scaffold_coverage:.4f})")


if __name__ == "__main__":
    main()
