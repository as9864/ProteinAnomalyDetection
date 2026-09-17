#!/usr/bin/env python
"""CLI: run the property/scaffold/AVE bias diagnosis over every target in a
DUD-E or LIT-PCBA download, and write a summary CSV.

Usage:
    python scripts/run_diagnosis.py --dataset dude --root data/raw/dude --out results/dude_bias.csv
    python scripts/run_diagnosis.py --dataset litpcba --root data/raw/litpcba --out results/litpcba_bias.csv
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from proteinanomaly.datasets import (  # noqa: E402
    list_dude_targets,
    list_litpcba_targets,
    load_dude_target,
    load_litpcba_target,
)
from proteinanomaly.pipeline import diagnose_target, summarize  # noqa: E402

LOADERS = {
    "dude": (list_dude_targets, load_dude_target),
    "litpcba": (list_litpcba_targets, load_litpcba_target),
}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset", choices=sorted(LOADERS), required=True)
    parser.add_argument("--root", type=Path, required=True, help="dataset root directory")
    parser.add_argument("--out", type=Path, required=True, help="output summary CSV path")
    parser.add_argument("--targets", nargs="*", default=None, help="subset of target names to run (default: all)")
    parser.add_argument("--ave-repeats", type=int, default=10)
    parser.add_argument("--ave-max-reference", type=int, default=2000)
    args = parser.parse_args()

    list_targets_fn, load_target_fn = LOADERS[args.dataset]
    targets = args.targets or list_targets_fn(args.root)
    if not targets:
        raise SystemExit(f"No targets found under {args.root} — check the directory layout in the module docstring.")

    reports = []
    for target in targets:
        try:
            actives, decoys = load_target_fn(args.root, target)
        except FileNotFoundError as exc:
            print(f"[skip] {target}: {exc}", file=sys.stderr)
            continue
        if len(actives) < 2 or len(decoys) < 2:
            print(f"[skip] {target}: too few parsed molecules (actives={len(actives)}, decoys={len(decoys)})", file=sys.stderr)
            continue
        print(f"[run] {target}: {len(actives)} actives, {len(decoys)} decoys", file=sys.stderr)
        report = diagnose_target(
            target,
            actives,
            decoys,
            ave_n_repeats=args.ave_repeats,
            ave_max_reference=args.ave_max_reference,
        )
        reports.append(report)

    if not reports:
        raise SystemExit("No targets produced a report — nothing to write.")

    summary = summarize(reports)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(args.out)
    print(f"Wrote summary for {len(reports)} targets to {args.out}")

    detail_path = args.out.with_suffix(".detail.json")
    detail = {
        r.target: {
            "property_bias_score": r.property_bias_score,
            "scaffold_bias_score": r.scaffold_bias_score,
            "ave_bias_score": r.ave_bias_score,
            "property_detail": r.property_detail.to_dict(orient="index"),
            "scaffold_detail": r.scaffold_detail,
        }
        for r in reports
    }
    detail_path.write_text(json.dumps(detail, indent=2))
    print(f"Wrote per-target detail to {detail_path}")


if __name__ == "__main__":
    main()
