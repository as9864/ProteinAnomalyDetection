#!/usr/bin/env python
"""Fetch the free RCSB-based substitute for a PDBbind reference set (see
proteinanomaly.datasets.rcsb_affinity for why) and cache it to CSV.

Usage:
    python scripts/fetch_rcsb_binding_reference.py --limit 300 --out data/raw/rcsb_affinity/reference.csv
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from proteinanomaly.datasets.rcsb_affinity import fetch_binding_affinity_reference  # noqa: E402


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--limit", type=int, default=300)
    parser.add_argument("--batch-size", type=int, default=50)
    parser.add_argument("--out", type=Path, default=Path("data/raw/rcsb_affinity/reference.csv"))
    args = parser.parse_args()

    print(f"Fetching up to {args.limit} RCSB entries with curated binding affinity...", file=sys.stderr)
    df = fetch_binding_affinity_reference(limit=args.limit, batch_size=args.batch_size)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.out, index=False)
    print(f"Wrote {len(df)} affinity records ({df['comp_id'].nunique()} unique ligands) to {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
