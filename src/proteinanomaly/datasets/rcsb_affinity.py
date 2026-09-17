"""Free substitute for PDBbind: real protein-ligand complexes with curated
experimental binding affinity, fetched from the public RCSB PDB Data API.

PDBbind itself (http://www.pdbbind.org.cn/) requires a registered account to
download, so it can't be fetched unattended here. RCSB aggregates the same
kind of curated experimental affinity annotations (its own affinity
"validation" report field, ``rcsb_binding_affinity``, cites BindingDB,
Binding MOAD, and PDBbind as sources depending on the entry) via a fully
public REST/GraphQL API with no login, so it serves the same role this
project needs it for (§3.3 of the proposal): a "normal", experimentally
grounded reference distribution of real binder physicochemistry/scaffolds
against which DUD-E/LIT-PCBA actives can be compared.

Two-step fetch, kept as separate functions so the (slow, network-bound)
fetch can be cached to disk and the (fast, local) load can be reused freely:

    df = fetch_binding_affinity_reference(limit=300)
    df.to_csv("data/raw/rcsb_affinity/reference.csv", index=False)
    ...
    ref = load_rcsb_affinity_reference("data/raw/rcsb_affinity/reference.csv")
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import requests

from ..features import MoleculeSet

_SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
_GRAPHQL_URL = "https://data.rcsb.org/graphql"

_UNIT_TO_MOLAR = {"M": 1.0, "mM": 1e-3, "uM": 1e-6, "µM": 1e-6, "nM": 1e-9, "pM": 1e-12, "fM": 1e-15}

_GRAPHQL_QUERY = """
query($ids: [String!]!) {
  entries(entry_ids: $ids) {
    rcsb_id
    rcsb_binding_affinity { comp_id type unit value provenance_code }
    nonpolymer_entities {
      nonpolymer_comp {
        chem_comp { id }
        rcsb_chem_comp_descriptor { SMILES }
      }
    }
  }
}
"""


def _to_neg_log_molar(value: float | None, unit: str | None) -> float | None:
    if value is None or unit not in _UNIT_TO_MOLAR or value <= 0:
        return None
    return -np.log10(value * _UNIT_TO_MOLAR[unit])


def _search_entry_ids(limit: int, timeout: float) -> list[str]:
    query = {
        "query": {
            "type": "terminal",
            "service": "text",
            "parameters": {"attribute": "rcsb_binding_affinity.value", "operator": "exists"},
        },
        "return_type": "entry",
        "request_options": {"paginate": {"start": 0, "rows": limit}},
    }
    resp = requests.post(_SEARCH_URL, json=query, timeout=timeout)
    resp.raise_for_status()
    return [hit["identifier"] for hit in resp.json().get("result_set", [])]


def _fetch_batch(entry_ids: list[str], timeout: float) -> list[dict]:
    resp = requests.post(
        _GRAPHQL_URL, json={"query": _GRAPHQL_QUERY, "variables": {"ids": entry_ids}}, timeout=timeout
    )
    resp.raise_for_status()
    payload = resp.json()
    if "errors" in payload:
        raise RuntimeError(f"RCSB GraphQL error: {payload['errors']}")
    return payload["data"]["entries"]


def fetch_binding_affinity_reference(limit: int = 300, batch_size: int = 50, timeout: float = 30.0) -> pd.DataFrame:
    """Fetch up to ``limit`` PDB entries carrying curated binding-affinity
    annotations, with each affinity record's ligand resolved to a SMILES.

    Network calls: one search request plus ``ceil(limit / batch_size)``
    GraphQL batch requests. Rows with a unit we don't recognize, a
    non-positive value, or a ligand with no resolvable SMILES are dropped.
    """
    entry_ids = _search_entry_ids(limit, timeout=timeout)
    rows = []
    for i in range(0, len(entry_ids), batch_size):
        batch = entry_ids[i : i + batch_size]
        for entry in _fetch_batch(batch, timeout=timeout):
            affinities = entry.get("rcsb_binding_affinity") or []
            if not affinities:
                continue
            comp_to_smiles = {}
            for np_entity in entry.get("nonpolymer_entities") or []:
                comp = (np_entity or {}).get("nonpolymer_comp") or {}
                comp_id = (comp.get("chem_comp") or {}).get("id")
                smiles = (comp.get("rcsb_chem_comp_descriptor") or {}).get("SMILES")
                if comp_id and smiles:
                    comp_to_smiles[comp_id] = smiles

            for aff in affinities:
                comp_id = aff.get("comp_id")
                smiles = comp_to_smiles.get(comp_id)
                neg_log = _to_neg_log_molar(aff.get("value"), aff.get("unit"))
                if smiles is None or neg_log is None:
                    continue
                rows.append(
                    {
                        "pdb_code": entry["rcsb_id"],
                        "comp_id": comp_id,
                        "smiles": smiles,
                        "affinity_type": aff.get("type"),
                        "unit": aff.get("unit"),
                        "value": aff.get("value"),
                        "neg_log_affinity": neg_log,
                        "provenance_code": aff.get("provenance_code"),
                    }
                )
    return pd.DataFrame(rows)


def load_rcsb_affinity_reference(path: Path, dedupe_by_comp_id: bool = True) -> MoleculeSet:
    """Load a CSV produced by ``fetch_binding_affinity_reference`` into a
    MoleculeSet. Dedupes to one row per unique ligand by default, since the
    same ligand often recurs across many PDB entries/targets and we want a
    reference *chemical space*, not a frequency-weighted resample of it.
    """
    df = pd.read_csv(path)
    if dedupe_by_comp_id:
        df = df.drop_duplicates(subset="comp_id")
    return MoleculeSet.from_records(list(zip(df["comp_id"], df["smiles"])))
