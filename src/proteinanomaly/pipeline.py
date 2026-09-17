"""Orchestrates the three bias metrics over one or more (actives, decoys)
target pairs and assembles a summary report.
"""
from __future__ import annotations

from dataclasses import dataclass

import pandas as pd

from .features import MoleculeSet
from .metrics import ave_bias, property_bias, scaffold_bias


@dataclass
class TargetBiasReport:
    target: str
    n_actives: int
    n_decoys: int
    property_bias_score: float
    scaffold_bias_score: float
    ave_bias_score: float
    property_detail: pd.DataFrame
    scaffold_detail: dict
    ave_detail: pd.DataFrame

    def to_summary_row(self) -> dict:
        return {
            "target": self.target,
            "n_actives": self.n_actives,
            "n_decoys": self.n_decoys,
            "property_bias": self.property_bias_score,
            "scaffold_bias": self.scaffold_bias_score,
            "ave_bias": self.ave_bias_score,
        }


def diagnose_target(
    target: str,
    actives: MoleculeSet,
    decoys: MoleculeSet,
    ave_n_repeats: int = 10,
    ave_max_reference: int = 2000,
    scaffold_n_repeats: int = 200,
    random_state: int = 0,
) -> TargetBiasReport:
    """Run all three bias metrics for one target's actives/decoys pair."""
    actives_desc = actives.descriptor_frame()
    decoys_desc = decoys.descriptor_frame()
    prop_result = property_bias(actives_desc, decoys_desc)

    scaf_result = scaffold_bias(
        actives.scaffolds(), decoys.scaffolds(), n_repeats=scaffold_n_repeats, random_state=random_state
    )

    ave_result = ave_bias(
        actives.fingerprints(),
        decoys.fingerprints(),
        n_repeats=ave_n_repeats,
        max_reference=ave_max_reference,
        random_state=random_state,
    )

    return TargetBiasReport(
        target=target,
        n_actives=len(actives),
        n_decoys=len(decoys),
        property_bias_score=prop_result.score,
        scaffold_bias_score=scaf_result.score,
        ave_bias_score=ave_result.score,
        property_detail=prop_result.per_descriptor,
        scaffold_detail=scaf_result.to_dict(),
        ave_detail=ave_result.per_repeat,
    )


def summarize(reports: list[TargetBiasReport]) -> pd.DataFrame:
    """Collapse a list of per-target reports into one summary DataFrame."""
    return pd.DataFrame([r.to_summary_row() for r in reports]).set_index("target")
