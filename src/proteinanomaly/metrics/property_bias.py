"""Metric 1: physicochemical property-matching bias.

DUD-E/LIT-PCBA claim decoys are property-matched to actives (similar MW, logP,
HBD/HBA, rotatable bonds, ...) so a classifier can't trivially separate the two
classes using those properties alone. This module quantifies how well that
claim actually holds for a given actives/decoys pair, per descriptor and
aggregated.

Per descriptor we report:
  - KS statistic (0 = identical distributions, 1 = fully separated) — scale-free,
    used as the primary bias signal.
  - Wasserstein distance in standardized units (divided by the pooled std), a
    signed-magnitude-free effect-size analogue — useful for ranking which
    descriptor drives the bias.

The aggregate score is the mean KS statistic across descriptors, in [0, 1].
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy.stats import ks_2samp, wasserstein_distance

from ..features import DESCRIPTOR_NAMES


@dataclass
class PropertyBiasResult:
    per_descriptor: pd.DataFrame  # index=descriptor, columns=[ks_stat, ks_pvalue, standardized_wasserstein]
    score: float  # mean KS statistic across descriptors, in [0, 1]

    def to_dict(self) -> dict:
        return {
            "score": self.score,
            "per_descriptor": self.per_descriptor.to_dict(orient="index"),
        }


def property_bias(
    actives_desc: pd.DataFrame,
    decoys_desc: pd.DataFrame,
    descriptor_names: tuple[str, ...] = DESCRIPTOR_NAMES,
) -> PropertyBiasResult:
    """Compare per-descriptor distributions between actives and decoys.

    Parameters
    ----------
    actives_desc, decoys_desc:
        DataFrames as returned by ``MoleculeSet.descriptor_frame()``.
    """
    rows = {}
    for name in descriptor_names:
        a = actives_desc[name].to_numpy(dtype=float)
        d = decoys_desc[name].to_numpy(dtype=float)
        a = a[np.isfinite(a)]
        d = d[np.isfinite(d)]
        if len(a) < 2 or len(d) < 2:
            rows[name] = {"ks_stat": np.nan, "ks_pvalue": np.nan, "standardized_wasserstein": np.nan}
            continue

        ks = ks_2samp(a, d)
        pooled_std = np.sqrt((a.var(ddof=1) + d.var(ddof=1)) / 2.0)
        w = wasserstein_distance(a, d)
        std_w = w / pooled_std if pooled_std > 0 else np.nan

        rows[name] = {
            "ks_stat": ks.statistic,
            "ks_pvalue": ks.pvalue,
            "standardized_wasserstein": std_w,
        }

    per_descriptor = pd.DataFrame(rows).T
    score = float(per_descriptor["ks_stat"].dropna().mean()) if per_descriptor["ks_stat"].notna().any() else float("nan")
    return PropertyBiasResult(per_descriptor=per_descriptor, score=score)
