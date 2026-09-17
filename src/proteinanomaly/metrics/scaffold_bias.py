"""Metric 2: scaffold-concentration / scaffold-leakage bias.

Two distinct failure modes are conflated under "scaffold bias" in the
literature, so we report them separately:

1. Scaffold leakage: actives and decoys occupy near-disjoint scaffold sets, so
   a trivial "which Bemis-Murcko scaffold is this" lookup would separate the
   classes almost perfectly. This is the classic DUD-E criticism: decoys are
   chosen to be topologically dissimilar to actives.
2. Scaffold concentration: within one class, how many distinct scaffolds
   actually carry the samples. Reported as ``1 - normalized_entropy`` per
   class (0 = every scaffold equally represented across as many buckets as
   there are molecules, 1 = a single scaffold dominates).

Naive symmetric Jaccard overlap on the raw scaffold sets is confounded by
class size: DUD-E/LIT-PCBA decoy sets outnumber actives 50-60:1, so the
decoy side alone accumulates far more *distinct* scaffolds simply by having
more molecules (the same species-accumulation effect rarefaction curves
correct for in ecology — more individuals sampled reveals more species,
independent of the true underlying diversity/overlap rate). That inflates
the Jaccard union and pushes the "leakage" score toward 1 even when the
smaller class's scaffolds are, proportionally, reasonably well represented
in the larger class. This was caught empirically: TP53 (79 actives / 4168
decoys) scored 0.985 on naive Jaccard, but 18 of the 46 active scaffolds
(39%) do occur somewhere in the 4168 decoys — the naive score mostly
reflects the 53x size gap, not a 39%-vs-100% distinction.

We therefore report the leakage score as a **rarefied** Jaccard: repeatedly
subsample the larger class's scaffold multiset down to the smaller class's
molecule count (without replacement), recompute Jaccard on each subsample,
and average. This matches denominators across targets/datasets with very
different active:decoy ratios or absolute sizes, so leakage scores are
comparable across e.g. AMPC (48 actives) and ALDH1 (5363 actives) rather than
mostly encoding "how many molecules did this target have."  The naive
(un-rarefied) Jaccard and the asymmetric active-scaffold coverage are kept
as diagnostic detail since they're each informative in their own right.
"""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass

import numpy as np


def _normalized_entropy(scaffolds: list[str]) -> float:
    if len(scaffolds) == 0:
        return float("nan")
    counts = np.array(list(Counter(scaffolds).values()), dtype=float)
    n_unique = len(counts)
    if n_unique <= 1:
        return 0.0
    probs = counts / counts.sum()
    entropy = -np.sum(probs * np.log(probs))
    return float(entropy / np.log(n_unique))


def _jaccard(a: set, b: set) -> float:
    if not a and not b:
        return float("nan")
    union = a | b
    if not union:
        return float("nan")
    return len(a & b) / len(union)


def _rarefied_jaccard(
    actives_scaffolds: list[str],
    decoys_scaffolds: list[str],
    n_repeats: int,
    random_state: int,
) -> tuple[float, float]:
    """Mean/std Jaccard overlap over repeated subsamples matching the
    smaller class's molecule count, removing the class-size confound."""
    n_a, n_d = len(actives_scaffolds), len(decoys_scaffolds)
    if n_a == 0 or n_d == 0:
        return float("nan"), float("nan")

    n_min = min(n_a, n_d)
    rng = np.random.default_rng(random_state)
    a_arr = np.array(actives_scaffolds)
    d_arr = np.array(decoys_scaffolds)

    scores = []
    for _ in range(n_repeats):
        a_sample = rng.choice(a_arr, size=n_min, replace=False) if n_a > n_min else a_arr
        d_sample = rng.choice(d_arr, size=n_min, replace=False) if n_d > n_min else d_arr
        scores.append(_jaccard(set(a_sample), set(d_sample)))
    scores = np.asarray(scores, dtype=float)
    return float(np.nanmean(scores)), float(np.nanstd(scores))


@dataclass
class ScaffoldBiasResult:
    score: float  # 1 - rarefied jaccard overlap, in [0, 1] — size-imbalance-controlled leakage score
    jaccard_overlap: float  # rarefied jaccard overlap mean (the value `score` is derived from)
    jaccard_overlap_std: float
    raw_jaccard_overlap: float  # naive jaccard on full (un-rarefied) scaffold sets — confounded by class size, kept for reference
    active_scaffold_coverage: float  # |A ∩ D| / |A|, asymmetric: fraction of active scaffolds seen anywhere among decoys
    decoy_scaffold_coverage: float  # |A ∩ D| / |D|, asymmetric: fraction of decoy scaffolds seen anywhere among actives
    n_unique_actives: int
    n_unique_decoys: int
    actives_concentration: float  # 1 - normalized entropy, in [0, 1]
    decoys_concentration: float

    def to_dict(self) -> dict:
        return {
            "score": self.score,
            "jaccard_overlap": self.jaccard_overlap,
            "jaccard_overlap_std": self.jaccard_overlap_std,
            "raw_jaccard_overlap": self.raw_jaccard_overlap,
            "active_scaffold_coverage": self.active_scaffold_coverage,
            "decoy_scaffold_coverage": self.decoy_scaffold_coverage,
            "n_unique_actives": self.n_unique_actives,
            "n_unique_decoys": self.n_unique_decoys,
            "actives_concentration": self.actives_concentration,
            "decoys_concentration": self.decoys_concentration,
        }


def scaffold_bias(
    actives_scaffolds: list[str],
    decoys_scaffolds: list[str],
    n_repeats: int = 200,
    random_state: int = 0,
) -> ScaffoldBiasResult:
    a_set, d_set = set(actives_scaffolds), set(decoys_scaffolds)
    intersection = len(a_set & d_set)
    raw_overlap = _jaccard(a_set, d_set)

    rarefied_mean, rarefied_std = _rarefied_jaccard(
        actives_scaffolds, decoys_scaffolds, n_repeats=n_repeats, random_state=random_state
    )
    score = float("nan") if np.isnan(rarefied_mean) else 1.0 - rarefied_mean

    return ScaffoldBiasResult(
        score=score,
        jaccard_overlap=rarefied_mean,
        jaccard_overlap_std=rarefied_std,
        raw_jaccard_overlap=raw_overlap,
        active_scaffold_coverage=(intersection / len(a_set)) if a_set else float("nan"),
        decoy_scaffold_coverage=(intersection / len(d_set)) if d_set else float("nan"),
        n_unique_actives=len(a_set),
        n_unique_decoys=len(d_set),
        actives_concentration=1.0 - _normalized_entropy(actives_scaffolds),
        decoys_concentration=1.0 - _normalized_entropy(decoys_scaffolds),
    )
