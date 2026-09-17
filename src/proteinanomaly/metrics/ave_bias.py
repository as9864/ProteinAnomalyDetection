"""Metric 3: AVE (Asymmetric Validation Embedding) bias — the inter-class
"nearest-neighbor memorizability" distance.

Follows Wallach & Heifets, "Most Ligand-Based Classification Benchmarks Reward
Memorization Rather Than Generalization" (J. Chem. Inf. Model. 2018). Under a
random train/test split of actives (A) and decoys/inactives (I):

    AVE = (NN_AA - NN_AI) + (NN_II - NN_IA)

where NN_XY is the mean, over molecules in the *test* half of class X, of the
Tanimoto similarity to their nearest neighbor in the *train* half of class Y.

AVE > 0 means a pure nearest-neighbor memorizer (no actual model) would score
above chance: test actives sit closer to train actives than to train decoys,
and vice versa for decoys — i.e. the benchmark rewards memorizing training
chemistry rather than learning real structure-activity signal. AVE ~ 0 means
the split is not trivially memorizable from fingerprint similarity alone.

We report the mean (and std) over repeated random splits, since a single
split is a noisy estimate.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from ..features import tanimoto_matrix


def _subsample(idx: np.ndarray, max_n: int, rng: np.random.Generator) -> np.ndarray:
    if len(idx) <= max_n:
        return idx
    return rng.choice(idx, size=max_n, replace=False)


def _mean_nn_similarity(query_fps: list, ref_fps: list) -> float:
    if len(query_fps) == 0 or len(ref_fps) == 0:
        return float("nan")
    sim = tanimoto_matrix(query_fps, ref_fps)
    return float(sim.max(axis=1).mean())


@dataclass
class AveBiasResult:
    score: float  # mean AVE bias across repeats
    std: float
    per_repeat: pd.DataFrame  # columns: nn_aa, nn_ai, nn_ii, nn_ia, bias

    def to_dict(self) -> dict:
        return {
            "score": self.score,
            "std": self.std,
            "per_repeat": self.per_repeat.to_dict(orient="records"),
        }


def ave_bias(
    actives_fps: list,
    decoys_fps: list,
    test_fraction: float = 0.5,
    n_repeats: int = 10,
    max_reference: int = 2000,
    random_state: int = 0,
) -> AveBiasResult:
    """Estimate AVE bias via repeated random train/test splits.

    ``max_reference`` caps the train-side reference set size per repeat so
    the O(n_test * n_train) Tanimoto matrix stays tractable on large decoy
    sets (DUD-E decoy sets can run into the tens of thousands per target).
    """
    rng = np.random.default_rng(random_state)
    n_a, n_d = len(actives_fps), len(decoys_fps)
    if n_a < 2 or n_d < 2:
        empty = pd.DataFrame(columns=["nn_aa", "nn_ai", "nn_ii", "nn_ia", "bias"])
        return AveBiasResult(score=float("nan"), std=float("nan"), per_repeat=empty)

    rows = []
    for _ in range(n_repeats):
        a_idx = rng.permutation(n_a)
        d_idx = rng.permutation(n_d)
        n_a_test = max(1, int(round(n_a * test_fraction)))
        n_d_test = max(1, int(round(n_d * test_fraction)))
        test_a_idx, train_a_idx = a_idx[:n_a_test], a_idx[n_a_test:]
        test_d_idx, train_d_idx = d_idx[:n_d_test], d_idx[n_d_test:]

        train_a_idx = _subsample(train_a_idx, max_reference, rng)
        train_d_idx = _subsample(train_d_idx, max_reference, rng)

        test_a_fps = [actives_fps[i] for i in test_a_idx]
        test_d_fps = [decoys_fps[i] for i in test_d_idx]
        train_a_fps = [actives_fps[i] for i in train_a_idx]
        train_d_fps = [decoys_fps[i] for i in train_d_idx]

        nn_aa = _mean_nn_similarity(test_a_fps, train_a_fps)
        nn_ai = _mean_nn_similarity(test_a_fps, train_d_fps)
        nn_ii = _mean_nn_similarity(test_d_fps, train_d_fps)
        nn_ia = _mean_nn_similarity(test_d_fps, train_a_fps)
        bias = (nn_aa - nn_ai) + (nn_ii - nn_ia)

        rows.append({"nn_aa": nn_aa, "nn_ai": nn_ai, "nn_ii": nn_ii, "nn_ia": nn_ia, "bias": bias})

    per_repeat = pd.DataFrame(rows)
    return AveBiasResult(
        score=float(per_repeat["bias"].mean()),
        std=float(per_repeat["bias"].std(ddof=1)) if len(per_repeat) > 1 else 0.0,
        per_repeat=per_repeat,
    )
