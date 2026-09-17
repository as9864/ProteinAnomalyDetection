from proteinanomaly.metrics import ave_bias, property_bias, scaffold_bias


def test_property_bias_low_for_matched_high_for_biased(actives, decoys_matched, decoys_biased):
    matched = property_bias(actives.descriptor_frame(), decoys_matched.descriptor_frame())
    biased = property_bias(actives.descriptor_frame(), decoys_biased.descriptor_frame())

    assert 0.0 <= matched.score <= 1.0
    assert 0.0 <= biased.score <= 1.0
    assert matched.score < biased.score
    assert biased.score > 0.5  # long alkanes vs. aromatic amines should be starkly separable


def test_scaffold_bias_zero_when_scaffolds_fully_shared(actives, decoys_matched):
    result = scaffold_bias(actives.scaffolds(), decoys_matched.scaffolds())
    assert result.jaccard_overlap == 1.0
    assert result.score == 0.0


def test_scaffold_bias_high_when_scaffolds_disjoint(actives, decoys_biased):
    result = scaffold_bias(actives.scaffolds(), decoys_biased.scaffolds())
    assert result.jaccard_overlap == 0.0
    assert result.score == 1.0


def test_scaffold_bias_rarefied_score_controls_for_class_size_imbalance():
    # 8 actives, all sharing one scaffold. 8 decoys share that same scaffold
    # exactly (full overlap), diluted by 200 decoys with 200 distinct never-
    # seen scaffolds -- simulating a large, diverse screening library where
    # the active chemotype is still fully represented, just outnumbered.
    active_scaffolds = ["shared"] * 8
    matching_decoys = ["shared"] * 8
    diverse_decoys = [f"unique_{i}" for i in range(200)]
    diluted_decoy_scaffolds = matching_decoys + diverse_decoys

    matched = scaffold_bias(active_scaffolds, matching_decoys, n_repeats=100, random_state=1)
    diluted = scaffold_bias(active_scaffolds, diluted_decoy_scaffolds, n_repeats=100, random_state=1)

    # every active scaffold is still present somewhere among the decoys
    assert diluted.active_scaffold_coverage == 1.0

    # the naive (un-rarefied) score collapses toward "fully separated" purely
    # because the decoy pool accumulated 200 extra distinct scaffolds
    naive_score = 1.0 - diluted.raw_jaccard_overlap
    assert naive_score > 0.9

    # the size-matched (rarefied) score is substantially less extreme than the
    # naive score, since it stops penalizing the target just for having a
    # larger, more diverse decoy pool
    assert diluted.score < naive_score

    # and the fully-overlapping, equal-size case is unaffected by rarefaction
    assert matched.score == 0.0


def test_ave_bias_higher_for_biased_than_matched(actives, decoys_matched, decoys_biased):
    matched = ave_bias(actives.fingerprints(), decoys_matched.fingerprints(), n_repeats=20, random_state=1)
    biased = ave_bias(actives.fingerprints(), decoys_biased.fingerprints(), n_repeats=20, random_state=1)

    assert matched.score < biased.score
    assert biased.score > 0  # nearest-neighbor memorization should trivially separate the biased case


def test_ave_bias_nan_when_class_too_small():
    result = ave_bias([], [], n_repeats=5)
    import math

    assert math.isnan(result.score)
