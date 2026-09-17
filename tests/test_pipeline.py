from proteinanomaly.pipeline import diagnose_target, summarize


def test_diagnose_target_end_to_end(actives, decoys_biased):
    report = diagnose_target("toy_target", actives, decoys_biased, ave_n_repeats=5)
    assert report.n_actives == 8
    assert report.n_decoys == 8
    assert report.property_bias_score > 0
    assert report.scaffold_bias_score == 1.0
    row = report.to_summary_row()
    assert row["target"] == "toy_target"


def test_summarize_multiple_targets(actives, decoys_matched, decoys_biased):
    r1 = diagnose_target("matched_target", actives, decoys_matched, ave_n_repeats=5)
    r2 = diagnose_target("biased_target", actives, decoys_biased, ave_n_repeats=5)
    summary = summarize([r1, r2])
    assert list(summary.index) == ["matched_target", "biased_target"]
    assert summary.loc["biased_target", "scaffold_bias"] > summary.loc["matched_target", "scaffold_bias"]
