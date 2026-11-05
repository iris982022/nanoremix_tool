import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from nanoremix.cli import main, parser
from nanoremix.labels import BASE, NATIVE
from nanoremix.pipeline import run


def options(path, output, chemistry="RNA002", task="multi", mode="both", extra=()):
    return parser().parse_args([
        "predict", "--input", str(path), "--output", str(output), "--task", task,
        "--chemistry", chemistry, "--mods", "m6A,m5C,pseudoU", "--mode", mode,
        "--strategy", "both", "--input-kind", "probabilities", "--truth-column", "true_label",
        *extra,
    ])


@pytest.mark.parametrize("chemistry", ["RNA002", "RNA004"])
@pytest.mark.parametrize("task", ["binary", "multi"])
@pytest.mark.parametrize("mode", ["thirdgen", "fused", "both", "secondgen"])
def test_all_four_model_routes_and_modes(input_factory, tmp_path, chemistry, task, mode):
    path, _ = input_factory(chemistry, task)
    output = tmp_path / "out"
    report = run(options(path, output, chemistry, task, mode))
    expected = {"A", "m6A", "C", "m5C", "U", "pseudoU"}
    assert set(report["selected_labels"]) == expected
    for file in output.rglob("predictions.tsv"):
        predictions = pd.read_csv(file, sep="\t")
        cols = [c for c in predictions if c.startswith("prob_")]
        assert set(cols) == {"prob_" + c for c in expected}
        valid = predictions.status.eq("ok")
        np.testing.assert_allclose(predictions.loc[valid, cols].sum(1), 1, atol=1e-12)
        for c in expected:
            assert predictions.loc[valid & predictions.center_ref_base.ne(BASE[c]), "prob_" + c].eq(0).all()
        assert predictions.loc[predictions.center_ref_base.eq("G"), "status"].eq("outside_selected_family").all()
    assert report["model_retrained"] is False


@pytest.mark.parametrize("task", ["binary", "multi"])
def test_chunks_do_not_change_predictions(input_factory, tmp_path, task):
    path, _ = input_factory(task=task)
    run(options(path, tmp_path / "small", task=task, extra=["--chunk-size", "1"]))
    run(options(path, tmp_path / "large", task=task, extra=["--chunk-size", "100"]))
    for relative in ["thirdgen/s1", "thirdgen/s2", "fused/s1", "fused/s2"]:
        a = pd.read_csv(tmp_path / "small" / relative / "predictions.tsv", sep="\t")
        b = pd.read_csv(tmp_path / "large" / relative / "predictions.tsv", sep="\t")
        pd.testing.assert_frame_equal(a, b)


def test_third_only_does_not_read_missing_second_table(input_factory, tmp_path):
    path, _ = input_factory()
    args = options(path, tmp_path / "out", mode="thirdgen",
                   extra=["--second-probabilities", str(tmp_path / "not_here.csv")])
    report = run(args)
    assert report["second_source_used"] is False


def test_no_silent_fused_fallback(input_factory, tmp_path):
    path, frame = input_factory()
    frame.drop(columns=[c for c in frame if c.startswith("second_prob_")]).to_csv(path, index=False)
    with pytest.raises(ValueError, match="requires second"):
        run(options(path, tmp_path / "out", mode="fused"))
    assert not (tmp_path / "out").exists()
    assert not list(tmp_path.glob(".nanoremix-*"))


def test_rejects_mixed_chemistry(input_factory, tmp_path):
    path, frame = input_factory()
    frame.loc[8, "chemistry"] = "RNA004"
    frame.to_csv(path, index=False)
    with pytest.raises(ValueError, match="chemistry"):
        run(options(path, tmp_path / "out", extra=["--chunk-size", "3"]))


def test_duplicate_key_across_chunks_is_atomic_failure(input_factory, tmp_path):
    path, frame = input_factory()
    frame.loc[10, "record_key"] = frame.loc[0, "record_key"]
    frame.to_csv(path, index=False)
    with pytest.raises(ValueError, match="Duplicate"):
        run(options(path, tmp_path / "out", extra=["--chunk-size", "3"]))
    assert not (tmp_path / "out").exists()


def test_selection_does_not_leak_unselected_input_probabilities(input_factory, tmp_path):
    path, _ = input_factory()
    args = options(path, tmp_path / "out", mode="thirdgen", extra=["--keep-features"])
    args.mods = "m6A"
    run(args)
    output = pd.read_csv(tmp_path / "out/thirdgen/s1/predictions.tsv", sep="\t")
    assert {c for c in output if "prob_" in c} == {"prob_A", "prob_m6A"}


def test_truth_does_not_change_inference(input_factory, tmp_path):
    path, frame = input_factory()
    run(options(path, tmp_path / "a"))
    frame["true_label"] = "m1A"
    frame.to_csv(path, index=False)
    run(options(path, tmp_path / "b"))
    a = pd.read_csv(tmp_path / "a/fused/s2/predictions.tsv", sep="\t")
    b = pd.read_csv(tmp_path / "b/fused/s2/predictions.tsv", sep="\t")
    cols = ["predicted_label"] + [c for c in a if c.startswith("prob_")]
    pd.testing.assert_frame_equal(a[cols], b[cols])


def test_rejects_self_evaluating_truth_column(input_factory, tmp_path):
    path, _ = input_factory()
    args = options(path, tmp_path / "out")
    args.truth_column = "predicted_label"
    with pytest.raises(ValueError, match="Truth"):
        run(args)


def test_bad_probability_can_be_marked_and_other_mode_preserved(input_factory, tmp_path):
    path, frame = input_factory()
    frame.loc[0, "second_prob_A"] = -1
    frame.to_csv(path, index=False)
    run(options(path, tmp_path / "out", extra=["--on-invalid", "mark"]))
    pure = pd.read_csv(tmp_path / "out/thirdgen/s2/predictions.tsv", sep="\t")
    fused = pd.read_csv(tmp_path / "out/fused/s2/predictions.tsv", sep="\t")
    assert pure.loc[0, "status"] == "ok"
    assert fused.loc[0, "status"] == "invalid_second_probability"
    assert np.isnan(fused.loc[0, "prob_A"])


def test_low_confidence_retains_audit(input_factory, tmp_path):
    path, _ = input_factory()
    args = options(path, tmp_path / "out", extra=["--min-confidence", "1", "--output-filter", "modified"])
    run(args)
    audit = pd.read_csv(tmp_path / "out/thirdgen/s1/selection_audit.tsv", sep="\t")
    output = pd.read_csv(tmp_path / "out/thirdgen/s1/predictions.tsv", sep="\t")
    assert len(audit) == 12
    assert output.empty


def test_second_join_uses_keys_not_row_order(input_factory, tmp_path):
    path, frame = input_factory(task="binary")
    second_cols = ["site_id"] + [c for c in frame if c.startswith("second_prob_")]
    table = tmp_path / "second.tsv"
    frame[second_cols].iloc[::-1].to_csv(table, sep="\t", index=False)
    run(options(path, tmp_path / "inline", task="binary"))
    frame.drop(columns=second_cols[1:]).to_csv(path, index=False)
    run(options(path, tmp_path / "joined", task="binary",
                extra=["--second-probabilities", str(table), "--join-on", "site_id"]))
    a = pd.read_csv(tmp_path / "inline/fused/s2/predictions.tsv", sep="\t")
    b = pd.read_csv(tmp_path / "joined/fused/s2/predictions.tsv", sep="\t")
    pd.testing.assert_frame_equal(a, b)


def test_second_join_duplicate_fails(input_factory, tmp_path):
    path, frame = input_factory()
    second_cols = ["site_id"] + [c for c in frame if c.startswith("second_prob_")]
    table = tmp_path / "second.csv"
    pd.concat([frame[second_cols], frame[second_cols].iloc[[0]]]).to_csv(table, index=False)
    with pytest.raises(ValueError, match="unique"):
        run(options(path, tmp_path / "out", extra=["--second-probabilities", str(table), "--join-on", "site_id"]))


def test_existing_output_is_not_overwritten(input_factory, tmp_path):
    path, _ = input_factory()
    (tmp_path / "out").mkdir()
    with pytest.raises(ValueError, match="already exists"):
        run(options(path, tmp_path / "out"))


def test_validate_writes_no_output(input_factory, tmp_path):
    path, _ = input_factory()
    args = options(path, tmp_path / "out")
    args.command = "validate"
    run(args)
    assert not (tmp_path / "out").exists()


def test_s3_alias_and_auto_chemistry(input_factory, tmp_path):
    path, _ = input_factory("RNA004")
    args = options(path, tmp_path / "out", chemistry="auto", extra=["--strategy", "s3"])
    report = run(args)
    assert report["strategy"] == "s2" and report["chemistry"] == "RNA004"


def test_fused_exclusive_class_is_labelled_honestly(input_factory, tmp_path):
    path, _ = input_factory("RNA004")
    args = options(path, tmp_path / "out", chemistry="RNA004")
    args.mods = "hm5C"
    run(args)
    pure = pd.read_csv(tmp_path / "out/thirdgen/s1/predictions.tsv", sep="\t")
    fused = pd.read_csv(tmp_path / "out/fused/s1/predictions.tsv", sep="\t")
    rows = pure.status.eq("ok")
    assert fused.loc[rows, "evidence_source"].eq("thirdgen_only_in_selected_family").all()
    np.testing.assert_allclose(pure.loc[rows, ["prob_C", "prob_hm5C"]], fused.loc[rows, ["prob_C", "prob_hm5C"]])


def test_cli_unknown_mod_is_friendly_error(input_factory, tmp_path, capsys):
    path, _ = input_factory()
    code = main(["predict", "--input", str(path), "--output", str(tmp_path / "out"),
                 "--task", "multi", "--mods", "unknown", "--input-kind", "probabilities",
                 "--data-type", "synthetic"])
    assert code == 2
    assert "Unknown" in capsys.readouterr().err


def test_secondgen_mode_does_not_require_signal_features(input_factory, tmp_path):
    path, frame = input_factory("RNA004")
    frame.drop(columns=["9mer_final"]).to_csv(path, index=False)
    args = options(path, tmp_path / "out", chemistry="RNA004", mode="secondgen")
    args.input_kind = "features"
    report = run(args)
    assert report["second_source_used"]
    assert (tmp_path / "out/secondgen/none/predictions.tsv").is_file()


def test_changed_input_digest_aborts_atomically(input_factory, tmp_path, monkeypatch):
    import nanoremix.pipeline as pipeline
    from nanoremix.tables import sha256

    path, _ = input_factory()
    hashes = iter([sha256(path), "changed-during-inference"])
    monkeypatch.setattr(pipeline, "sha256", lambda _: next(hashes))
    with pytest.raises(ValueError, match="Input changed"):
        run(options(path, tmp_path / "out"))
    assert not (tmp_path / "out").exists()
