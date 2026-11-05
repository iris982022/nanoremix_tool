import json

import numpy as np
import pandas as pd
import pytest

from nanoremix.cli import main, parser
from nanoremix.engines import Engines
from nanoremix.pipeline import run


def arguments(path, output, task="multi", routing=("--data-type", "biological"), extra=()):
    return [
        "predict", "--input", str(path), "--output", str(output),
        "--task", task, "--mods", "m6A,m5C,pseudoU",
        "--input-kind", "probabilities", "--strategy", "both",
        "--truth-column", "true_label", *routing, *extra,
    ]


POLICIES = [
    (["--data-type", "biological"], ["fused"]),
    (["--data-type", "biological", "--also-thirdgen"], ["fused", "thirdgen"]),
    (["--data-type", "biological", "--also-secondgen"], ["fused", "secondgen"]),
    (["--data-type", "biological", "--also-thirdgen", "--also-secondgen"],
     ["fused", "thirdgen", "secondgen"]),
    (["--data-type", "synthetic"], ["thirdgen"]),
]


@pytest.mark.parametrize("chemistry", ["RNA002", "RNA004"])
@pytest.mark.parametrize("task", ["binary", "multi"])
@pytest.mark.parametrize("routing,modes", POLICIES)
def test_data_type_output_matrix(input_factory, tmp_path, chemistry, task, routing, modes):
    path, _ = input_factory(chemistry, task)
    output = tmp_path / "out"
    report = run(parser().parse_args(arguments(path, output, task, routing)))
    expected = {
        f"{mode}/{strategy}" for mode in modes
        for strategy in (("none",) if mode == "secondgen" else ("s1", "s2"))
    }
    assert {p.parent.relative_to(output).as_posix() for p in output.rglob("predictions.tsv")} == expected
    assert report["data_type"] == routing[1]
    assert report["output_modes"] == tuple(modes)
    assert report["primary_mode"] == modes[0]
    assert report["mode"] is None
    assert report["second_source_used"] == ("fused" in modes)
    saved = json.loads((output / "run_manifest.json").read_text())
    assert saved["output_modes"] == modes
    assert saved["configuration"]["data_type"] == routing[1]
    for file in output.rglob("predictions.tsv"):
        frame = pd.read_csv(file, sep="\t")
        assert frame.data_type.eq(routing[1]).all()
        assert len(frame) == 12
        columns = [c for c in frame if c.startswith("prob_")]
        np.testing.assert_allclose(frame.loc[frame.status.eq("ok"), columns].sum(1), 1)


@pytest.mark.parametrize("routing", [
    [],
    ["--data-type", "biological", "--mode", "thirdgen"],
    ["--data-type", "synthetic", "--mode", "fused"],
])
def test_data_type_or_explicit_legacy_mode_is_required(tmp_path, routing):
    with pytest.raises(SystemExit):
        parser().parse_args(arguments(tmp_path / "input.csv", tmp_path / "out", routing=routing))


@pytest.mark.parametrize("routing", [
    ["--data-type", "synthetic", "--also-thirdgen"],
    ["--data-type", "synthetic", "--also-secondgen"],
    ["--mode", "fused", "--also-thirdgen"],
    ["--mode", "thirdgen", "--also-secondgen"],
])
def test_optional_generation_outputs_are_biological_only(input_factory, tmp_path, routing):
    path, _ = input_factory()
    with pytest.raises(ValueError, match="require --data-type biological"):
        run(parser().parse_args(arguments(path, tmp_path / "out", routing=routing)))
    assert not (tmp_path / "out").exists()


def test_synthetic_never_reads_second_inputs(input_factory, tmp_path, monkeypatch):
    path, frame = input_factory()
    frame.drop(columns=[c for c in frame if c.startswith("second_prob_")]).to_csv(path, index=False)

    def unexpected_second(*args):
        raise AssertionError("Synthetic workflow must not run the second generation")

    monkeypatch.setattr(Engines, "second", unexpected_second)
    options = arguments(path, tmp_path / "out", routing=["--data-type", "synthetic"],
                        extra=["--second-probabilities", str(tmp_path / "missing-second.csv")])
    report = run(parser().parse_args(options))
    assert report["output_modes"] == ("thirdgen",)
    assert report["second_source_used"] is False
    assert "second_probabilities_sha256" not in report


def test_biological_missing_second_input_does_not_fall_back(input_factory, tmp_path):
    path, frame = input_factory()
    frame.drop(columns=[c for c in frame if c.startswith("second_prob_")]).to_csv(path, index=False)
    options = arguments(path, tmp_path / "out", extra=["--also-thirdgen"])
    with pytest.raises(ValueError, match="requires second"):
        run(parser().parse_args(options))
    assert not (tmp_path / "out").exists()
    assert not list(tmp_path.glob(".nanoremix-*"))


@pytest.mark.parametrize("task", ["binary", "multi"])
def test_optional_outputs_preserve_each_generation_probability(input_factory, tmp_path, task):
    path, _ = input_factory(task=task)
    combined = tmp_path / "combined"
    run(parser().parse_args(arguments(path, combined, task, extra=["--also-thirdgen", "--also-secondgen"])))
    for mode in ("fused", "thirdgen", "secondgen"):
        standalone = tmp_path / mode
        run(parser().parse_args(arguments(path, standalone, task, routing=["--mode", mode])))
        for strategy in (("none",) if mode == "secondgen" else ("s1", "s2")):
            a = pd.read_csv(combined / mode / strategy / "predictions.tsv", sep="\t")
            b = pd.read_csv(standalone / mode / strategy / "predictions.tsv", sep="\t")
            pd.testing.assert_frame_equal(a.drop(columns="data_type"), b.drop(columns="data_type"))


def test_extra_second_output_does_not_repeat_sequence_inference(input_factory, tmp_path, monkeypatch):
    path, _ = input_factory()
    original = Engines.second
    calls = []

    def counted(engine, frame):
        calls.append(len(frame))
        return original(engine, frame)

    monkeypatch.setattr(Engines, "second", counted)
    run(parser().parse_args(arguments(path, tmp_path / "out", extra=["--also-thirdgen", "--also-secondgen"])))
    assert calls == [9]


@pytest.mark.parametrize("chemistry", ["RNA002", "RNA004"])
def test_all_modifications_keep_honest_second_output_scope(input_factory, tmp_path, chemistry):
    path, _ = input_factory(chemistry)
    args = parser().parse_args(arguments(path, tmp_path / "out", extra=["--also-secondgen", "--mods", "all"]))
    report = run(args)
    labels = report["selected_labels_by_mode"]["secondgen"]
    output = pd.read_csv(tmp_path / "out/secondgen/none/predictions.tsv", sep="\t")
    assert {c for c in output if c.startswith("prob_")} == {"prob_" + c for c in labels}
    assert "hm5C" in report["unsupported_secondgen_modifications"]
    assert "hm5C" in report["excluded_labels_by_mode"]["secondgen"]
    assert "prob_hm5C" not in output
    np.testing.assert_allclose(output.loc[output.status.eq("ok"), ["prob_" + c for c in labels]].sum(1), 1)


def test_unsupported_second_family_is_not_called_canonical(input_factory, tmp_path):
    path, _ = input_factory()
    args = parser().parse_args(arguments(path, tmp_path / "out", extra=["--also-secondgen", "--mods", "hm5C"]))
    report = run(args)
    output = pd.read_csv(tmp_path / "out/secondgen/none/predictions.tsv", sep="\t")
    rows = output.center_ref_base.eq("C")
    assert output.loc[rows, "status"].eq("unsupported_by_second_generation").all()
    assert output.predicted_label.isna().all()
    assert not any(c.startswith("prob_") for c in output)
    assert report["selected_labels_by_mode"]["secondgen"] == ()
    assert report["reports"]["secondgen/none"]["labels"] == ()
    assert report["reports"]["secondgen/none"]["conditional_accuracy"] is None
    fused = pd.read_csv(tmp_path / "out/fused/s2/predictions.tsv", sep="\t")
    assert fused.loc[rows, "status"].eq("ok").all()
    assert "prob_hm5C" in fused


def test_second_subset_excludes_unsupported_truth_from_evaluation(input_factory, tmp_path):
    path, frame = input_factory()
    frame.loc[frame.center_ref_base.eq("C"), "true_label"] = "hm5C"
    frame.to_csv(path, index=False)
    args = parser().parse_args(arguments(path, tmp_path / "out",
                                        extra=["--also-secondgen", "--mods", "m5C,hm5C"]))
    report = run(args)
    counts = report["reports"]["secondgen/none"]["status_counts_and_evaluation_counts"]
    assert counts["truth_outside_selected_scope"] == 12
    assert report["selected_labels_by_mode"]["secondgen"] == ("C", "m5C")
    output = pd.read_csv(tmp_path / "out/secondgen/none/predictions.tsv", sep="\t")
    assert "prob_hm5C" not in output
    assert output.loc[output.status.eq("ok"), "evidence_source"].eq("secondgen_supported_selected_subset").all()


def test_second_subset_cannot_leave_one_modified_competitor(input_factory, tmp_path):
    path, _ = input_factory()
    args = parser().parse_args(arguments(path, tmp_path / "out", extra=["--also-secondgen"]))
    args.mods, args.targets = None, "m5C,hm5C"
    report = run(args)
    assert report["selected_labels_by_mode"]["secondgen"] == ()
    output = pd.read_csv(tmp_path / "out/secondgen/none/predictions.tsv", sep="\t")
    assert output.predicted_label.isna().all()


def test_input_data_type_conflict_is_rejected_across_chunks(input_factory, tmp_path):
    path, frame = input_factory()
    frame["data_type"] = "biological"
    frame.loc[10, "data_type"] = "synthetic"
    frame.to_csv(path, index=False)
    options = arguments(path, tmp_path / "out", extra=["--chunk-size", "3"])
    with pytest.raises(ValueError, match="Input data_type conflicts"):
        run(parser().parse_args(options))
    assert not (tmp_path / "out").exists()


def test_validate_obeys_policy_without_saving_predictions(input_factory, tmp_path):
    path, _ = input_factory()
    args = parser().parse_args(arguments(path, tmp_path / "out", extra=["--also-secondgen"]))
    args.command = "validate"
    report = run(args)
    assert set(report["reports"]) == {"fused/s1", "fused/s2", "secondgen/none"}
    assert not (tmp_path / "out").exists()


def test_cli_reports_resolved_output_modes(input_factory, tmp_path, capsys):
    path, _ = input_factory()
    assert main(arguments(path, tmp_path / "out", extra=["--also-secondgen"])) == 0
    result = json.loads(capsys.readouterr().out)
    assert result["data_type"] == "biological"
    assert result["primary_mode"] == "fused"
    assert result["output_modes"] == ["fused", "secondgen"]
