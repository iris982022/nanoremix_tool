import json
from pathlib import Path

import numpy as np
import pandas as pd

from nanoremix.cli import main


def quickstart(output, *extra):
    source = Path(__file__).resolve().parents[1] / "examples/minimal_binary_probabilities.csv"
    return main([
        "predict", "--input", str(source), "--task", "binary",
        "--chemistry", "RNA004", "--input-kind", "probabilities",
        "--data-type", "synthetic", "--mods", "m6A", "--strategy", "s2",
        "--output", str(output), *extra,
    ])


def test_readme_quickstart_probabilities_and_labels(tmp_path):
    output = tmp_path / "out"
    assert quickstart(output) == 0
    predictions = pd.read_csv(output / "thirdgen/s2/predictions.tsv", sep="\t")
    assert predictions.record_key.tolist() == ["demo-001", "demo-002"]
    assert predictions.predicted_label.tolist() == ["m6A", "A"]
    assert predictions.status.tolist() == ["ok", "ok"]
    np.testing.assert_allclose(
        predictions[["prob_A", "prob_m6A", "confidence"]],
        [[0.1, 0.9, 0.9], [0.8, 0.2, 0.8]],
    )
    report = json.loads((output / "run_manifest.json").read_text())
    assert report["input_rows"] == 2
    assert report["second_source_used"] is False
    assert report["model_assets_sha256"] == {}
    assert report["output_modes"] == ["thirdgen"]


def test_readme_example_rejections_remain_auditable(tmp_path):
    output = tmp_path / "out"
    assert quickstart(output, "--min-confidence", "0.85", "--output-filter", "modified") == 0
    predictions = pd.read_csv(output / "thirdgen/s2/predictions.tsv", sep="\t")
    audit = pd.read_csv(output / "thirdgen/s2/selection_audit.tsv", sep="\t")
    assert predictions.record_key.tolist() == ["demo-001"]
    assert len(audit) == 2
    assert audit.loc[1, "status"] == "below_confidence"
    assert audit.loc[1, "candidate_label"] == "A"
    assert pd.isna(audit.loc[1, "predicted_label"])
    assert not any(column.startswith("prob_") for column in audit.columns)


def test_step6_positive_retains_center_identity_and_signal(tmp_path):
    source = Path(__file__).resolve().parents[1] / "docs/guide_20261104/examples/RNA004_step6_positive.csv"
    raw = pd.read_csv(source)
    for i in range(1, 10):
        raw[f"mean_output{i}"] = np.arange(len(raw)) * 100 + i
    path, output = tmp_path / "step6.csv", tmp_path / "prepared.csv"
    raw.to_csv(path, index=False)
    assert main(["prepare", "--input", str(path), "--output", str(output),
                 "--chemistry", "RNA004", "--format", "step6"]) == 0
    prepared = pd.read_csv(output)
    assert len(prepared) == len(raw) - 8
    for _, candidate in prepared.iterrows():
        original = raw.loc[raw.ref_start.eq(candidate.ref_start)].iloc[0]
        assert candidate.read_id == original.read_id
        assert candidate["9mer_final"].replace("_", "")[4:13].replace("U", "T") == original.sequence_input
        for i in range(1, 10):
            assert candidate[f"mean_output{i}"] == original[f"mean_output{i}"]


def test_mixed_step6_strands_abort_without_partial_output(tmp_path, capsys):
    examples = Path(__file__).resolve().parents[1] / "docs/guide_20261104/examples"
    positive = pd.read_csv(examples / "RNA004_step6_positive.csv")
    negative = pd.read_csv(examples / "RNA004_step6_negative_limitation.csv")
    path, output = tmp_path / "mixed.csv", tmp_path / "prepared.csv"
    pd.concat([positive, negative], ignore_index=True).to_csv(path, index=False)
    assert main(["prepare", "--input", str(path), "--output", str(output),
                 "--chemistry", "RNA004", "--format", "step6"]) == 2
    assert "negative-strand rows" in capsys.readouterr().err
    assert not output.exists()


def test_step6_requires_explicit_valid_strand(tmp_path, capsys):
    source = Path(__file__).resolve().parents[1] / "docs/guide_20261104/examples/RNA004_step6_positive.csv"
    raw = pd.read_csv(source)
    path, output = tmp_path / "step6.csv", tmp_path / "prepared.csv"
    for frame in (raw.drop(columns="strand"), raw.assign(strand="unknown")):
        frame.to_csv(path, index=False)
        assert main(["prepare", "--input", str(path), "--output", str(output),
                     "--chemistry", "RNA004", "--format", "step6"]) == 2
        assert "strand" in capsys.readouterr().err
        assert not output.exists()


def test_step6_multiple_references_for_one_read_are_rejected(tmp_path, capsys):
    source = Path(__file__).resolve().parents[1] / "docs/guide_20261104/examples/RNA004_step6_positive.csv"
    raw = pd.read_csv(source)
    path, output = tmp_path / "multiple.csv", tmp_path / "prepared.csv"
    pd.concat([raw, raw.assign(ref_name="other-reference")], ignore_index=True).to_csv(path, index=False)
    assert main(["prepare", "--input", str(path), "--output", str(output),
                 "--chemistry", "RNA004", "--format", "step6"]) == 2
    assert "multiple references" in capsys.readouterr().err
    assert not output.exists()
