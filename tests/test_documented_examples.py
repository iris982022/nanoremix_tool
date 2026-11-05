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
