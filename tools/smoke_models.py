"""Check real pretrained assets with artificial inputs; not a biological benchmark."""

import argparse
import json
from pathlib import Path
import subprocess
import sys

import numpy as np
import pandas as pd

from nanoremix.assets import inspect_models
from nanoremix.pipeline import software_versions


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model-root", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--trust-model-files", required=True, action="store_true")
    args = parser.parse_args()
    root = Path(args.model_root).resolve()
    output = Path(args.output).resolve()
    if output.exists():
        raise ValueError("Smoke-test output exists; choose a new directory")
    inventory = inspect_models(root, model_checksums=root / "SHA256.json")
    if not inventory["preflight_passed"]:
        raise ValueError("Model preflight failed; run nanoremix doctor --strict first")
    output.mkdir(parents=True)
    repo = Path(__file__).resolve().parents[1]
    subprocess.run([sys.executable, str(repo / "examples/make_feature_demo.py"),
                    "--output", str(output / "inputs")], check=True)
    checks = []
    for chemistry in ("RNA002", "RNA004"):
        for task in ("binary", "multi"):
            for data_type in ("synthetic", "biological"):
                name = f"{chemistry}-{task}-{data_type}"
                target = output / name
                command = [
                    sys.executable, "-m", "nanoremix", "predict",
                    "--input", str(output / f"inputs/{chemistry}_features.csv"),
                    "--output", str(target), "--task", task, "--chemistry", chemistry,
                    "--mods", "all", "--data-type", data_type, "--strategy", "both",
                    "--model-root", str(root), "--model-checksums", str(root / "SHA256.json"),
                    "--trust-model-files",
                ]
                if data_type == "biological":
                    command.extend(["--sequence-column", "sequence81", "--also-thirdgen", "--also-secondgen"])
                run = subprocess.run(command, capture_output=True, text=True)
                (output / f"{name}.log").write_text(run.stdout + run.stderr)
                if run.returncode:
                    raise RuntimeError(f"Model smoke failed: {name}; inspect {output / (name + '.log')}")
                predictions = sorted(target.rglob("predictions.tsv"))
                assert len(predictions) == (5 if data_type == "biological" else 2)
                for path in predictions:
                    frame = pd.read_csv(path, sep="\t")
                    assert len(frame) == 8
                    expected = np.where((frame.center_ref_base == "G") & (task == "binary"),
                                        "outside_selected_family", "ok")
                    np.testing.assert_array_equal(frame.status, expected)
                    valid = frame.status.eq("ok")
                    cols = [c for c in frame if c.startswith("prob_")]
                    probs = frame.loc[valid, cols].to_numpy()
                    assert np.isfinite(probs).all() and (probs >= 0).all() and (probs <= 1).all()
                    np.testing.assert_allclose(probs.sum(1), 1, atol=1e-6)
                if data_type == "biological":
                    for strategy in ("s1", "s2"):
                        third = pd.read_csv(target / f"thirdgen/{strategy}/predictions.tsv", sep="\t")
                        independent = pd.read_csv(output / f"{chemistry}-{task}-synthetic/thirdgen/{strategy}/predictions.tsv", sep="\t")
                        cols = [c for c in third if c.startswith("prob_")]
                        np.testing.assert_allclose(third[cols], independent[cols], equal_nan=True, atol=1e-12)
                checks.append({"route": name, "passed": True, "prediction_tables": len(predictions)})
                print(f"passed: {name}", flush=True)
    report = {
        "pretrained_model_loading_and_prediction_passed": True,
        "input_type": "artificial features and artificial 81-nt sequences",
        "biological_accuracy_evaluated": False, "paper_equivalence_verified": False,
        "environment": software_versions(), "asset_count": len(inventory["files"]),
        "commands_passed": len(checks), "prediction_tables": sum(x["prediction_tables"] for x in checks),
        "checks": checks,
    }
    (output / "model_smoke_report.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
