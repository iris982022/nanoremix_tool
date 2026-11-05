"""Exercise an installed package in a separate directory with synthetic inputs."""

import argparse
import importlib.metadata
import json
from pathlib import Path
import platform
import subprocess
import sys
import tempfile

import numpy as np
import pandas as pd

import nanoremix


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report", required=True)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    scenarios = [
        (mode, ["--mode", mode], count, "unspecified")
        for mode, count in (("thirdgen", 2), ("fused", 2), ("both", 4), ("secondgen", 1))
    ] + [
        ("biological", ["--data-type", "biological"], 2, "biological"),
        ("biological-third", ["--data-type", "biological", "--also-thirdgen"], 4, "biological"),
        ("biological-second", ["--data-type", "biological", "--also-secondgen"], 3, "biological"),
        ("biological-all", ["--data-type", "biological", "--also-thirdgen", "--also-secondgen"], 5, "biological"),
        ("synthetic", ["--data-type", "synthetic"], 2, "synthetic"),
    ]
    results = []
    with tempfile.TemporaryDirectory(prefix="nanoremix-wheel-smoke-") as scratch:
        scratch = Path(scratch)
        subprocess.run([sys.executable, str(root / "examples/make_demo.py"),
                        "--output", str(scratch / "inputs")], cwd=scratch, check=True,
                       capture_output=True, text=True)
        for chemistry in ("RNA002", "RNA004"):
            for task in ("binary", "multi"):
                for scenario, routing, expected_count, data_type in scenarios:
                    output = scratch / f"{chemistry}-{task}-{scenario}"
                    command = [
                        sys.executable, "-m", "nanoremix", "predict",
                        "--input", str(scratch / "inputs" / f"{chemistry}_{task}.csv"),
                        "--output", str(output), "--task", task, "--chemistry", "auto",
                        "--input-kind", "probabilities", "--mods", "m6A,m5C,pseudoU",
                        *routing, "--strategy", "both", "--chunk-size", "3",
                        "--truth-column", "true_label",
                    ]
                    subprocess.run(command, cwd=scratch, check=True, capture_output=True, text=True)
                    files = list(output.rglob("predictions.tsv"))
                    assert len(files) == expected_count
                    for file in files:
                        data = pd.read_csv(file, sep="\t")
                        columns = [c for c in data if c.startswith("prob_")]
                        assert set(columns) == {"prob_" + c for c in ("A", "m6A", "C", "m5C", "U", "pseudoU")}
                        assert len(data) == 16
                        assert data.data_type.eq(data_type).all()
                        valid = data.status.eq("ok")
                        np.testing.assert_allclose(data.loc[valid, columns].sum(axis=1), 1)
                    results.append({"chemistry": chemistry, "task": task, "scenario": scenario,
                                    "prediction_tables": len(files), "passed": True})
    report = {
        "software_only": True, "synthetic_data": True, "original_weights_available": False,
        "python": platform.python_version(),
        "versions": {x: importlib.metadata.version(x) for x in ("nanoremix", "numpy", "pandas")},
        "installed_wheel_not_source_tree": "/src/nanoremix/" not in str(nanoremix.__file__),
        "commands_passed": len(results), "prediction_tables_checked": sum(x["prediction_tables"] for x in results),
        "results": results,
    }
    Path(args.report).write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({k: v for k, v in report.items() if k != "results"}, indent=2))


if __name__ == "__main__":
    main()
