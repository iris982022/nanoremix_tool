"""Synthetic trained estimators test plumbing, not author-model equivalence."""

import importlib
import json
import pickle

import numpy as np
import pandas as pd
import pytest

pytest.importorskip("sklearn")
xgboost = pytest.importorskip("xgboost")

from nanoremix.cli import parser
from nanoremix.assets import inspect_models
from nanoremix.engines import Engines
from nanoremix.labels import NATIVE, code_strategy
from nanoremix.pipeline import run
from nanoremix.probability import binary_product, condition
from nanoremix.tables import sha256


def prepared(chemistry, count=72):
    rng = np.random.default_rng(981)
    k = 5 if chemistry == "RNA002" else 9
    sequence = rng.choice(list("ACGU"), (count, 2 * k - 1))
    sequence[:, k - 1] = "A"
    data = pd.DataFrame({
        "9mer_final": ["_".join(s) for s in sequence], "center_ref_base": "A",
        "record_key": [f"synthetic-{i}" for i in range(count)], "chemistry": chemistry,
    })
    for i in range(1, k + 1):
        for stat in ("mean", "median", "std", "length"):
            data[f"{stat}_output{i}"] = rng.uniform(.1, 2, count)
    return data


def options(tmp_path, chemistry, task, strategy):
    return parser().parse_args([
        "predict", "--input", str(tmp_path / "input.csv"), "--output", str(tmp_path / "output"),
        "--task", task, "--chemistry", chemistry, "--mods", "m6A", "--mode", "thirdgen",
        "--strategy", strategy, "--model-root", str(tmp_path / "models"), "--trust-model-files",
    ])


@pytest.mark.parametrize("chemistry", ["RNA002", "RNA004"])
@pytest.mark.parametrize("strategy", ["s1", "s2"])
def test_binary_trained_model_adapter(tmp_path, chemistry, strategy):
    from sklearn.linear_model import LogisticRegression

    k = 5 if chemistry == "RNA002" else 9
    features = importlib.import_module(f"nanoremix._vendor.binary{k}.features")
    data = prepared(chemistry)
    args = options(tmp_path, chemistry, "binary", strategy)
    data.to_csv(args.input, index=False)
    expanded = features.add_nine_mer_positions(data) if k == 5 else features.add_seventeen_mer_positions(data)
    directory = tmp_path / "models" / chemistry / "binary/binary_classify/A_m6A" / f"strategy{code_strategy(strategy)[-1]}"
    directory.mkdir(parents=True)
    expected = []
    for spec in features.binary_model_specs(f"strategy{code_strategy(strategy)[-1]}", "full"):
        matrix = features.build_feature_matrix(expanded, spec)
        y = np.arange(len(data)) % 2
        model = LogisticRegression(max_iter=1000).fit(matrix, y)
        with (directory / f"{spec.model_name}_calibrated.pkl").open("wb") as handle:
            pickle.dump(model, handle)
        expected.append(model.predict_proba(matrix)[:, 1])
    report = run(args)
    output = pd.read_csv(tmp_path / f"output/thirdgen/{strategy}/predictions.tsv", sep="\t")
    probability = binary_product(expected) if len(expected) > 1 else expected[0]
    np.testing.assert_allclose(output.prob_m6A, probability, atol=1e-12)
    assert len(report["model_assets_sha256"]) == (k if strategy == "s1" else 1)


@pytest.mark.parametrize("chemistry", ["RNA002", "RNA004"])
@pytest.mark.parametrize("strategy", ["s1", "s2"])
def test_multi_trained_xgboost_manifest_adapter(tmp_path, chemistry, strategy):
    k = 5 if chemistry == "RNA002" else 9
    module_name = f"train_5mer_{code_strategy(strategy)}" if k == 5 else f"train_9mer_limit_{code_strategy(strategy)}"
    module = importlib.import_module(f"nanoremix._vendor.multi{k}.{module_name}")
    data = prepared(chemistry, len(NATIVE[chemistry]) * 4)
    args = options(tmp_path, chemistry, "multi", strategy)
    data.to_csv(args.input, index=False)
    directory = tmp_path / "models" / chemistry / "multi" / strategy
    directory.mkdir(parents=True)
    schemas = module.schemas(k, code_strategy(strategy))
    expanded = module.prepare(data, k, list(NATIVE[chemistry]), training=False)
    components, predictions = [], []
    for schema in schemas:
        matrix = module.features(expanded, schema)
        y = np.arange(len(data)) % len(NATIVE[chemistry])
        model = xgboost.XGBClassifier(
            n_estimators=2, max_depth=2, n_jobs=1, tree_method="hist",
            objective="multi:softprob", num_class=len(NATIVE[chemistry]),
        ).fit(matrix, y)
        path = directory / f"{schema['name']}.pkl"
        with path.open("wb") as handle:
            pickle.dump(model, handle)
        components.append({"pickle": path.name, "sha256": sha256(path)})
        predictions.append(model.predict_proba(matrix))
    fusion = module.fixed_config(len(components), "arithmetic", temperature=.85)
    manifest = {
        "format": "nanoremix-thirdgen-standalone-v1", "module_alias": module.MODULE_ALIAS,
        "k": k, "strategy": code_strategy(strategy), "classes": list(NATIVE[chemistry]),
        "schemas": schemas, "components": components, "fusion": fusion,
    }
    (directory / "manifest.json").write_text(json.dumps(manifest))
    inventory = inspect_models(tmp_path / "models", chemistry=chemistry, task="multi",
                               strategy=strategy, mods="m6A", data_type="synthetic")
    assert inventory["files_complete"], inventory["groups"]
    run(args)
    raw = module.aggregate(np.asarray(predictions), np.array(["A"] * len(data)), list(NATIVE[chemistry]), fusion)
    expected, _, _ = condition(raw, NATIVE[chemistry], ("A", "m6A"), np.array(["A"] * len(data)))
    output = pd.read_csv(tmp_path / f"output/thirdgen/{strategy}/predictions.tsv", sep="\t")
    np.testing.assert_allclose(output[["prob_A", "prob_m6A"]], expected, atol=1e-7)


def test_model_loading_requires_explicit_trust(tmp_path):
    args = options(tmp_path, "RNA002", "binary", "s2")
    args.trust_model_files = False
    engine = Engines(args, "RNA002")
    with pytest.raises(ValueError, match="trust-model-files"):
        engine.trusted(tmp_path / "missing.pkl")


def test_checksum_mismatch_is_fatal(tmp_path):
    root = tmp_path / "models"
    root.mkdir()
    model = root / "asset.pkl"
    model.write_bytes(b"not a pickle; checksum verification must happen first")
    checksum = tmp_path / "sha.json"
    checksum.write_text(json.dumps({"asset.pkl": "0" * 64}))
    args = options(tmp_path, "RNA002", "binary", "s2")
    args.model_checksums = str(checksum)
    with pytest.raises(ValueError, match="SHA256"):
        Engines(args, "RNA002").trusted(model)
