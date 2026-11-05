import json
from pathlib import Path

import pytest

from nanoremix.assets import inspect_models
from nanoremix.assets import dependency_import_check
from nanoremix.cli import main
from nanoremix.labels import NATIVE, code_strategy
from nanoremix.tables import sha256


def scope(**overrides):
    return {"chemistry": "RNA004", "task": "binary", "mods": "m6A",
            "strategy": "s2", "data_type": "synthetic", **overrides}


def binary_file(root, data=b"test bytes, NOT pretrained weights"):
    path = root / "RNA004/binary/binary_classify/A_m6A/strategy3/m_all_calibrated.pkl"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(data)
    return path


def multi_files(root, chemistry="RNA004", strategy="s2"):
    k = 5 if chemistry == "RNA002" else 9
    count = k if strategy == "s1" else 1
    directory = root / chemistry / "multi" / strategy
    directory.mkdir(parents=True)
    components = []
    for i in range(count):
        component = directory / f"component{i}.pkl"
        component.write_bytes(b"not a model; static checks must not unpickle")
        components.append({"pickle": component.name, "sha256": sha256(component)})
    manifest = {
        "format": "nanoremix-thirdgen-standalone-v1",
        "module_alias": f"nr_six_v3_{'5mer' if k == 5 else '9mer_limit'}_{code_strategy(strategy)}",
        "k": k, "strategy": code_strategy(strategy), "classes": list(NATIVE[chemistry]),
        "schemas": [{} for _ in range(count)], "components": components,
        "fusion": {"gate": True, "families": {
            base: {"method": "arithmetic", "weights": [1 / count] * count, "temperature": 1}
            for base in "ACGU"}},
    }
    path = directory / "manifest.json"
    path.write_text(json.dumps(manifest))
    return path, manifest


def test_empty_root_is_not_ready(tmp_path):
    result = inspect_models(tmp_path)
    assert not result["files_complete"]
    assert not result["preflight_passed"]
    assert not result["real_model_inference_verified"]
    assert len(result["files"]) == 78  # Multi component paths are not known yet.
    assert all(item["status"] == "missing" for item in result["files"])


def test_directory_without_weight_is_missing(tmp_path):
    (tmp_path / "RNA004/binary/binary_classify").mkdir(parents=True)
    result = inspect_models(tmp_path, **scope())
    assert len(result["files"]) == 1
    assert result["files"][0]["status"] == "missing"


@pytest.mark.parametrize("data,status", [
    (b"", "empty"),
    (b"version https://git-lfs.github.com/spec/v1\noid sha256:abc\nsize 50\n", "git_lfs_pointer"),
])
def test_unusable_file_is_not_a_weight(tmp_path, data, status):
    binary_file(tmp_path, data)
    result = inspect_models(tmp_path, **scope())
    assert not result["files_complete"]
    assert result["files"][0]["status"] == status


def test_present_file_is_not_claimed_as_verified_inference(tmp_path):
    binary_file(tmp_path)
    result = inspect_models(tmp_path, **scope())
    assert result["files_complete"]
    assert not result["real_model_inference_verified"]
    assert all("secondgen" not in item["path"] for item in result["files"])


@pytest.mark.parametrize("chemistry,count", [("RNA002", 5), ("RNA004", 9)])
def test_binary_s1_checks_all_requested_positions(tmp_path, chemistry, count):
    result = inspect_models(tmp_path, **scope(chemistry=chemistry, strategy="s1"))
    assert len(result["files"]) == count
    assert all("A_m6A/strategy1/" in item["path"] for item in result["files"])


def test_binary_ablation_does_not_accept_full_weights(tmp_path):
    binary_file(tmp_path)
    result = inspect_models(tmp_path, **scope(condition="mean"))
    assert not result["files_complete"]
    assert "/only_mean/" in result["files"][0]["path"]


def test_biological_sequence_requires_all_checkpoints(tmp_path):
    binary_file(tmp_path)
    result = inspect_models(tmp_path, **scope(data_type="biological"))
    assert not result["files_complete"]
    assert sum(item["path"].endswith("best_model.pt") for item in result["files"]) == 25
    assert any(item["path"].endswith("scripts/train.py") for item in result["files"])


def test_existing_second_scores_skip_sequence_network(tmp_path):
    binary_file(tmp_path)
    result = inspect_models(tmp_path, **scope(data_type="biological", second_source="probabilities"))
    assert result["files_complete"]
    assert len(result["files"]) == 1


def test_trusted_checksum_index_is_enforced(tmp_path):
    path = binary_file(tmp_path)
    index = tmp_path / "SHA256.json"
    relative = path.relative_to(tmp_path).as_posix()
    index.write_text(json.dumps({relative: sha256(path)}))
    assert inspect_models(tmp_path, **scope(model_checksums=index))["files_complete"]
    index.write_text(json.dumps({relative: "0" * 64}))
    assert inspect_models(tmp_path, **scope(model_checksums=index))["files"][0]["status"] == "sha256_mismatch"
    index.write_text(json.dumps({"other.pkl": "0" * 64}))
    assert inspect_models(tmp_path, **scope(model_checksums=index))["files"][0]["status"] == "not_covered_by_checksums"


@pytest.mark.parametrize("chemistry", ["RNA002", "RNA004"])
@pytest.mark.parametrize("strategy", ["s1", "s2"])
def test_multi_component_inventory_without_unpickling(tmp_path, chemistry, strategy):
    _, manifest = multi_files(tmp_path, chemistry, strategy)
    result = inspect_models(tmp_path, **scope(task="multi", chemistry=chemistry, strategy=strategy))
    assert result["files_complete"]
    assert len(result["files"]) == 1 + len(manifest["components"])
    assert not result["real_model_inference_verified"]


@pytest.mark.parametrize("change", ["classes", "temperature", "weights", "components", "escape", "duplicate"])
def test_incompatible_multi_metadata_rejected(tmp_path, change):
    path, manifest = multi_files(tmp_path, strategy="s1")
    if change == "classes":
        manifest["classes"].reverse()
    elif change == "temperature":
        manifest["fusion"]["families"]["A"]["temperature"] = 0
    elif change == "weights":
        manifest["fusion"]["families"]["A"]["weights"] = [float("nan")] * 9
    elif change == "components":
        manifest["components"].pop()
    elif change == "escape":
        manifest["components"][0]["pickle"] = "../../../outside.pkl"
    else:
        manifest["components"][1] = manifest["components"][0]
    path.write_text(json.dumps(manifest))
    result = inspect_models(tmp_path, **scope(task="multi", strategy="s1"))
    assert not result["files_complete"]
    assert result["groups"][0]["issues"]


def test_multi_component_digest_must_match(tmp_path):
    path, manifest = multi_files(tmp_path)
    (path.parent / manifest["components"][0]["pickle"]).write_bytes(b"modified")
    result = inspect_models(tmp_path, **scope(task="multi"))
    assert not result["files_complete"]
    assert result["files"][1]["status"] == "component_sha256_mismatch"


def test_symbolic_link_outside_model_root_rejected(tmp_path):
    root = tmp_path / "models"
    path = binary_file(root)
    outside = tmp_path / "outside.pkl"
    outside.write_bytes(b"external")
    path.unlink()
    path.symlink_to(outside)
    result = inspect_models(root, **scope())
    assert result["files"][0]["status"] == "outside_model_root"


def test_extra_checkpoint_rejected(tmp_path):
    binary_file(tmp_path)
    directory = tmp_path / "secondgen/models/final_25fold_ensemble_81bp_cnn_gru_attention"
    for repeat in range(1, 6):
        for fold in range(1, 6):
            checkpoint = directory / f"repeat_{repeat:02d}/fold_{fold:02d}/best_model.pt"
            checkpoint.parent.mkdir(parents=True)
            checkpoint.write_bytes(b"static test")
    extra = directory / "unexpected/best_model.pt"
    extra.parent.mkdir()
    extra.write_bytes(b"extra")
    source = tmp_path / "secondgen/scripts/train.py"
    source.parent.mkdir(parents=True)
    source.write_text("raise RuntimeError('doctor must never execute this')")
    result = inspect_models(tmp_path, **scope(data_type="biological"))
    assert not result["files_complete"]
    assert result["groups"][-1]["issues"]


def test_strict_cli_reports_failure_as_json(tmp_path, capsys):
    assert main(["doctor", "--model-root", str(tmp_path), "--strict"]) == 2
    report = json.loads(capsys.readouterr().out)
    assert not report["preflight_passed"]
    assert main(["doctor", "--model-root", str(tmp_path)]) == 0


def test_invalid_checksum_file_gives_friendly_error(tmp_path, capsys):
    index = tmp_path / "SHA256.json"
    index.write_text('{"../outside.pkl": "not-a-digest"}')
    assert main(["doctor", "--model-root", str(tmp_path), "--model-checksums", str(index)]) == 2
    assert "Invalid relative model path or SHA256" in capsys.readouterr().err


def test_dependency_presence_is_not_importability():
    report = dependency_import_check(["json", "nanoremix_nonexistent_module"])
    assert report["performed"] and not report["passed"]
    assert report["modules"]["json"]["loaded"]
    assert not report["modules"]["nanoremix_nonexistent_module"]["loaded"]


def test_doctor_refuses_a_broken_native_dependency(tmp_path, monkeypatch, capsys):
    binary_file(tmp_path)
    monkeypatch.setattr("nanoremix.assets.importlib.util.find_spec", lambda name: object())
    monkeypatch.setattr("nanoremix.assets.dependency_import_check", lambda names: {
        "performed": True, "passed": False, "error": "OpenMP unavailable"})
    assert main(["doctor", "--model-root", str(tmp_path), "--task", "binary",
                 "--chemistry", "RNA004", "--mods", "m6A", "--strategy", "s2",
                 "--data-type", "synthetic", "--strict"]) == 2
    report = json.loads(capsys.readouterr().out)
    assert report["files_complete"]
    assert not report["preflight_passed"]
