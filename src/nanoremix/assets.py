"""Inspect required inference assets without deserializing models or importing them."""

import importlib.util
import json
import math
from pathlib import Path
import subprocess
import sys

from . import __version__
from .labels import BASE, BINARY, NATIVE, code_strategy, resolve_labels
from .tables import sha256


def read_checksums(path):
    values = json.loads(Path(path).read_text(encoding="utf-8"))
    if not isinstance(values, dict) or not values:
        raise ValueError("Model checksums must be a nonempty JSON object")
    for name, digest in values.items():
        relative = Path(name)
        if (not name or relative.is_absolute() or ".." in relative.parts
                or relative.as_posix() != name or "\\" in name
                or not isinstance(digest, str) or len(digest) != 64
                or any(c not in "0123456789abcdef" for c in digest)):
            raise ValueError(f"Invalid relative model path or SHA256: {name!r}")
    return values


def dependency_import_check(names):
    code = (
        "import importlib, json, sys\n"
        "result = {}\n"
        "for name in sys.argv[1:]:\n"
        "    try:\n"
        "        module = importlib.import_module(name)\n"
        "        result[name] = {'loaded': True, 'version': getattr(module, '__version__', None)}\n"
        "    except Exception as error:\n"
        "        result[name] = {'loaded': False, 'error': str(error)}\n"
        "print(json.dumps(result))\n"
    )
    try:
        process = subprocess.run([sys.executable, "-c", code, *names],
                                 capture_output=True, text=True, timeout=120)
        if process.returncode:
            return {"performed": True, "passed": False, "error": process.stderr[-4000:]}
        modules = json.loads(process.stdout)
        return {"performed": True, "passed": all(info["loaded"] for info in modules.values()),
                "modules": modules}
    except (subprocess.TimeoutExpired, ValueError, OSError) as error:
        return {"performed": True, "passed": False, "error": str(error)}


def _multi_metadata(manifest, chemistry, strategy):
    k = 5 if chemistry == "RNA002" else 9
    count = k if strategy == "s1" else 1
    original = code_strategy(strategy)
    expected = {
        "format": "nanoremix-thirdgen-standalone-v1",
        "module_alias": f"nr_six_v3_{'5mer' if k == 5 else '9mer_limit'}_{original}",
        "k": k, "strategy": original, "classes": list(NATIVE[chemistry]),
    }
    if not isinstance(manifest, dict) or any(manifest.get(key) != value for key, value in expected.items()):
        raise ValueError("Multi manifest format, chemistry, strategy or class order mismatch")
    components = manifest.get("components")
    schemas = manifest.get("schemas")
    if (not isinstance(components, list) or len(components) != count
            or not isinstance(schemas, list) or len(schemas) != count):
        raise ValueError(f"Multi manifest requires {count} components and feature schemas")
    names = []
    for component in components:
        if not isinstance(component, dict):
            raise ValueError("Invalid Multi component record")
        name, digest = component.get("pickle"), component.get("sha256")
        if (not isinstance(name, str) or not name or Path(name).is_absolute()
                or ".." in Path(name).parts or "\\" in name
                or not isinstance(digest, str) or len(digest) != 64
                or any(c not in "0123456789abcdef" for c in digest)):
            raise ValueError("Multi components require relative pickle paths and SHA256 values")
        names.append(name)
    if len(set(names)) != count:
        raise ValueError("Multi components must have distinct file paths")
    fusion = manifest.get("fusion", {})
    if (not isinstance(fusion, dict) or not isinstance(fusion.get("gate"), bool)
            or not isinstance(fusion.get("families"), dict)
            or set(fusion["families"]) != set("ACGU")):
        raise ValueError("Multi manifest requires frozen family-fusion parameters")
    for config in fusion["families"].values():
        if not isinstance(config, dict):
            raise ValueError("Invalid Multi family-fusion record")
        weights = config.get("weights", [])
        temperature = config.get("temperature")
        if (not isinstance(weights, list) or len(weights) != count
                or any(type(w) not in (int, float) or not math.isfinite(w) or w < 0 for w in weights)
                or not math.isclose(sum(weights), 1.0, rel_tol=1e-5, abs_tol=1e-8)
                or type(temperature) not in (int, float)
                or not math.isfinite(temperature) or temperature <= 0
                or config.get("method") not in ("arithmetic", "geometric")):
            raise ValueError("Invalid Multi frozen weights, method or temperature")
    return components


def inspect_models(model_root="models", *, chemistry=None, task=None, mods="all",
                   strategy="both", data_type="biological", second_source="sequence",
                   condition="full", model_checksums=None, check_imports=False):
    root = Path(model_root or "models").expanduser().resolve()
    checksums = read_checksums(model_checksums) if model_checksums else None
    chemistries = (chemistry,) if chemistry else tuple(NATIVE)
    tasks = (task,) if task else ("binary", "multi")
    strategies = ("s1", "s2") if strategy == "both" else (strategy,)
    if condition != "full" and tasks != ("binary",):
        raise ValueError("Non-full conditions require --task binary")
    files, groups = {}, []

    def inspect(relative, expected=None):
        key = Path(relative).as_posix()
        if key in files:
            return files[key]
        path = root / relative
        entry = {"path": key, "status": "missing", "bytes": None, "sha256": None}
        files[key] = entry
        if not path.resolve().is_relative_to(root):
            entry["status"] = "outside_model_root"
        elif path.is_file():
            try:
                entry["bytes"] = path.stat().st_size
                with path.open("rb") as handle:
                    header = handle.read(128)
                if not header:
                    entry["status"] = "empty"
                elif header.startswith(b"version https://git-lfs.github.com/spec/v1"):
                    entry["status"] = "git_lfs_pointer"
                else:
                    entry["sha256"] = sha256(path)
                    if expected and entry["sha256"] != expected:
                        entry["status"] = "component_sha256_mismatch"
                    elif checksums is not None and key not in checksums:
                        entry["status"] = "not_covered_by_checksums"
                    elif checksums is not None and entry["sha256"] != checksums[key]:
                        entry["status"] = "sha256_mismatch"
                    else:
                        entry["status"] = "ok"
            except OSError as error:
                entry.update(status="unreadable", error=str(error))
        return entry

    for chem in chemistries:
        for current_task in tasks:
            selected = resolve_labels(chem, current_task, ("thirdgen",), "extended", mods=mods)
            for current_strategy in strategies:
                group = {"chemistry": chem, "task": current_task, "strategy": current_strategy,
                         "required_files": [], "issues": []}
                groups.append(group)
                if current_task == "binary":
                    k = 5 if chem == "RNA002" else 9
                    names = [f"m{i}" for i in range(1, k + 1)] if current_strategy == "s1" else ["m_all"]
                    for mod in BINARY:
                        if mod not in selected:
                            continue
                        directory = Path(chem) / "binary/binary_classify" / f"{BASE[mod]}_{mod}"
                        directory /= f"strategy{code_strategy(current_strategy)[-1]}"
                        if condition != "full":
                            directory /= f"only_{condition}"
                        for name in names:
                            entry = inspect(directory / f"{name}_calibrated.pkl")
                            group["required_files"].append(entry["path"])
                else:
                    directory = Path(chem) / "multi" / current_strategy
                    entry = inspect(directory / "manifest.json")
                    group["required_files"].append(entry["path"])
                    if entry["status"] != "ok":
                        group["issues"].append("Component paths cannot be checked until manifest.json is valid")
                        continue
                    try:
                        manifest = json.loads((root / entry["path"]).read_text(encoding="utf-8"))
                        components = _multi_metadata(manifest, chem, current_strategy)
                        for component in components:
                            if not (root / directory / component["pickle"]).resolve().is_relative_to(
                                    (root / directory).resolve()):
                                raise ValueError("Multi component escapes its model directory")
                            entry = inspect(directory / component["pickle"], component["sha256"])
                            group["required_files"].append(entry["path"])
                    except (ValueError, OSError) as error:
                        group["issues"].append(str(error))
    needs_second = data_type == "biological" and second_source == "sequence"
    if needs_second:
        paths = [Path("secondgen/scripts/train.py")]
        ensemble = Path("secondgen/models/final_25fold_ensemble_81bp_cnn_gru_attention")
        paths.extend(ensemble / f"repeat_{r:02d}/fold_{f:02d}/best_model.pt"
                     for r in range(1, 6) for f in range(1, 6))
        group = {"task": "secondgen", "strategy": "none", "required_files": [], "issues": []}
        groups.append(group)
        for path in paths:
            group["required_files"].append(inspect(path)["path"])
        if len(list((root / ensemble).rglob("best_model.pt"))) > 25:
            group["issues"].append("Second-generation loader requires exactly 25 checkpoints")
    for group in groups:
        group["complete"] = not group["issues"] and all(
            files[name]["status"] == "ok" for name in group["required_files"])
    required_dependencies = ["numpy", "pandas", "scipy", "sklearn", "xgboost", "joblib"]
    if needs_second:
        required_dependencies.append("torch")
    dependencies = {name: importlib.util.find_spec(name) is not None for name in required_dependencies}
    files_complete = all(group["complete"] for group in groups)
    preflight = files_complete and all(dependencies.values())
    imports = {"performed": False, "passed": None}
    if check_imports and preflight:
        imports = dependency_import_check(required_dependencies)
        preflight = imports["passed"]
    return {
        "version": __version__, "model_root": str(root),
        "scope": {"chemistry": list(chemistries), "task": list(tasks), "mods": mods,
                  "strategy": list(strategies), "data_type": data_type, "condition": condition,
                  "second_source": second_source if data_type == "biological" else "not_used"},
        "dependencies_present": dependencies,
        "dependency_import_check": imports,
        "files_complete": files_complete,
        "checksums_supplied": checksums is not None,
        "preflight_passed": preflight,
        "real_model_inference_verified": False,
        "note": "Static inventory only; no model was executed. Run validate with trusted assets and compare frozen expected predictions before release.",
        "groups": groups, "files": list(files.values()),
    }
