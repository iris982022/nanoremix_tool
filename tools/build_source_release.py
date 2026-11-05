"""Build a minimal source-only upload package, excluding author working documents."""

import argparse
import hashlib
import json
from pathlib import Path
import shutil
import zipfile


ROOT_FILES = (
    "README.md", "pyproject.toml", "constraints-models.txt", ".gitignore", "LICENSE", "NOTICE.md", "SECURITY.md",
    ".github/workflows/tests.yml",
    "docs/README.zh-CN.md", "docs/INPUTS.md", "docs/MODELS.md", "docs/OPTIONS.md",
    "examples/make_demo.py", "examples/make_feature_demo.py", "examples/minimal_binary_probabilities.csv",
    "tools/build_source_release.py", "tools/smoke_installed.py", "tools/smoke_models.py",
)


def source_files(root):
    files = [root / name for name in ROOT_FILES]
    files.extend(sorted((root / "src/nanoremix").rglob("*.py")))
    files.extend(sorted((root / "tests").glob("*.py")))
    if not (root / "src/nanoremix/cli.py").is_file() or not (root / "tests/test_pipeline.py").is_file():
        raise ValueError("Source tree lacks the inference package or regression tests")
    for path in files:
        if (not path.is_file() or path.is_symlink()
                or not path.resolve().is_relative_to(root.resolve())):
            raise ValueError(f"Missing or nonlocal release file: {path}")
    return sorted(files)


def build(root, output):
    root, output = Path(root).resolve(), Path(output).resolve()
    files = source_files(root)
    if output.exists():
        raise ValueError("Release output already exists; choose a new directory")
    output.mkdir(parents=True)
    project = output / "NanoRemix"
    archive = output / "NanoRemix_GitHub_source.zip"
    manifest = []
    with zipfile.ZipFile(archive, "x", compression=zipfile.ZIP_DEFLATED) as bundle:
        for source in files:
            relative = source.relative_to(root)
            target = project / relative
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source, target)
            data = target.read_bytes()
            info = zipfile.ZipInfo((Path("NanoRemix") / relative).as_posix())
            info.compress_type = zipfile.ZIP_DEFLATED
            info.external_attr = 0o100644 << 16
            bundle.writestr(info, data)
            manifest.append({"path": relative.as_posix(), "bytes": len(data),
                             "sha256": hashlib.sha256(data).hexdigest()})
    report = {
        "package_kind": "source_only", "pretrained_weights_included": False,
        "real_model_inference_verified": False, "files": manifest,
        "file_count": len(manifest), "archive": archive.name,
        "archive_sha256": hashlib.sha256(archive.read_bytes()).hexdigest(),
        "excluded": ["author manuscript audits", "historical validation reports",
                     "release drafts", "paper/figures", "raw data", "model weights",
                     "Git metadata", "build products and caches"],
    }
    (output / "source_package_manifest.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8")
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", required=True, help="New directory for upload folder, ZIP and manifest")
    args = parser.parse_args()
    try:
        report = build(Path(__file__).resolve().parents[1], args.output_dir)
    except (ValueError, OSError) as error:
        parser.exit(2, f"{error}\n")
    print(json.dumps({k: v for k, v in report.items() if k not in ("files", "excluded")}, indent=2))


if __name__ == "__main__":
    main()
