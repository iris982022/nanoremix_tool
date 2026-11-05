"""Build the documented source distribution with refreshed integrity metadata."""

import argparse
import hashlib
import json
from pathlib import Path, PurePosixPath
import shutil
import zipfile


def source_files(root):
    # The published manifest defines the source selection, including every manual
    # linked by README. Never glob a working tree containing data or credentials.
    manifest = json.loads((root / "MANIFEST.json").read_text(encoding="utf-8"))
    names = [row["path"] for row in manifest["files"]]
    if len(names) != len(set(names)):
        raise ValueError("Duplicate source manifest paths")
    for name in names:
        relative = PurePosixPath(name)
        if (not name or not relative.parts or relative.is_absolute() or ".." in relative.parts
                or str(relative) != name or "\\" in name
                or name in {"MANIFEST.json", "SHA256SUMS.txt"}
                or relative.parts[0] in {".git", "models", "data", "raw_data"}
                or relative.suffix in {".pkl", ".pickle", ".pt", ".pth", ".pyc"}):
            raise ValueError(f"Invalid source manifest path: {name}")
    files = [root / name for name in names]
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
    for source in files:
        target = project / source.relative_to(root)
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, target)
    distribution_file = project / "DISTRIBUTION.json"
    distribution = json.loads(distribution_file.read_text(encoding="utf-8"))
    distribution.update(project_file_count=len(files) + 2,
                        source_zip_contains_project_files=len(files) + 2,
                        source_archive=archive.name)
    distribution_file.write_text(json.dumps(distribution, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    selected = sorted(project / path.relative_to(root) for path in files)
    manifest = [{"path": path.relative_to(project).as_posix(), "bytes": path.stat().st_size,
                 "sha256": hashlib.sha256(path.read_bytes()).hexdigest()} for path in selected]
    (project / "MANIFEST.json").write_text(json.dumps({
        "profile": "source", "excluded_self_referential_files": ["MANIFEST.json", "SHA256SUMS.txt"],
        "files": manifest,
    }, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    selected = sorted([*selected, project / "MANIFEST.json"])
    (project / "SHA256SUMS.txt").write_text("".join(
        f"{hashlib.sha256(path.read_bytes()).hexdigest()}  {path.relative_to(project).as_posix()}\n"
        for path in selected), encoding="utf-8")
    selected = sorted([*selected, project / "SHA256SUMS.txt"])
    archive_manifest = []
    with zipfile.ZipFile(archive, "x", compression=zipfile.ZIP_DEFLATED) as bundle:
        for source in selected:
            relative = source.relative_to(project)
            data = source.read_bytes()
            info = zipfile.ZipInfo((Path("NanoRemix") / relative).as_posix())
            info.compress_type = zipfile.ZIP_DEFLATED
            info.external_attr = 0o100644 << 16
            bundle.writestr(info, data)
            archive_manifest.append({"path": relative.as_posix(), "bytes": len(data),
                             "sha256": hashlib.sha256(data).hexdigest()})
    report = {
        "package_kind": "source_only", "pretrained_weights_included": False,
        "real_model_inference_verified": False, "files": archive_manifest,
        "file_count": len(archive_manifest), "archive": archive.name,
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
