#!/usr/bin/env python3
"""Verify the distribution's file hashes using only the Python standard library."""

import argparse
import hashlib
import json
from pathlib import Path, PurePosixPath
import sys

sys.dont_write_bytecode = True


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def verify(root):
    root = Path(root).resolve()
    checksum_file = root / "SHA256SUMS.txt"
    expected = {}
    errors = []
    for line in checksum_file.read_text(encoding="utf-8").splitlines():
        digest, separator, relative = line.partition("  ")
        path = PurePosixPath(relative)
        if (not separator or len(digest) != 64
                or any(c not in "0123456789abcdef" for c in digest)
                or not relative or path.is_absolute() or ".." in path.parts
                or relative == "SHA256SUMS.txt" or relative in expected):
            raise ValueError(f"Invalid checksum entry: {line!r}")
        expected[relative] = digest
        target = root / relative
        if not target.resolve().is_relative_to(root) or target.is_symlink():
            errors.append({"path": relative, "error": "nonlocal_or_symbolic_link"})
        elif not target.is_file():
            errors.append({"path": relative, "error": "missing_file"})
        elif sha256(target) != digest:
            errors.append({"path": relative, "error": "sha256_mismatch"})
    if not expected:
        raise ValueError("Empty checksum list")
    manifest = json.loads((root / "MANIFEST.json").read_text(encoding="utf-8"))
    manifest_files = {row["path"]: row for row in manifest["files"]}
    if len(manifest_files) != len(manifest["files"]):
        errors.append({"error": "duplicate_manifest_entries"})
    if set(expected) != set(manifest_files) | {"MANIFEST.json"}:
        errors.append({"error": "checksum_manifest_file_sets_differ"})
    for relative, row in manifest_files.items():
        target = root / relative
        if expected.get(relative) != row["sha256"]:
            errors.append({"path": relative, "error": "manifest_sha256_mismatch"})
        if target.is_file() and target.stat().st_size != row["bytes"]:
            errors.append({"path": relative, "error": "manifest_size_mismatch"})
    distribution = json.loads((root / "DISTRIBUTION.json").read_text(encoding="utf-8"))
    weights = [name for name in expected if name.startswith("models/") and Path(name).suffix in (".pt", ".pkl")]
    if len(weights) != (89 if distribution["profile"] == "full" else 0):
        errors.append({"error": "unexpected_weight_count", "count": len(weights)})
    actual = {p.relative_to(root).as_posix() for p in root.rglob("*") if p.is_file()}
    extras = sorted(actual - set(expected) - {"SHA256SUMS.txt"})
    return {"status": "FAIL" if errors else "PASS", "profile": distribution["profile"],
            "nanoremix_version": distribution["nanoremix_version"],
            "checksummed_files": len(expected), "model_weights": len(weights),
            "extra_files_count": len(extras), "extra_files": extras[:20],
            "errors": errors, "models_deserialized": False,
            "scope": "File integrity only; not preprocessing or biological validation."}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    args = parser.parse_args()
    try:
        report = verify(args.root)
    except (OSError, ValueError, KeyError) as error:
        parser.exit(2, f"Verification failed: {error}\n")
    print(json.dumps(report, ensure_ascii=False, indent=2))
    return 0 if report["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
