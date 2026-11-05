import importlib.util
import json
from pathlib import Path
import re
import zipfile

import pytest


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("source_release", ROOT / "tools/build_source_release.py")
RELEASE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(RELEASE)


def test_release_excludes_author_documents_and_weights(tmp_path):
    report = RELEASE.build(ROOT, tmp_path / "release")
    package = tmp_path / "release/NanoRemix"
    assert report["package_kind"] == "source_only"
    assert not report["pretrained_weights_included"]
    assert sorted(path.name for path in (package / "docs").iterdir()) == [
        "INPUTS.md", "MODELS.md", "OPTIONS.md", "README.zh-CN.md"]
    assert not (package / "models").exists()
    assert not (package / ".git").exists()
    assert (package / ".github/workflows/tests.yml").is_file()
    assert (package / "src/nanoremix/_vendor/multi9/train_9mer_limit_s3.py").is_file()
    assert not any(path.suffix in (".pdf", ".pkl", ".pt", ".pyc") for path in package.rglob("*"))
    with zipfile.ZipFile(tmp_path / "release" / report["archive"]) as archive:
        assert len(archive.namelist()) == report["file_count"]
        assert all(name.startswith("NanoRemix/") for name in archive.namelist())
    assert json.loads((tmp_path / "release/source_package_manifest.json").read_text()) == report
    with pytest.raises(ValueError, match="already exists"):
        RELEASE.build(ROOT, tmp_path / "release")


def test_public_relative_document_links_resolve(tmp_path):
    RELEASE.build(ROOT, tmp_path / "release")
    root = tmp_path / "release/NanoRemix"
    for path in root.rglob("*.md"):
        for destination in re.findall(r"\[[^\]]+\]\(([^)]+)\)", path.read_text()):
            if destination.startswith(("#", "https://", "http://")):
                continue
            target = (path.parent / destination.split("#", 1)[0]).resolve()
            assert target.is_relative_to(root)
            assert target.is_file(), (path, destination)


def test_source_archive_is_reproducible(tmp_path):
    a = RELEASE.build(ROOT, tmp_path / "a")
    b = RELEASE.build(ROOT, tmp_path / "b")
    assert a["archive_sha256"] == b["archive_sha256"]
