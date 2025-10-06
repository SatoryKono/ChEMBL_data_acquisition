"""Tests for dictionary resource utilities."""

from __future__ import annotations

import logging
from pathlib import Path, PureWindowsPath

import pytest
import yaml

from library.resources import dictionaries


def _create_sample_dictionary(tmp_path: Path) -> Path:
    root = tmp_path / "dictionary"
    (root / "subdir").mkdir(parents=True)
    (root / "subdir" / "data.csv").write_text("value\n", encoding="utf-8")
    (root / "root.txt").write_text("root\n", encoding="utf-8")
    return root


@pytest.mark.unit
def test_compute_sha256__stable_across_path_separators(tmp_path, monkeypatch):
    base = _create_sample_dictionary(tmp_path)
    expected = dictionaries._compute_sha256(base)

    original_relative_to = Path.relative_to

    def _relative_to_windows(self: Path, other: Path):  # pragma: no cover - behaviour covered via assertion
        return PureWindowsPath(original_relative_to(self, other))

    monkeypatch.setattr(Path, "relative_to", _relative_to_windows)

    assert dictionaries._compute_sha256(base) == expected


@pytest.mark.unit
def test_compute_sha256__independent_of_rglob_order(tmp_path, monkeypatch):
    base = _create_sample_dictionary(tmp_path)
    expected = dictionaries._compute_sha256(base)

    original_rglob = Path.rglob

    def _reversed_rglob(self: Path, pattern: str):  # pragma: no cover - behaviour asserted below
        return iter(list(original_rglob(self, pattern))[::-1])

    monkeypatch.setattr(Path, "rglob", _reversed_rglob)

    assert dictionaries._compute_sha256(base) == expected


@pytest.fixture()
def _temporary_manifest(tmp_path: Path) -> Path:
    resource_path = tmp_path / "data.csv"
    resource_path.write_text("value\n", encoding="utf-8")
    manifest = {
        "version": 1,
        "resources": {
            "example": {
                "path": "data.csv",
                "version": "test",
                "sha256": "deadbeef",
                "generator": "tests",
            }
        },
    }
    (tmp_path / "manifest.yaml").write_text(yaml.safe_dump(manifest), encoding="utf-8")
    return tmp_path


def _clear_manifest_cache() -> None:
    dictionaries._load_manifest.cache_clear()


@pytest.mark.unit
def test_manifest_checksum_mismatch__raises_without_override(
    _temporary_manifest: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    _clear_manifest_cache()
    monkeypatch.delenv("CHEMBL_DICTIONARY_ALLOW_MISMATCH", raising=False)
    with pytest.raises(dictionaries.DictionaryManifestError):
        dictionaries.get_resource("example", base_dir=_temporary_manifest)


@pytest.mark.unit
def test_manifest_checksum_mismatch__logs_warning_with_override(
    _temporary_manifest: Path, monkeypatch: pytest.MonkeyPatch, caplog: pytest.LogCaptureFixture
) -> None:
    _clear_manifest_cache()
    monkeypatch.setenv("CHEMBL_DICTIONARY_ALLOW_MISMATCH", "1")
    with caplog.at_level(logging.WARNING):
        resource = dictionaries.get_resource("example", base_dir=_temporary_manifest)
    assert resource.name == "example"
    assert resource.path == (_temporary_manifest / "data.csv").resolve()
    assert any(
        record.message == "dictionary_checksum_mismatch" and record.resource == "example"
        for record in caplog.records
    )
    monkeypatch.delenv("CHEMBL_DICTIONARY_ALLOW_MISMATCH", raising=False)
    _clear_manifest_cache()
