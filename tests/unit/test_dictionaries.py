"""Tests for dictionary manifest validation utilities."""

from __future__ import annotations

import textwrap
from pathlib import Path, PureWindowsPath

import pytest

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


def _write_manifest(tmp_path, body: str) -> None:
    manifest_path = tmp_path / "manifest.yaml"
    manifest_path.write_text(textwrap.dedent(body), encoding="utf-8")


@pytest.mark.unit
@pytest.mark.parametrize(
    ("resource_name", "checksum"),
    [
        *[
            ("dictionary_root", checksum)
            for checksum in dictionaries.WINDOWS_DICTIONARY_ROOT_CHECKSUMS
        ],
        (
            "target_uniprot_cache",
            "c86b314b5d8a0906f1174c8e9f494cf9dde6841be2cb1e8b66c5772976afb5ca",
        ),
    ],
)
def test_parse_manifest__accepts_known_checksum_variants(
    tmp_path, monkeypatch, resource_name, checksum
):
    # Arrange
    _write_manifest(
        tmp_path,
        f"""
        version: 1
        resources:
          {resource_name}:
            path: .
            version: "test"
            sha256:
              - "deadbeef"
            generator: tools/build_dictionary_resources.py
        """,
    )
    (tmp_path / "placeholder.txt").write_text("data", encoding="utf-8")
    monkeypatch.setattr(dictionaries, "_compute_sha256", lambda path, value=checksum: value)

    # Act
    dictionaries._env_checksum_allowlist.cache_clear()
    resources = dictionaries._parse_manifest(base_dir=tmp_path)

    # Assert
    try:
        assert resources[resource_name].sha256 == checksum
    finally:
        dictionaries._env_checksum_allowlist.cache_clear()


@pytest.mark.unit
def test_parse_manifest__env_allowlist_accepts_unknown_checksum(tmp_path, monkeypatch):
    # Arrange
    resource_name = "custom_resource"
    _write_manifest(
        tmp_path,
        f"""
        version: 1
        resources:
          {resource_name}:
            path: .
            version: "test"
            sha256: "baseline"
            generator: generator.py
        """,
    )
    monkeypatch.setattr(dictionaries, "_compute_sha256", lambda path: "override")
    dictionaries._env_checksum_allowlist.cache_clear()

    # Act & Assert
    with pytest.raises(dictionaries.DictionaryManifestError):
        dictionaries._parse_manifest(base_dir=tmp_path)

    monkeypatch.setenv("CHEMBL_DICTIONARY_CHECKSUM_ALLOWLIST", f"{resource_name}=override")
    dictionaries._env_checksum_allowlist.cache_clear()
    resources = dictionaries._parse_manifest(base_dir=tmp_path)

    # Assert
    try:
        assert resources[resource_name].sha256 == "override"
    finally:
        dictionaries._env_checksum_allowlist.cache_clear()


@pytest.mark.unit
def test_list_resources__loads_bundled_manifest():
    dictionaries._load_manifest.cache_clear()
    try:
        resources = dictionaries.list_resources()
        assert "dictionary_root" in resources
        root_resource = resources["dictionary_root"]
        assert root_resource.path.exists()
        assert root_resource.sha256
    finally:
        dictionaries._load_manifest.cache_clear()
