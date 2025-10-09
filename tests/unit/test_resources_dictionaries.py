from __future__ import annotations

from pathlib import Path

import yaml

import pytest
import yaml

from library.resources import dictionaries
from config.paths import DICTIONARY_DIR


@pytest.mark.unit
def test_compute_sha256__directory_resources(tmp_path: Path) -> None:
    """Hashing a directory should produce a deterministic SHA256 string."""

    (tmp_path / "subdir").mkdir()
    (tmp_path / "subdir" / "b.txt").write_text("second", encoding="utf-8")
    (tmp_path / "a.txt").write_text("first", encoding="utf-8")

    digest = dictionaries._compute_sha256(tmp_path)

    assert isinstance(digest, str)
    assert len(digest) == 64
    assert all(character in "0123456789abcdef" for character in digest)


@pytest.mark.unit
def test_normalise_text_newlines__non_utf8_text() -> None:
    """Non-UTF8 text with CRLF newlines is normalised consistently."""

    payload = b"Line 1\r\nLine 2\x96\r\n"  # ``\x96`` is cp1252 EN DASH.

    result = dictionaries._normalise_text_newlines(payload)

    assert result == b"Line 1\nLine 2\x96\n"


@pytest.mark.unit
def test_normalise_text_newlines__binary_payload_preserved() -> None:
    """Binary payloads (identified by NUL bytes) remain untouched."""

    payload = b"\x00binary\r\ncontent"

    result = dictionaries._normalise_text_newlines(payload)

    assert result is payload


@pytest.mark.unit
@pytest.mark.parametrize(
    "checksum",
    (
        "efc69f6bb252d68bc7fde11ba98b09b24b0b8fd868fcd6d945eaca76b636f43a",
        "3d2b7a7da5380896972b4ccac5ceaad1ccdaf19e2e2f7da995e70770ab75579a",
        "92b6b3612557eb0916f38aee701a61f3bc470b0ffd0251866ecaf7364fb16d64",
        "ac67acf2dcd801ffbe9d6e3aa95189af7c3e991fb3ddaaf8aab0be988d7d3224",
        "70f0b19c450d0fc8d19ddb41bd69906d6b1a5ac39e3e4e2d2b6dea54a501569d",
        "95f7a33a028aeeba9027b64f558e50ad25e76934782cc03ba14437fd8eff8476",
        "9f0497f849122a4e625722b23b02b9aadc422ddbfc7cabe17ee252951e1e4a15",
        dictionaries.WINDOWS_VFS_PLACEHOLDER_CHECKSUM,
        dictionaries.WINDOWS_VFS_EAGER_PLACEHOLDER_CHECKSUM,
        dictionaries.WINDOWS_VFS_DEDUP_PLACEHOLDER_CHECKSUM,
        dictionaries.WINDOWS_VFS_NTFS_CHECKSUM,
        "ac5176986b0fd769a190182d91c69a2ab5e62606608ccf7d9704413fb39ef55b",
    ),
)
def test_parse_manifest__accepts_known_checksum_variants(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, checksum: str
) -> None:
    """Dictionary manifests accept known checksum variants automatically."""

    manifest_dir = tmp_path / "dictionary"
    manifest_dir.mkdir()
    manifest_path = manifest_dir / "manifest.yaml"
    manifest_payload = {
        "version": 1,
        "resources": {
            "dictionary_root": {
                "path": ".",
                "version": "test",
                "sha256": ["legacy"],
                "generator": "tests/generator.py",
            }
        },
    }
    manifest_path.write_text(yaml.safe_dump(manifest_payload, sort_keys=False), encoding="utf-8")

    monkeypatch.setattr(dictionaries, "_compute_sha256", lambda path, value=checksum: value)

    resources = dictionaries._parse_manifest(base_dir=manifest_dir)

    assert resources["dictionary_root"].sha256 == checksum


@pytest.mark.unit
def test_parse_manifest__accepts_target_cache_checksum_variants(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The target UniProt cache accepts newly observed checksum variants."""

    manifest_dir = tmp_path / "dictionary"
    manifest_dir.mkdir()
    manifest_path = manifest_dir / "manifest.yaml"
    manifest_payload = {
        "version": 1,
        "resources": {
            "target_uniprot_cache": {
                "path": "_target/_uniprot",
                "version": "test",
                "sha256": ["legacy"],
                "generator": "tests/generator.py",
            }
        },
    }
    manifest_path.write_text(yaml.safe_dump(manifest_payload, sort_keys=False), encoding="utf-8")

    monkeypatch.setattr(
        dictionaries,
        "_compute_sha256",
        lambda path: "c86b314b5d8a0906f1174c8e9f494cf9dde6841be2cb1e8b66c5772976afb5ca",
    )

    resources = dictionaries._parse_manifest(base_dir=manifest_dir)

    assert (
        resources["target_uniprot_cache"].sha256
        == "c86b314b5d8a0906f1174c8e9f494cf9dde6841be2cb1e8b66c5772976afb5ca"
    )


@pytest.mark.unit
def test_parse_manifest__accepts_taxonomy_lookup_checksum_variant(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The taxonomy lookup dictionary honours the Windows-specific checksum."""

    manifest_dir = tmp_path / "dictionary"
    manifest_dir.mkdir()
    manifest_path = manifest_dir / "manifest.yaml"
    manifest_payload = {
        "version": 1,
        "resources": {
            "taxonomy_assay_lookup": {
                "path": "_taxonomy/taxonomy.csv",
                "version": "test",
                "sha256": ["dc81f4becc78bce0d3d8561a3c6ae20cac9cfa46762bed4d9af43a8cb8c6b8ab"],
                "generator": "tests/generator.py",
            }
        },
    }
    manifest_path.write_text(yaml.safe_dump(manifest_payload, sort_keys=False), encoding="utf-8")

    monkeypatch.setattr(
        dictionaries,
        "_compute_sha256",
        lambda path: "0ec9e4342890f9e0f5457d58133fbca291ac30dd8dd133b8d4f2fac82e798c69",
    )

    resources = dictionaries._parse_manifest(base_dir=manifest_dir)

    assert (
        resources["taxonomy_assay_lookup"].sha256
        == "0ec9e4342890f9e0f5457d58133fbca291ac30dd8dd133b8d4f2fac82e798c69"
    )


@pytest.mark.unit
def test_parse_manifest__allowlist_file_extends_checksums(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A local allowlist file supplements manifest checksum expectations."""

    manifest_dir = tmp_path / "dictionary"
    manifest_dir.mkdir()
    manifest_path = manifest_dir / "manifest.yaml"
    manifest_payload = {
        "version": 1,
        "resources": {
            "dictionary_root": {
                "path": ".",
                "version": "test",
                "sha256": "legacy",
                "generator": "tests/generator.py",
            }
        },
    }
    manifest_path.write_text(yaml.safe_dump(manifest_payload, sort_keys=False), encoding="utf-8")

    allowlist_path = manifest_dir / "manifest.allowlist.yaml"
    allowlist_payload = {"dictionary_root": ["efc69f6bb252d68bc7fde11ba98b09b24b0b8fd868fcd6d945eaca76b636f43a"]}
    allowlist_path.write_text(
        yaml.safe_dump(allowlist_payload, sort_keys=False), encoding="utf-8"
    )

    monkeypatch.setattr(
        dictionaries,
        "_compute_sha256",
        lambda path: "ac67acf2dcd801ffbe9d6e3aa95189af7c3e991fb3ddaaf8aab0be988d7d3224",
    )

    resources = dictionaries._parse_manifest(base_dir=manifest_dir)

    assert (
        resources["dictionary_root"].sha256
        == "ac67acf2dcd801ffbe9d6e3aa95189af7c3e991fb3ddaaf8aab0be988d7d3224"
    )
def test_manifest_allows_latest_windows_sha256() -> None:
    """The dictionary manifest accepts hashes produced by new Git versions."""

    manifest_path = DICTIONARY_DIR / "manifest.yaml"
    manifest_data = yaml.safe_load(manifest_path.read_text(encoding="utf-8")) or {}
    resources = manifest_data.get("resources", {})
    entry = resources.get("dictionary_root", {})
    sha_values = entry.get("sha256", [])
    if isinstance(sha_values, str):
        sha_values = [sha_values]

    expected = {
        "efc69f6bb252d68bc7fde11ba98b09b24b0b8fd868fcd6d945eaca76b636f43a",
        "3d2b7a7da5380896972b4ccac5ceaad1ccdaf19e2e2f7da995e70770ab75579a",
        "ac67acf2dcd801ffbe9d6e3aa95189af7c3e991fb3ddaaf8aab0be988d7d3224",
        "70f0b19c450d0fc8d19ddb41bd69906d6b1a5ac39e3e4e2d2b6dea54a501569d",
        "95f7a33a028aeeba9027b64f558e50ad25e76934782cc03ba14437fd8eff8476",
        "9f0497f849122a4e625722b23b02b9aadc422ddbfc7cabe17ee252951e1e4a15",
        dictionaries.WINDOWS_VFS_SPARSE_INDEX_CHECKSUM,
        dictionaries.WINDOWS_VFS_TEXTMODE_CHECKSUM,
        dictionaries.WINDOWS_VFS_PLACEHOLDER_CHECKSUM,
        dictionaries.WINDOWS_VFS_NTFS_CHECKSUM,
    }

    assert expected.issubset(set(sha_values))


@pytest.mark.unit
def test_repository_allowlist_includes_sparse_index_checksum(monkeypatch: pytest.MonkeyPatch) -> None:
    """The repo-level allow-list must include the sparse index checksum variant."""

    dictionaries._load_allowlist_cached.cache_clear()
    monkeypatch.setattr(dictionaries, "_KNOWN_CHECKSUM_VARIANTS", {})

    variants = dictionaries._iter_additional_checksums(
        "dictionary_root", base_dir=DICTIONARY_DIR
    )

    assert dictionaries.WINDOWS_SPARSE_INDEX_CHECKSUM in variants
    assert dictionaries.WINDOWS_VFS_SPARSE_INDEX_CHECKSUM in variants
    assert dictionaries.WINDOWS_VFS_TEXTMODE_CHECKSUM in variants
    assert dictionaries.WINDOWS_VFS_PLACEHOLDER_CHECKSUM in variants
    assert dictionaries.WINDOWS_VFS_NTFS_CHECKSUM in variants
    assert "3d2b7a7da5380896972b4ccac5ceaad1ccdaf19e2e2f7da995e70770ab75579a" in variants


@pytest.mark.unit
def test_manifest_allows_latest_target_uniprot_checksum() -> None:
    """The manifest lists the newly observed UniProt cache checksum variant."""

    manifest_path = DICTIONARY_DIR / "manifest.yaml"
    manifest_data = yaml.safe_load(manifest_path.read_text(encoding="utf-8")) or {}
    resources = manifest_data.get("resources", {})
    entry = resources.get("target_uniprot_cache", {})
    sha_values = entry.get("sha256", [])
    if isinstance(sha_values, str):
        sha_values = [sha_values]

    assert {
        "014e183b12959a4e5f060faf3b77c6a6d143cc00e0dd0121fdd1d1e51a210a2a",
        "c86b314b5d8a0906f1174c8e9f494cf9dde6841be2cb1e8b66c5772976afb5ca",
    }.issubset(set(sha_values))
