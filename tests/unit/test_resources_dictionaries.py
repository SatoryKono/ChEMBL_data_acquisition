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
        "ac67acf2dcd801ffbe9d6e3aa95189af7c3e991fb3ddaaf8aab0be988d7d3224",
        "70f0b19c450d0fc8d19ddb41bd69906d6b1a5ac39e3e4e2d2b6dea54a501569d",
        "95f7a33a028aeeba9027b64f558e50ad25e76934782cc03ba14437fd8eff8476",
        dictionaries.WINDOWS_GIT_248_SPARSE_INDEX_CHECKSUM,
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
        "ac67acf2dcd801ffbe9d6e3aa95189af7c3e991fb3ddaaf8aab0be988d7d3224",
        "70f0b19c450d0fc8d19ddb41bd69906d6b1a5ac39e3e4e2d2b6dea54a501569d",
        "95f7a33a028aeeba9027b64f558e50ad25e76934782cc03ba14437fd8eff8476",
        dictionaries.WINDOWS_GIT_248_SPARSE_INDEX_CHECKSUM,
    }

    assert expected.issubset(set(sha_values))


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
