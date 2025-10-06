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
def test_parse_manifest__accepts_known_checksum_variants(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
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

    monkeypatch.setattr(
        dictionaries,
        "_compute_sha256",
        lambda path: "efc69f6bb252d68bc7fde11ba98b09b24b0b8fd868fcd6d945eaca76b636f43a",
    )

    resources = dictionaries._parse_manifest(base_dir=manifest_dir)

    assert resources["dictionary_root"].sha256 == "efc69f6bb252d68bc7fde11ba98b09b24b0b8fd868fcd6d945eaca76b636f43a"


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
        lambda path: "efc69f6bb252d68bc7fde11ba98b09b24b0b8fd868fcd6d945eaca76b636f43a",
    )

    resources = dictionaries._parse_manifest(base_dir=manifest_dir)

    assert resources["dictionary_root"].sha256 == "efc69f6bb252d68bc7fde11ba98b09b24b0b8fd868fcd6d945eaca76b636f43a"
def test_manifest_allows_latest_windows_sha256() -> None:
    """The dictionary manifest accepts the hash produced by new Git versions."""

    manifest_path = DICTIONARY_DIR / "manifest.yaml"
    manifest_data = yaml.safe_load(manifest_path.read_text(encoding="utf-8")) or {}
    resources = manifest_data.get("resources", {})
    entry = resources.get("dictionary_root", {})
    sha_values = entry.get("sha256", [])
    if isinstance(sha_values, str):
        sha_values = [sha_values]

    assert (
        "efc69f6bb252d68bc7fde11ba98b09b24b0b8fd868fcd6d945eaca76b636f43a"
        in sha_values
    )
