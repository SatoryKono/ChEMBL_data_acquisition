from __future__ import annotations

import json
from pathlib import Path

import pytest

from library.resources import dictionaries


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
def test_get_resource__accepts_extra_checksum(monkeypatch: pytest.MonkeyPatch) -> None:
    """Extra hashes supplied via the environment authorise new dictionaries."""

    dictionaries._load_manifest.cache_clear()
    custom_hash = "deadbeef" * 8
    monkeypatch.setenv(
        "CHEMBL_DICTIONARY_EXTRA_HASHES",
        json.dumps({"dictionary_root": [custom_hash]}),
    )
    original_compute = dictionaries._compute_sha256

    def _fake_compute(path: Path) -> str:
        if Path(path).resolve() == dictionaries.DICTIONARY_DIR.resolve():
            return custom_hash
        return original_compute(path)

    monkeypatch.setattr(dictionaries, "_compute_sha256", _fake_compute)

    resource = dictionaries.get_resource("dictionary_root")

    assert resource.sha256 == custom_hash

    dictionaries._load_manifest.cache_clear()


@pytest.mark.unit
def test_get_resource__supports_wildcard(monkeypatch: pytest.MonkeyPatch) -> None:
    """The wildcard value disables checksum validation for a resource."""

    dictionaries._load_manifest.cache_clear()
    computed_hash = "cafebabe" * 8
    monkeypatch.setenv(
        "CHEMBL_DICTIONARY_EXTRA_HASHES",
        json.dumps({"dictionary_root": "*"}),
    )
    original_compute = dictionaries._compute_sha256

    def _fake_compute(path: Path) -> str:
        if Path(path).resolve() == dictionaries.DICTIONARY_DIR.resolve():
            return computed_hash
        return original_compute(path)

    monkeypatch.setattr(dictionaries, "_compute_sha256", _fake_compute)

    resource = dictionaries.get_resource("dictionary_root")

    assert resource.sha256 == computed_hash

    dictionaries._load_manifest.cache_clear()


@pytest.mark.unit
def test_list_resources__invalid_extra_hashes(monkeypatch: pytest.MonkeyPatch) -> None:
    """Invalid JSON payloads raise a manifest error immediately."""

    dictionaries._load_manifest.cache_clear()
    monkeypatch.setenv("CHEMBL_DICTIONARY_EXTRA_HASHES", "not json")

    with pytest.raises(dictionaries.DictionaryManifestError):
        dictionaries.list_resources()

    dictionaries._load_manifest.cache_clear()
