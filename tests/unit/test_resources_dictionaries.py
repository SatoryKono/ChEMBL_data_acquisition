from __future__ import annotations

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
