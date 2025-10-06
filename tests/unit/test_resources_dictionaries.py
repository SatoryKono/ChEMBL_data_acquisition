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
