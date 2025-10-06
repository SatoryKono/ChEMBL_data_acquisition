"""Tests for dictionary resource utilities."""

from __future__ import annotations

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
