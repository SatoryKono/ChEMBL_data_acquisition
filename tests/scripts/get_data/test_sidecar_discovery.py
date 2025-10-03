"""Regression tests for sidecar discovery pruning."""

from __future__ import annotations

from pathlib import Path

import pytest

from scripts import get_data


def test_discover_sidecars_skips_irrelevant_directories(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Large unrelated directories should not be traversed during discovery."""

    final_dir = tmp_path / "outputs"
    final_dir.mkdir()
    final_output = final_dir / "output.activities_20240101.csv"
    final_output.touch()
    working_output = final_dir / ".output.activities_20240101.csv.tmp"
    working_output.touch()

    quality_dir = final_dir / f"{final_output.stem}_quality"
    reports_dir = quality_dir / "reports"
    reports_dir.mkdir(parents=True)
    quality_file = reports_dir / f"{final_output.name}.quality.json"
    quality_file.write_text("{}", encoding="utf-8")

    unrelated_root = final_dir / "unrelated"
    for idx in range(5):
        branch = unrelated_root / f"branch_{idx}"
        nested = branch / "payload"
        nested.mkdir(parents=True, exist_ok=True)
        for part in range(10):
            artifact = nested / f"artifact_{idx}_{part}.txt"
            artifact.write_text("payload", encoding="utf-8")

    path_type = type(final_dir)
    original_iterdir = path_type.iterdir

    def guarded_iterdir(self: Path):  # type: ignore[override]
        if self == unrelated_root:
            raise AssertionError("discover_sidecars should not scan unrelated trees")
        return original_iterdir(self)

    monkeypatch.setattr(path_type, "iterdir", guarded_iterdir)

    sidecars = get_data._discover_sidecars(final_output, working_output)

    rel_quality = quality_file.relative_to(final_dir)
    assert rel_quality in sidecars
    assert sidecars[rel_quality].destination == final_dir / rel_quality


def test_discover_sidecars_respects_depth_limit(tmp_path: Path) -> None:
    """Depth limits restrict traversal to the configured hierarchy."""

    final_dir = tmp_path / "outputs"
    final_dir.mkdir()
    final_output = final_dir / "output.targets_20240101.csv"
    final_output.touch()
    working_output = final_dir / ".output.targets_20240101.csv.tmp"
    working_output.touch()

    nested_dir = final_dir / f"{final_output.stem}_quality" / "nested" / "leaf"
    nested_dir.mkdir(parents=True)
    nested_file = nested_dir / f"{final_output.name}.quality.json"
    nested_file.write_text("{}", encoding="utf-8")

    shallow_sidecars = get_data._discover_sidecars(final_output, working_output, max_depth=1)
    nested_rel = nested_file.relative_to(final_dir)
    assert nested_rel not in shallow_sidecars

    deep_sidecars = get_data._discover_sidecars(final_output, working_output, max_depth=5)
    assert nested_rel in deep_sidecars
