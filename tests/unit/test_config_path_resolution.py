from __future__ import annotations

from pathlib import Path

import pytest

from library.config import _absolutise_path_value


@pytest.mark.unit
@pytest.mark.parametrize(
    "value",
    [
        "config/dictionary/_testitem/molecule_hierarchy.csv",
        Path("config/dictionary/_testitem/molecule_hierarchy.csv"),
    ],
)
def test_absolutise_path_value__avoids_duplicate_config_segment(tmp_path, value):
    project_root = tmp_path / "project"
    base_dir = project_root / "config"
    base_dir.mkdir(parents=True)

    resolved = _absolutise_path_value(value, base_dir)

    expected = (
        project_root / "config" / "dictionary" / "_testitem" / "molecule_hierarchy.csv"
    ).resolve()
    if isinstance(resolved, Path):
        assert resolved == expected
    else:
        assert resolved == str(expected)


@pytest.mark.unit
@pytest.mark.parametrize(
    "value",
    [
        "data/output/test.csv",
        Path("data/output/test.csv"),
    ],
)
def test_absolutise_path_value__keeps_standard_relative_paths(tmp_path, value):
    project_root = tmp_path / "workspace"
    base_dir = project_root / "config"
    base_dir.mkdir(parents=True)

    resolved = _absolutise_path_value(value, base_dir)

    expected = (base_dir / value).resolve()
    if isinstance(resolved, Path):
        assert resolved == expected
    else:
        assert resolved == str(expected)
