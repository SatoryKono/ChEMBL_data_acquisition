from __future__ import annotations

from library.io.path_utils import OUTPUT_DIR, ROOT


def test_root_points_to_project_root() -> None:
    assert (ROOT / "pyproject.toml").is_file()


def test_output_dir_is_under_data_output() -> None:
    expected = ROOT / "data" / "output"
    assert OUTPUT_DIR == expected
    assert OUTPUT_DIR.is_dir()
    assert OUTPUT_DIR.samefile(expected)
