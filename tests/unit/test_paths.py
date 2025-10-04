from __future__ import annotations

from pathlib import Path

import pytest

from config.paths import DICTIONARY_DIR


@pytest.mark.unit
def test_dictionary_dir__points_to_project_resource() -> None:
    """DICTIONARY_DIR should resolve to ``config/dictionary`` within the repo."""

    project_root = Path(__file__).resolve().parents[2]
    expected = project_root / "config" / "dictionary"

    assert DICTIONARY_DIR == expected
    assert DICTIONARY_DIR.is_absolute()
    assert DICTIONARY_DIR.relative_to(project_root) == Path("config/dictionary")
