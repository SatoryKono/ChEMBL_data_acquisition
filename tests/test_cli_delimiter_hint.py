from __future__ import annotations

from pathlib import Path

import pytest

from library.config import IoCfg
from library.io import read_ids


def test_read_ids_tab_separated_hint(tmp_path: Path) -> None:
    """Suggest delimiter when tab-separated file is parsed as comma-separated."""
    path = tmp_path / "ids.tsv"
    path.write_text("id\tname\n1\tfoo\n", encoding="utf-8")
    with pytest.raises(ValueError) as excinfo:
        list(read_ids(path, column="id", cfg=IoCfg()))
    msg = str(excinfo.value)
    assert "--sep '\t'" in msg
    assert "io.csv_sep" in msg
