from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path

import pandas as pd


def read_csv_header(path: Path) -> list[str]:
    """Return the header of ``path`` without loading the full dataset."""

    frame = pd.read_csv(path, nrows=0, dtype=str, keep_default_na=False)
    return frame.columns.tolist()


def assert_csv_header(csv_path: Path, expected_header: Sequence[str]) -> None:
    """Ensure ``csv_path`` matches ``expected_header``."""

    header = read_csv_header(csv_path)
    assert list(expected_header) == header
