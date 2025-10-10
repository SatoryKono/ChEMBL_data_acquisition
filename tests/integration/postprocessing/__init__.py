from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path

import pandas as pd


def read_csv_header(path: Path) -> list[str]:
    """Return the header of ``path`` without loading the full dataset."""

    frame = pd.read_csv(path, nrows=0, dtype=str, keep_default_na=False)
    return frame.columns.tolist()


def assert_meta_and_header(csv_path: Path, expected_header: Sequence[str]) -> None:
    """Ensure ``csv_path`` has a sidecar and matches ``expected_header``."""

    meta_path = Path(f"{csv_path}.meta.yaml")
    assert meta_path.is_file(), f"Expected metadata sidecar for {csv_path}"
    header = read_csv_header(csv_path)
    assert list(expected_header) == header
