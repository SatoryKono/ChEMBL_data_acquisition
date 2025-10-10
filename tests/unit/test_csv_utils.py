from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from library.common import csv_utils


@pytest.mark.unit
def test_write_csv_deterministic__orders_columns_and_rows(tmp_path: Path) -> None:
    frame = pd.DataFrame(
        {
            "name": ["b", "a"],
            "id": [2, 1],
            "active": [True, False],
            "timestamp": pd.to_datetime(["2024-01-02", "2023-12-31"]),
        }
    )

    out_path = tmp_path / "ordered.csv"
    csv_utils.write_csv_deterministic(
        frame,
        out_path,
        col_order=["id", "name"],
        key_cols=["id"],
    )

    with out_path.open("r", encoding="utf-8-sig") as handle:
        lines = handle.read().splitlines()

    assert lines == [
        "id,name,active,timestamp",
        "1,a,false,2023-12-31",
        "2,b,true,2024-01-02",
    ]


@pytest.mark.unit
def test_write_csv_deterministic__drops_unexpected_columns(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    frame = pd.DataFrame(
        {
            "id": [2, 1],
            "value": ["x", "y"],
            "extra": ["drop", "me"],
        }
    )

    events: list[tuple[str, dict[str, object]]] = []

    def capture_warning(event: str, **fields: object) -> None:
        events.append((event, fields))

    monkeypatch.setattr(csv_utils.logger, "warning", capture_warning)

    out_path = tmp_path / "trimmed.csv"
    csv_utils.write_csv_deterministic(
        frame,
        out_path,
        col_order=["id", "value"],
        key_cols=["id"],
        drop_unexpected_cols=True,
    )

    with out_path.open("r", encoding="utf-8-sig") as handle:
        lines = handle.read().splitlines()

    assert lines == [
        "id,value",
        "1,y",
        "2,x",
    ]
    assert events == [("unexpected_columns_dropped", {"columns": ["extra"]})]
