from __future__ import annotations

import csv

import pytest

from library.common.sidecar import SidecarErrors


@pytest.mark.unit
def test_sidecar_errors__flushes_buffer_when_chunk_full(tmp_path) -> None:
    errors = SidecarErrors(chunk_size=2)

    first = {"index": 1, "column": "value-1"}
    second = {"index": 2, "column": "value-2"}
    third = {"index": 3, "column": "value-3"}

    errors.add_error(first)
    assert errors._overflow_path is None  # type: ignore[attr-defined]

    errors.add_error(second)

    overflow_path = errors._overflow_path  # type: ignore[attr-defined]
    assert overflow_path is not None
    assert overflow_path.exists()
    assert errors._errors == []  # type: ignore[attr-defined]

    with overflow_path.open("r", encoding="utf8", newline="") as fh:
        reader = csv.DictReader(fh)
        rows = list(reader)

    assert [row["index"] for row in rows] == ["1", "2"]

    errors.add_error(third)

    assert len(errors._errors) == 1  # type: ignore[attr-defined]

    failure_path = tmp_path / "failure_cases.csv"
    errors.save(failure_path)

    with failure_path.open("r", encoding="utf8", newline="") as fh:
        reader = csv.DictReader(fh)
        saved_rows = list(reader)

    assert [row["index"] for row in saved_rows] == ["1", "2", "3"]
    assert errors._overflow_path is None  # type: ignore[attr-defined]

