"""Unit tests for :mod:`library.common.sidecar`."""

from __future__ import annotations

import csv

import pytest

from library.common import sidecar


@pytest.mark.unit
def test_sidecar_errors__flushes_buffer_when_chunk_full(tmp_path, monkeypatch):
    """Buffered failures spill to disk once the configured chunk fills up."""

    monkeypatch.setattr(sidecar, "write_meta_yaml", lambda *args, **kwargs: None)

    errors = sidecar.SidecarErrors(chunk_size=2)

    first = {"column": "a", "failure_case": "missing"}
    second = {"column": "b", "failure_case": "invalid"}
    third = {"column": "c", "failure_case": "mismatch"}

    errors.add_error(first)
    assert errors._overflow_path is None  # type: ignore[attr-defined]

    errors.add_error(second)
    overflow_path = errors._overflow_path  # type: ignore[attr-defined]
    assert overflow_path is not None
    assert overflow_path.exists()

    with overflow_path.open("r", encoding="utf8", newline="") as fh:
        flushed_rows = list(csv.DictReader(fh))
    assert [row["failure_case"] for row in flushed_rows] == [
        first["failure_case"],
        second["failure_case"],
    ]

    errors.add_error(third)

    destination = tmp_path / "failures.csv"
    errors.save(destination)

    with destination.open("r", encoding="utf8", newline="") as fh:
        written_rows = list(csv.DictReader(fh))
    assert [row["failure_case"] for row in written_rows] == [
        first["failure_case"],
        second["failure_case"],
        third["failure_case"],
    ]

    assert errors._overflow_path is None  # type: ignore[attr-defined]
