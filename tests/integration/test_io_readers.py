from __future__ import annotations

import csv
from pathlib import Path

import pytest

import library.io.readers as io_readers
from library.io import read_ids


@pytest.mark.integration
def test_read_ids__fallback_encoding_recovers_values(
    tmp_path: Path, cfg, monkeypatch: pytest.MonkeyPatch
) -> None:
    path = tmp_path / "ids_latin.csv"
    rows = [
        ["chembl_id"],
        ["CHEMBL1"],
        ["CHEMBLÑ"],
    ]
    with path.open("w", newline="", encoding="latin-1") as handle:
        writer = csv.writer(handle, delimiter=",")
        writer.writerows(rows)

    cfg.io.csv_encoding = "utf-8"
    cfg.io.csv_fallback_encodings = ["latin-1"]

    events: list[tuple[str, dict[str, object]]] = []

    def capture_warning(event: str, **fields: object) -> None:
        events.append((event, fields))

    monkeypatch.setattr(io_readers.logger, "warning", capture_warning)

    identifiers = list(read_ids(path, column="chembl_id", cfg=cfg.io))

    assert identifiers == ["CHEMBL1", "CHEMBLÑ"]
    assert any(event == "csv_decode_failed" for event, _ in events)


@pytest.mark.integration
def test_read_ids__fallback_separator_emits_info(
    tmp_path: Path, cfg, monkeypatch: pytest.MonkeyPatch
) -> None:
    path = tmp_path / "ids_semicolon.csv"
    path.write_text(
        "chembl_id;other\nCHEMBL10;foo\nCHEMBL11;bar\n",
        encoding="utf-8",
    )

    cfg.io.csv_sep = ","
    cfg.io.csv_fallback_separators = [";"]

    events: list[tuple[str, dict[str, object]]] = []

    def capture_info(event: str, **fields: object) -> None:
        events.append((event, fields))

    monkeypatch.setattr(io_readers.logger, "info", capture_info)

    identifiers = list(read_ids(path, column="chembl_id", cfg=cfg.io))

    assert identifiers == ["CHEMBL10", "CHEMBL11"]
    assert any(event == "csv_separator_fallback_used" for event, _ in events)


@pytest.mark.integration
def test_read_ids__empty_file_returns_no_values(tmp_path: Path, cfg) -> None:
    path = tmp_path / "ids_empty.csv"
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["chembl_id"])

    cfg.io.csv_encoding = "utf-8"

    identifiers = list(read_ids(path, column="chembl_id", cfg=cfg.io))

    assert identifiers == []
