"""Tests for the document type CLI utilities."""

from __future__ import annotations

import csv
from pathlib import Path

import pandas as pd
import pytest
import yaml

from library.config import Config, _serialize_paths
from library.utils.cli_tools import get_document_type


class _StubLogger:
    def __init__(self) -> None:
        self.info_events: list[tuple[str, dict[str, object]]] = []
        self.warning_events: list[tuple[str, dict[str, object]]] = []
        self.debug_events: list[tuple[str, dict[str, object]]] = []
        self.error_events: list[tuple[str, dict[str, object]]] = []

    def info(self, event: str, **fields: object) -> None:  # pragma: no cover - simple recorder
        self.info_events.append((event, fields))

    def warning(self, event: str, **fields: object) -> None:  # pragma: no cover - simple recorder
        self.warning_events.append((event, fields))

    def debug(self, event: str, **fields: object) -> None:  # pragma: no cover - simple recorder
        self.debug_events.append((event, fields))

    def error(self, event: str, **fields: object) -> None:  # pragma: no cover - simple recorder
        self.error_events.append((event, fields))


def _write_config(path: Path, cfg: Config) -> Path:
    data = _serialize_paths(cfg.to_dict())
    path.write_text(yaml.safe_dump(data, sort_keys=False), encoding="utf-8")
    return path


@pytest.mark.unit
def test_main__fallback_encoding_and_separator(
    tmp_path: Path, cfg: Config, monkeypatch: pytest.MonkeyPatch
) -> None:
    input_path = tmp_path / "document_types.csv"
    rows = [
        [
            "chembl_id",
            "PubMed.PublicationType",
            "scholar.PublicationTypes",
            "OpenAlex.PublicationTypes",
        ],
        ["CHEMBL1", "Reseña;Ensayo", "Artículo Ñ", ""],
        ["CHEMBL2", "", "", ""],
    ]
    with input_path.open("w", newline="", encoding="latin-1") as handle:
        writer = csv.writer(handle, delimiter=";")
        writer.writerows(rows)

    output_dir = tmp_path / "out"
    cache_dir = tmp_path / "cache"
    output_dir.mkdir()
    cache_dir.mkdir()
    cfg.io.output_dir = output_dir
    cfg.io.cache_dir = cache_dir
    cfg.io.exist_ok = True
    cfg.io.csv_encoding = "utf-8"
    cfg.io.csv_sep = ","
    cfg.io.csv_fallback_encodings = ["latin-1"]
    cfg.io.csv_fallback_separators = [";"]
    cfg.doc_type.limit = 1

    config_path = _write_config(tmp_path / "config.yaml", cfg)

    stub_logger = _StubLogger()
    monkeypatch.setattr(get_document_type, "configure_logger", lambda _: stub_logger)

    output_path = tmp_path / "classified.csv"
    exit_code = get_document_type.main(
        [
            "--config",
            str(config_path),
            "--input",
            str(input_path),
            "--output",
            str(output_path),
            "--encoding",
            "utf-8",
            "--sep",
            ",",
        ]
    )

    assert exit_code == 0
    frame = pd.read_csv(output_path)
    assert list(frame["chembl_id"]) == ["CHEMBL1"]
    assert "class_label" in frame.columns

    info_events = {event for event, _fields in stub_logger.info_events}
    assert "csv_encoding_fallback_used" in info_events
    assert "csv_separator_fallback_used" in info_events
