"""Tests for the document type CLI utilities."""

from __future__ import annotations

import csv
import hashlib
from pathlib import Path

import pandas as pd
import pytest
import yaml

from library.config import Config, _serialize_paths
from library.common.run_context import RunContext, set_current
from library.cli.run_context import compute_generated_at
from library.utils.cli_tools import get_document_type


class _StubLogger:
    def __init__(self) -> None:
        self.info_events: list[tuple[str, dict[str, object]]] = []
        self.warning_events: list[tuple[str, dict[str, object]]] = []
        self.debug_events: list[tuple[str, dict[str, object]]] = []
        self.error_events: list[tuple[str, dict[str, object]]] = []

    def info(
        self, event: str, **fields: object
    ) -> None:  # pragma: no cover - simple recorder
        self.info_events.append((event, fields))

    def warning(
        self, event: str, **fields: object
    ) -> None:  # pragma: no cover - simple recorder
        self.warning_events.append((event, fields))

    def debug(
        self, event: str, **fields: object
    ) -> None:  # pragma: no cover - simple recorder
        self.debug_events.append((event, fields))

    def error(
        self, event: str, **fields: object
    ) -> None:  # pragma: no cover - simple recorder
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


@pytest.mark.unit
def test_main__meta_yaml_deterministic_across_runs(
    tmp_path: Path, cfg: Config, monkeypatch: pytest.MonkeyPatch
) -> None:
    set_current(None)
    input_path = tmp_path / "document_types.csv"
    rows = [
        [
            "chembl_id",
            "PubMed.PublicationType",
            "scholar.PublicationTypes",
            "OpenAlex.PublicationTypes",
        ],
        ["CHEMBL1", "Review", "Article", "Article"],
        ["CHEMBL2", "Experimental", "", ""],
    ]
    with input_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
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

    config_path = _write_config(tmp_path / "config.yaml", cfg)

    stub_logger = _StubLogger()

    def _configure_logger_stub(cfg: object) -> _StubLogger:
        run_id = getattr(cfg, "run_id", "")
        level = getattr(cfg, "level", "INFO")
        generated_at = getattr(cfg, "generated_at", "")
        if run_id:
            computed = compute_generated_at(
                date_token=None,
                run_id=run_id,
                seed_parts=("create_logger_config", str(level).upper()),
            )
            try:
                cfg.generated_at = computed  # type: ignore[attr-defined]
            except AttributeError:
                pass
            generated_at = computed
        set_current(RunContext(run_id=run_id, generated_at=generated_at))
        return stub_logger

    monkeypatch.setattr(get_document_type, "configure_logger", _configure_logger_stub)

    output_path = tmp_path / "document_types_output.csv"

    def _run_and_hash() -> str:
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
                "--run-id",
                "unit-doc-type",
            ]
        )

        assert exit_code == 0
        meta_path = output_path.with_name(output_path.name + ".meta.yaml")
        assert meta_path.exists()
        return hashlib.sha256(meta_path.read_bytes()).hexdigest()

    first_hash = _run_and_hash()
    second_hash = _run_and_hash()

    assert first_hash == second_hash
    set_current(None)
