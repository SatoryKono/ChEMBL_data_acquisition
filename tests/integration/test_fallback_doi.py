from __future__ import annotations

from contextlib import contextmanager
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest
import yaml

from library.config import Config, _serialize_paths
from library.pipelines.document.service import DocumentPipeline
from scripts import get_document_data


class _MemoryLogger:
    def __init__(self) -> None:
        self.events: list[tuple[str, str, dict[str, object]]] = []

    def info(self, event: str, **payload: object) -> None:
        self.events.append(("info", event, dict(payload)))

    def warning(self, event: str, **payload: object) -> None:
        self.events.append(("warning", event, dict(payload)))

    def error(self, event: str, **payload: object) -> None:
        self.events.append(("error", event, dict(payload)))

    def exception(self, event: str, **payload: object) -> None:
        self.events.append(("exception", event, dict(payload)))

    def debug(
        self, event: str, **payload: object
    ) -> None:  # pragma: no cover - optional
        self.events.append(("debug", event, dict(payload)))


@pytest.fixture()
def document_cli_logger(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> _MemoryLogger:
    logger = _MemoryLogger()

    def _configure_logger(cfg: object) -> _MemoryLogger:
        return logger

    @contextmanager
    def _setup_cli_logging(_script: str, log_cfg: object, _date: str | None = None):
        yield SimpleNamespace(log_cfg=log_cfg, log_path=tmp_path / "cli.log")

    monkeypatch.setattr(get_document_data, "logger", logger)
    monkeypatch.setattr(get_document_data, "configure_logger", _configure_logger)
    monkeypatch.setattr(get_document_data, "setup_cli_logging", _setup_cli_logging)
    return logger


def _base_config(tmp_path: Path) -> dict[str, object]:
    cfg = Config()
    output_dir = tmp_path / "output"
    cache_dir = tmp_path / "cache"
    output_dir.mkdir()
    cache_dir.mkdir()
    cfg.io.output_dir = output_dir
    cfg.io.cache_dir = cache_dir
    cfg.io.exist_ok = True
    cfg.system.doc_quality.enable = False
    return _serialize_paths(cfg.to_dict())


def _write_config(tmp_path: Path, data: dict[str, object]) -> Path:
    config_path = tmp_path / "config.yaml"
    config_path.write_text(
        yaml.safe_dump(data, sort_keys=False),
        encoding="utf-8",
    )
    return config_path


def _create_input(tmp_path: Path) -> Path:
    path = tmp_path / "pmids.csv"
    path.write_text("PMID\n111\n222\n", encoding="utf-8")
    return path


@pytest.mark.integration
def test_main_pubmed__without_fallback_csv(
    tmp_path: Path,
    document_cli_logger: _MemoryLogger,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    config = _base_config(tmp_path)
    config_path = _write_config(tmp_path, config)
    input_csv = _create_input(tmp_path)
    output_path = tmp_path / "output" / "pubmed.csv"
    output_path.parent.mkdir(parents=True, exist_ok=True)

    fallback_maps: list[object] = []

    def _fake_fetch(self, pmids, *args, **kwargs):  # type: ignore[override]
        fallback_maps.append(kwargs.get("fallback_doi_map"))
        frame = pd.DataFrame({"document_chembl_id": ["DOC1"], "PubMed.PMID": ["111"]})
        return iter([frame])

    emitted_frames: list[pd.DataFrame] = []

    def _fake_finalise(frames, output, cfg, **kwargs):  # type: ignore[override]
        emitted_frames.extend(list(frames))
        return 0

    monkeypatch.setattr(DocumentPipeline, "fetch_pubmed_records", _fake_fetch)
    monkeypatch.setattr(get_document_data, "_finalise_export", _fake_finalise)
    monkeypatch.setattr(get_document_data, "normalize_documents", lambda frame: frame)

    exit_code = get_document_data.main(
        [
            "pubmed",
            "--config",
            str(config_path),
            "--input",
            str(input_csv),
            "--output-dir",
            str(output_path.parent),
            "--final-out",
            str(output_path),
        ]
    )

    assert exit_code == 0
    assert fallback_maps == [None]
    assert emitted_frames and isinstance(emitted_frames[0], pd.DataFrame)
    start_event = next(
        payload
        for level, event, payload in document_cli_logger.events
        if level == "info" and event == "document_pubmed_start"
    )
    assert start_event["fallback_doi"]["value"] is False
    assert any(
        event == "document_pubmed_done"
        for _level, event, _payload in document_cli_logger.events
    )


@pytest.mark.integration
def test_main_pubmed__with_fallback_csv(
    tmp_path: Path,
    document_cli_logger: _MemoryLogger,
    make_fallback_doi_csv,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    config = _base_config(tmp_path)
    config_path = _write_config(tmp_path, config)
    input_csv = _create_input(tmp_path)
    output_path = tmp_path / "output" / "pubmed.csv"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fallback_csv = make_fallback_doi_csv([("111", "10.1000/fallback")])

    fallback_maps: list[object] = []

    def _fake_fetch(self, pmids, *args, **kwargs):  # type: ignore[override]
        fallback_maps.append(kwargs.get("fallback_doi_map"))
        frame = pd.DataFrame({"document_chembl_id": ["DOC1"], "PubMed.PMID": ["111"]})
        return iter([frame])

    monkeypatch.setattr(DocumentPipeline, "fetch_pubmed_records", _fake_fetch)
    monkeypatch.setattr(
        get_document_data,
        "_finalise_export",
        lambda frames, *_args, **_kwargs: 0,
    )
    monkeypatch.setattr(get_document_data, "normalize_documents", lambda frame: frame)

    exit_code = get_document_data.main(
        [
            "pubmed",
            "--config",
            str(config_path),
            "--input",
            str(input_csv),
            "--output-dir",
            str(output_path.parent),
            "--final-out",
            str(output_path),
            "--fallback-doi-csv",
            str(fallback_csv),
        ]
    )

    assert exit_code == 0
    assert fallback_maps == [{"111": "10.1000/fallback"}]
    start_event = next(
        payload
        for level, event, payload in document_cli_logger.events
        if level == "info" and event == "document_pubmed_start"
    )
    assert start_event["fallback_doi"]["value"] is True


@pytest.mark.integration
@pytest.mark.parametrize("force", [False, True])
def test_main_pubmed__skip_existing_conflict(
    tmp_path: Path,
    document_cli_logger: _MemoryLogger,
    force: bool,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    config = _base_config(tmp_path)
    config_path = _write_config(tmp_path, config)
    input_csv = _create_input(tmp_path)
    output_path = tmp_path / "output" / "pubmed.csv"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text("existing", encoding="utf-8")

    called = False

    def _fake_fetch(self, pmids, *args, **kwargs):  # type: ignore[override]
        nonlocal called
        called = True
        frame = pd.DataFrame({"document_chembl_id": ["DOC1"], "PubMed.PMID": ["111"]})
        return iter([frame])

    monkeypatch.setattr(DocumentPipeline, "fetch_pubmed_records", _fake_fetch)
    monkeypatch.setattr(
        get_document_data,
        "_finalise_export",
        lambda frames, *_args, **_kwargs: 0,
    )
    monkeypatch.setattr(get_document_data, "normalize_documents", lambda frame: frame)

    args = [
        "pubmed",
        "--config",
        str(config_path),
        "--input",
        str(input_csv),
        "--output-dir",
        str(output_path.parent),
        "--final-out",
        str(output_path),
        "--skip-existing",
    ]
    if force:
        args.append("--force")

    exit_code = get_document_data.main(args)

    assert exit_code == 0
    if force:
        assert called is True
    else:
        assert called is False
        assert any(
            event == "pipeline_skip_existing"
            for _level, event, _payload in document_cli_logger.events
        )
