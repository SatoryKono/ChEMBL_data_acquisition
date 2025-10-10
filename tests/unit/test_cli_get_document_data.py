from __future__ import annotations

import argparse
from collections.abc import Callable
from contextlib import contextmanager
from pathlib import Path
from types import SimpleNamespace

import pytest
import yaml

from library.config import Config, _serialize_paths
from scripts import get_document_data


class _RunRecorder:
    def __init__(self) -> None:
        self.called = False
        self.config = None
        self.args: argparse.Namespace | None = None

    def __call__(self, cfg, args: argparse.Namespace) -> int:  # type: ignore[override]
        self.called = True
        self.config = cfg
        self.args = args
        return 0


def _prepare_config(tmp_path: Path, mutator: Callable[[Config], None]) -> Path:
    cfg = Config()
    output_dir = tmp_path / "output"
    cache_dir = tmp_path / "cache"
    output_dir.mkdir()
    cache_dir.mkdir()
    cfg.io.output_dir = output_dir
    cfg.io.cache_dir = cache_dir
    cfg.io.exist_ok = True
    cfg.system.doc_quality.enable = False
    mutator(cfg)
    data = _serialize_paths(cfg.to_dict())
    config_path = tmp_path / "config.yaml"
    config_path.write_text(
        yaml.safe_dump(data, sort_keys=False),
        encoding="utf-8",
    )
    return config_path


@contextmanager
def _stub_logging(tmp_path: Path) -> Callable[[object], object]:
    logger_stub = get_document_data.logger

    def _configure_logger(cfg: object) -> object:
        return logger_stub

    @contextmanager
    def _setup_cli_logging(_name: str, log_cfg: object, _date: str | None = None):
        yield SimpleNamespace(log_cfg=log_cfg, log_path=tmp_path / "cli.log")

    original_configure = get_document_data.configure_logger
    original_setup = get_document_data.setup_cli_logging
    get_document_data.configure_logger = _configure_logger  # type: ignore[assignment]
    get_document_data.setup_cli_logging = _setup_cli_logging  # type: ignore[assignment]
    try:
        yield lambda _: logger_stub
    finally:
        get_document_data.configure_logger = original_configure  # type: ignore[assignment]
        get_document_data.setup_cli_logging = original_setup  # type: ignore[assignment]


@pytest.mark.unit
def test_build_parser__pubmed_defaults() -> None:
    parser, _log_cfg = get_document_data.build_parser()
    args = parser.parse_args(["pubmed"])

    assert args.command == "pubmed"
    assert args.batch_size == 100
    assert args.limit is None
    assert args.fallback_doi_csv is None
    assert args.offset == 0
    assert args.timeout is None


@pytest.mark.unit
def test_main__cli_overrides_config(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    def _mutator(cfg: Config) -> None:
        cfg.document.pubmed.batch_size = 33

    config_path = _prepare_config(tmp_path, _mutator)
    input_csv = tmp_path / "pmids.csv"
    input_csv.write_text("PMID\n12345\n", encoding="utf-8")
    recorder = _RunRecorder()
    monkeypatch.setattr(get_document_data, "run", recorder)

    with _stub_logging(tmp_path):
        exit_code = get_document_data.main(
            [
                "pubmed",
                "--config",
                str(config_path),
                "--input",
                str(input_csv),
                "--final-out",
                str(tmp_path / "pubmed.csv"),
                "--batch-size",
                "77",
            ]
        )

    assert exit_code == 0
    assert recorder.called is True
    assert recorder.config.document.pubmed.batch_size == 77
    assert recorder.args is not None
    assert recorder.args.batch_size == 77


@pytest.mark.unit
def test_main__config_overrides_defaults(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    def _mutator(cfg: Config) -> None:
        cfg.document.pubmed.batch_size = 31

    config_path = _prepare_config(tmp_path, _mutator)
    input_csv = tmp_path / "pmids.csv"
    input_csv.write_text("PMID\n67890\n", encoding="utf-8")
    recorder = _RunRecorder()
    monkeypatch.setattr(get_document_data, "run", recorder)

    with _stub_logging(tmp_path):
        exit_code = get_document_data.main(
            [
                "pubmed",
                "--config",
                str(config_path),
                "--input",
                str(input_csv),
                "--final-out",
                str(tmp_path / "pubmed.csv"),
            ]
        )

    assert exit_code == 0
    assert recorder.called is True
    assert recorder.config.document.pubmed.batch_size == 31
    assert recorder.args is not None
    assert recorder.args.batch_size == 31


@pytest.mark.unit
@pytest.mark.parametrize(
    "option, value, message",
    [
        ("--limit", "-1", "--limit must be zero or a positive integer"),
        ("--offset", "-5", "--offset must be zero or a positive integer"),
    ],
)
def test_main__negative_values_raise(
    option: str,
    value: str,
    message: str,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    called = False

    def _fail(*_args: object, **_kwargs: object) -> int:
        nonlocal called
        called = True
        return 0

    monkeypatch.setattr(get_document_data, "run_cli_command", _fail)

    with pytest.raises(SystemExit) as excinfo:
        get_document_data.main(["pubmed", option, value])

    assert excinfo.value.code == 2
    assert called is False
    captured = capsys.readouterr()
    assert message in captured.err
