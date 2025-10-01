from __future__ import annotations

from pathlib import Path

import pytest

from library.config import Config
from library.utils.cli_tools import mapper_batch_main as cli


class DummyLogger:
    """Capture log events emitted via :func:`cli.main`."""

    def __init__(self) -> None:
        self.events: list[tuple[str, dict[str, object]]] = []
        self.config: dict[str, object | None] = {}

    def info(self, event: str, **payload: object) -> None:
        self.events.append((event, payload))


def _configure_logger_factory(loggers: list[DummyLogger]):
    def _configure_logger(
        cfg: object, *, fmt: str | None = None, datefmt: str | None = None
    ) -> DummyLogger:
        if fmt is not None or datefmt is not None:
            raise AssertionError("format overrides should not be supplied")
        logger = DummyLogger()
        logger.config = {"fmt": fmt, "datefmt": datefmt, "cfg": cfg}
        loggers.append(logger)
        return logger

    return _configure_logger


def test_main_invokes_run_without_print_config(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """CLI entry point should configure logging and execute ``run``."""

    configured_loggers: list[DummyLogger] = []
    monkeypatch.setattr(
        cli, "configure_logger", _configure_logger_factory(configured_loggers)
    )

    cfg = Config()
    cfg.api.user_agent = "test@example.org"
    monkeypatch.setattr(cli.cli, "apply_config_overrides", lambda *_, **__: cfg)
    monkeypatch.setattr(cli, "ensure_dirs", lambda *_: None)

    printed: dict[str, Config] = {}
    monkeypatch.setattr(cli, "print_config", lambda c: printed.setdefault("cfg", c))

    called: dict[str, object] = {}

    def fake_run(config: Config, args) -> int:
        called["cfg"] = config
        called["args"] = args
        return 0

    monkeypatch.setattr(cli, "run", fake_run)

    input_csv = tmp_path / "input.csv"
    input_csv.write_text("chembl_id\nCHEMBL1\n", encoding="utf-8")
    output_csv = tmp_path / "out.csv"

    exit_code = cli.main(
        [
            "--input",
            str(input_csv),
            "--output",
            str(output_csv),
            "--column",
            "chembl_id",
            "--chunk-size",
            "2",
        ]
    )

    assert exit_code == 0
    assert called["cfg"] is cfg
    assert called["args"].input_csv == input_csv
    assert "cfg" not in printed
    assert len(configured_loggers) >= 2
    assert any(
        event == "pipeline_start" for event, _ in configured_loggers[0].events
    )
    assert any(
        event == "pipeline_end" and payload.get("exit_code") == 0
        for event, payload in configured_loggers[-1].events
    )
    assert configured_loggers[-1].config["fmt"] is None
    assert configured_loggers[-1].config["datefmt"] is None


def test_main_prints_config_when_requested(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """``--print-config`` should output configuration and skip execution."""

    configured_loggers: list[DummyLogger] = []
    monkeypatch.setattr(
        cli, "configure_logger", _configure_logger_factory(configured_loggers)
    )

    cfg = Config()
    monkeypatch.setattr(cli.cli, "apply_config_overrides", lambda *_, **__: cfg)
    monkeypatch.setattr(cli, "ensure_dirs", lambda *_: None)

    printed: dict[str, Config] = {}
    monkeypatch.setattr(cli, "print_config", lambda c: printed.setdefault("cfg", c))

    def fake_run(*_: object, **__: object) -> int:  # pragma: no cover - safety
        raise AssertionError("run should not be invoked when printing config")

    monkeypatch.setattr(cli, "run", fake_run)

    input_csv = tmp_path / "input.csv"
    input_csv.write_text("chembl_id\nCHEMBL1\n", encoding="utf-8")

    exit_code = cli.main(
        [
            "--input",
            str(input_csv),
            "--column",
            "chembl_id",
            "--chunk-size",
            "2",
            "--print-config",
        ]
    )

    assert exit_code == 0
    assert printed["cfg"] is cfg
    assert len(configured_loggers) >= 2
    assert any(
        event == "pipeline_end" and payload.get("exit_code") == 0
        for event, payload in configured_loggers[-1].events
    )
