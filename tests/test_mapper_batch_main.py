from __future__ import annotations

import argparse
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

    def error(self, event: str, **payload: object) -> None:
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
    assert any(event == "pipeline_start" for event, _ in configured_loggers[0].events)
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


@pytest.mark.parametrize(
    "keep_markers, expected",
    [
        (False, ["CHEMBL1", "CHEMBL3"]),
        (True, ["CHEMBL1", "CUSTOM_NA", "CHEMBL3"]),
    ],
)
def test_run_respects_na_marker_configuration(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    cfg: Config,
    keep_markers: bool,
    expected: list[str],
) -> None:
    cfg.io.na_markers = ("CUSTOM_NA",)
    cfg.io.keep_na_markers = keep_markers

    input_path = tmp_path / "ids.csv"
    input_path.write_text("chembl_id\nCHEMBL1\nCUSTOM_NA\nCHEMBL3\n", encoding="utf-8")
    output_path = tmp_path / "out.csv"

    captured: dict[str, list[str]] = {}

    def fake_map(
        ids: list[str],
        cfg_mapping: object,
        *,
        batch_size: int,
        rps: float,
        max_workers: int | None,
    ) -> dict[str, str | None]:
        captured["ids"] = list(ids)
        return {chembl_id: None for chembl_id in ids}

    monkeypatch.setattr(cli, "map_chembl_ids_to_uniprot", fake_map)

    written: dict[str, object] = {}

    def fake_write_csv(
        df_out,
        path: Path,
        *,
        cfg: Config,
        sep: str,
        encoding: str,
    ) -> Path:
        written["df"] = df_out.copy()
        written["path"] = path
        return path

    monkeypatch.setattr(cli.io, "write_csv", fake_write_csv)

    args = argparse.Namespace(
        input_csv=input_path,
        output_csv=output_path,
        column="chembl_id",
        sep=",",
        encoding="utf8",
        chunk_size=2,
        rps=5.0,
        workers=2,
    )

    exit_code = cli.run(cfg, args)
    assert exit_code == 0
    assert captured["ids"] == expected
    assert written["path"] == output_path
