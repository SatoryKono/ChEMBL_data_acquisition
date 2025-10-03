"""Regression tests for ``csv_utils_main`` configuration handling."""

from __future__ import annotations

import textwrap
from pathlib import Path

import yaml

from library.cli import LoggerConfig

from library.utils.cli_tools import csv_utils_main as cli


def test_csv_utils_cli_print_config_exits_without_writing(
    tmp_path: Path, capsys
) -> None:
    """Printing the configuration exits early and avoids writing CSV output."""

    input_csv = tmp_path / "input.csv"
    input_csv.write_text("a,b\n1,2\n", encoding="utf-8")
    output_csv = tmp_path / "out.csv"
    config_path = tmp_path / "config.yaml"
    config_path.write_text(
        textwrap.dedent(
            """
            system:
              log:
                level: "ERROR"
            local:
              io:
                csv_sep: ","
                csv_encoding: "utf-8"
                exist_ok: true
            """
        ),
        encoding="utf-8",
    )

    exit_code = cli.main(
        [
            "--config",
            str(config_path),
            "--input",
            str(input_csv),
            "--output",
            str(output_csv),
            "--log-level",
            "DEBUG",
            "--sep",
            "|",
            "--print-config",
        ]
    )

    assert exit_code == 0
    captured = capsys.readouterr()
    data = yaml.safe_load(captured.out)
    assert data["system"]["log"]["level"] == "DEBUG"
    assert data["local"]["io"]["csv_sep"] == "|"
    assert not output_csv.exists()


def test_csv_utils_cli_print_config_uses_logger_config_run_id(monkeypatch) -> None:
    """``--print-config`` configures logging with the generated run identifier."""

    created_config = LoggerConfig(run_id="sentinel", level="INFO")
    configured: list[LoggerConfig] = []

    monkeypatch.setattr(cli, "create_logger_config", lambda level: created_config)
    monkeypatch.setattr(cli, "configure_logger", lambda cfg: configured.append(cfg))
    monkeypatch.setattr(cli, "print_config", lambda cfg: None)

    exit_code = cli.main(["--print-config"])

    assert exit_code == 0
    assert configured == [created_config]
