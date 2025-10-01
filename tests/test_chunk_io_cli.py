"""CLI smoke tests for chunk_io_main."""

from __future__ import annotations

from pathlib import Path
import textwrap

import pandas as pd
import yaml

from library.utils.cli_tools import chunk_io_main as cli


def test_cli_limit(tmp_path: Path) -> None:
    """CLI processes only the requested number of rows."""
    input_path = tmp_path / "input.csv"
    df = pd.DataFrame({"a": range(200)})
    df.to_csv(input_path, index=False)
    output_path = tmp_path / "out.csv"
    checkpoint = tmp_path / "cp.json"
    cli.main(
        [
            "--input",
            str(input_path),
            "--output",
            str(output_path),
            "--chunk-size",
            "50",
            "--limit",
            "100",
            "--checkpoint",
            str(checkpoint),
            "--log-level",
            "WARNING",
        ]
    )
    result = pd.read_csv(output_path)
    assert len(result) == 100


def test_cli_creates_output_directory(tmp_path: Path) -> None:
    """CLI creates missing output directories when allowed by configuration."""
    input_path = tmp_path / "input.csv"
    df = pd.DataFrame({"a": range(10)})
    df.to_csv(input_path, index=False)
    output_path = tmp_path / "missing" / "sub" / "out.csv"
    checkpoint = tmp_path / "cp.json"
    assert not output_path.parent.exists()

    exit_code = cli.main(
        [
            "--input",
            str(input_path),
            "--output",
            str(output_path),
            "--chunk-size",
            "5",
            "--checkpoint",
            str(checkpoint),
            "--log-level",
            "WARNING",
        ]
    )

    assert exit_code == 0
    assert output_path.exists()
    result = pd.read_csv(output_path)
    assert len(result) == len(df)


def test_cli_print_config_honours_overrides(tmp_path: Path, capsys) -> None:
    """CLI prints the effective configuration with overrides applied."""

    input_path = tmp_path / "input.csv"
    input_path.write_text("a\n1\n", encoding="utf-8")
    config_path = tmp_path / "config.yaml"
    config_path.write_text(
        textwrap.dedent(
            """
            system:
              log:
                level: "ERROR"
            local:
              io:
                csv_sep: ";"
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
            str(input_path),
            "--output",
            str(tmp_path / "out.csv"),
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


def test_cli_invalid_config_logs_error(tmp_path: Path, monkeypatch) -> None:
    """CLI reports invalid configurations and exits with code 1."""

    input_path = tmp_path / "input.csv"
    input_path.write_text("a\n1\n", encoding="utf-8")
    config_path = tmp_path / "config.yaml"
    config_path.write_text(
        textwrap.dedent(
            """
            system:
              log:
                level: "NOT_A_LEVEL"
            """
        ),
        encoding="utf-8",
    )

    logged: list[tuple[str, dict[str, object]]] = []

    def fake_error(event: str, *args: object, **kwargs: object) -> None:
        logged.append((event, kwargs))

    monkeypatch.setattr(cli.logger, "error", fake_error)

    exit_code = cli.main(
        [
            "--config",
            str(config_path),
            "--input",
            str(input_path),
            "--output",
            str(tmp_path / "out.csv"),
            "--log-level",
            "INFO",
        ]
    )

    assert exit_code == 1
    assert logged
    event, payload = logged[0]
    assert event == "config_error"
    assert payload.get("config") == str(config_path)


def test_cli_invalid_chunk_size_logs_error(tmp_path: Path, monkeypatch) -> None:
    """CLI rejects non-positive chunk sizes with an error message."""

    input_path = tmp_path / "input.csv"
    input_path.write_text("a\n1\n", encoding="utf-8")
    output_path = tmp_path / "out.csv"

    logged: list[tuple[str, dict[str, object]]] = []

    def fake_error(event: str, *args: object, **kwargs: object) -> None:
        logged.append((event, kwargs))

    monkeypatch.setattr(cli.logger, "error", fake_error)

    exit_code = cli.main(
        [
            "--input",
            str(input_path),
            "--output",
            str(output_path),
            "--chunk-size",
            "0",
            "--log-level",
            "WARNING",
        ]
    )

    assert exit_code == 1
    event, payload = logged[0]
    assert event == "invalid_chunk_size"
    assert payload.get("chunk_size") == 0


def test_cli_invalid_limit_logs_error(tmp_path: Path, monkeypatch) -> None:
    """CLI rejects non-positive limits with an error message."""

    input_path = tmp_path / "input.csv"
    input_path.write_text("a\n1\n", encoding="utf-8")
    output_path = tmp_path / "out.csv"

    logged: list[tuple[str, dict[str, object]]] = []

    def fake_error(event: str, *args: object, **kwargs: object) -> None:
        logged.append((event, kwargs))

    monkeypatch.setattr(cli.logger, "error", fake_error)

    exit_code = cli.main(
        [
            "--input",
            str(input_path),
            "--output",
            str(output_path),
            "--chunk-size",
            "5",
            "--limit",
            "0",
            "--log-level",
            "WARNING",
        ]
    )

    assert exit_code == 1
    event, payload = logged[0]
    assert event == "invalid_limit"
    assert payload.get("limit") == 0


def test_cli_missing_output_directory_logs_error(tmp_path: Path, monkeypatch) -> None:
    """CLI reports a missing output directory when creation is disabled."""

    input_path = tmp_path / "input.csv"
    input_path.write_text("a\n1\n", encoding="utf-8")
    output_path = tmp_path / "missing" / "out.csv"
    config_path = tmp_path / "config.yaml"
    config_path.write_text(
        textwrap.dedent(
            """
            local:
              io:
                exist_ok: false
            """
        ),
        encoding="utf-8",
    )

    logged: list[tuple[str, dict[str, object]]] = []

    def fake_error(event: str, *args: object, **kwargs: object) -> None:
        logged.append((event, kwargs))

    monkeypatch.setattr(cli.logger, "error", fake_error)

    exit_code = cli.main(
        [
            "--config",
            str(config_path),
            "--input",
            str(input_path),
            "--output",
            str(output_path),
            "--chunk-size",
            "5",
            "--log-level",
            "WARNING",
        ]
    )

    assert exit_code == 1
    event, payload = logged[0]
    assert event == "output_directory_missing"
    assert payload.get("directory") == str(output_path.parent)


def test_cli_output_parent_not_directory_logs_error(tmp_path: Path, monkeypatch) -> None:
    """CLI reports when the output parent path is not a directory."""

    input_path = tmp_path / "input.csv"
    input_path.write_text("a\n1\n", encoding="utf-8")
    parent_path = tmp_path / "not_a_dir"
    parent_path.write_text("content", encoding="utf-8")
    output_path = parent_path / "out.csv"

    logged: list[tuple[str, dict[str, object]]] = []

    def fake_error(event: str, *args: object, **kwargs: object) -> None:
        logged.append((event, kwargs))

    monkeypatch.setattr(cli.logger, "error", fake_error)

    exit_code = cli.main(
        [
            "--input",
            str(input_path),
            "--output",
            str(output_path),
            "--chunk-size",
            "5",
            "--log-level",
            "WARNING",
        ]
    )

    assert exit_code == 1
    event, payload = logged[0]
    assert event == "output_directory_not_directory"
    assert payload.get("directory") == str(parent_path)
