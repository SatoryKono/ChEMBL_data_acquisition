"""CLI smoke tests for chunk_io_main."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

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


def test_cli_uses_config_separator(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """CLI honours ``io.csv_sep`` from the provided configuration file."""

    config_path = tmp_path / "config.yaml"
    output_dir = tmp_path / "output"
    cache_dir = tmp_path / "cache"
    config_path.write_text(
        "local:\n"
        "  io:\n"
        "    csv_sep: ';'\n"
        f"    output_dir: '{output_dir.as_posix()}'\n"
        f"    cache_dir: '{cache_dir.as_posix()}'\n"
        "system:\n"
        "  log:\n"
        "    level: INFO\n"
    )

    input_path = tmp_path / "input.csv"
    df = pd.DataFrame({"a": ["x", "y"], "b": [1, 2]})
    df.to_csv(input_path, index=False, sep=";")

    captured: dict[str, object] = {}

    def _stub_process(
        input_csv: Path,
        output_csv: Path,
        *,
        cfg,
        chunk_size: int,
        limit,
        checkpoint_path,
        sep,
        encoding,
    ) -> int:
        captured.update({"sep": sep, "cfg_sep": cfg.csv_sep})
        return 0

    monkeypatch.setattr(cli, "process_csv_chunks", _stub_process)

    exit_code = cli.main(
        [
            "--config",
            str(config_path),
            "--input",
            str(input_path),
            "--output",
            str(tmp_path / "result.csv"),
            "--checkpoint",
            str(tmp_path / "checkpoint.json"),
        ]
    )

    assert exit_code == 0
    assert captured["sep"] == ";"
    assert captured["cfg_sep"] == ";"
