"""Integration tests for the activity CLI entry point."""

from __future__ import annotations

import contextlib
import io
from pathlib import Path

import pytest

import library.cli.base as cli_base
import library.cli.entrypoints.activity as activity_entrypoint
from library.cli import parser as parser_module
from library.cli.entrypoints.activity import ActivityPipelineCLI
from library.cli.logging import CLILoggingContext
import library.cli.utils as cli_utils
from library.config.loader import DEFAULT_CONFIG_PATH


@pytest.mark.integration
def test_activity_cli__require_stamp_mode_without_date(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Config-driven stamp mode should append the date when CLI omits it."""

    input_csv = tmp_path / "activity.csv"
    input_csv.write_text("activity_chembl_id\nCHEMBL1\n", encoding="utf-8")

    config_path = DEFAULT_CONFIG_PATH

    original_load_config = parser_module.load_config

    def _fake_load_config(path, **kwargs):
        cfg, metadata = original_load_config(path, **kwargs)
        io_cfg = getattr(cfg, "io", None)
        if io_cfg is not None:
            setattr(io_cfg, "output_stamp_mode", "require")
            setattr(io_cfg, "default_date_prefix", "20240101")
        local_cfg = getattr(cfg, "local", None)
        if local_cfg is not None:
            local_io = getattr(local_cfg, "io", None)
            if local_io is not None:
                setattr(local_io, "output_stamp_mode", "require")
                setattr(local_io, "default_date_prefix", "20240101")
        metadata.snapshot.setdefault("local", {}).setdefault("io", {})["output_stamp_mode"] = "require"
        metadata.snapshot["local"]["io"]["default_date_prefix"] = "20240101"
        return cfg, metadata

    monkeypatch.setattr(parser_module, "load_config", _fake_load_config)
    monkeypatch.setattr(cli_utils, "ensure_dirs", lambda _cfg: None)

    captured: dict[str, object] = {}

    def _fake_run(cfg, args):
        final_out = getattr(args, "final_out", None)
        output_csv = getattr(args, "output_csv", None)
        captured["final_out_name"] = (
            Path(final_out).name if isinstance(final_out, Path) else str(final_out)
        )
        captured["output_csv_name"] = (
            Path(output_csv).name if isinstance(output_csv, Path) else str(output_csv)
        )
        captured["date"] = getattr(args, "date", None)
        captured["stamp_mode"] = getattr(args, "output_stamp_mode", None)
        return 0

    monkeypatch.setattr(activity_entrypoint, "run", _fake_run)

    @contextlib.contextmanager
    def _fake_setup_cli_logging(_program: str, log_cfg, _date_token):
        yield CLILoggingContext(
            log_path=tmp_path / "activity.log",
            log_cfg=log_cfg,
            console_stream=io.StringIO(),
        )

    monkeypatch.setattr(cli_base, "setup_cli_logging", _fake_setup_cli_logging)

    cli = ActivityPipelineCLI()
    exit_code = cli.main(
        [
            "--config",
            str(config_path),
            "--input",
            str(input_csv),
            "--output-dir",
            str(tmp_path),
        ]
    )

    assert exit_code == 0
    assert captured["stamp_mode"] == "require"
    assert captured["date"] == "20240101"
    assert captured["final_out_name"] == "output.activity_20240101.csv"
    assert captured["output_csv_name"] == "output.activity_20240101.csv"
