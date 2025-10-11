"""Unit tests for :mod:`library.cli_utils`."""

from __future__ import annotations

import argparse
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import Mock

import pytest

from library import cli_utils
from library.cli import LoggerConfig


@pytest.mark.unit
def test_run_cli_command__logs_exc_info_on_value_error(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Ensure ``exc_info`` is attached when configuration resolution fails."""

    parser = argparse.ArgumentParser()
    args = argparse.Namespace(config=object(), verbose=False, log_level=None)
    log_cfg = LoggerConfig(level="INFO", run_id="test-run", redact_secrets=False)

    logger = Mock()
    monkeypatch.setattr(cli_utils.cli, "configure_logger", lambda cfg: logger)

    exit_code = cli_utils.run_cli_command(
        args=args,
        parser=parser,
        log_cfg=log_cfg,
        mapping={},
        run=lambda *_: 0,
        logger=logger,
    )

    assert exit_code == 1
    error_call = logger.error.call_args
    assert error_call is not None
    exc_info = error_call.kwargs.get("exc_info")
    assert exc_info is not None
    assert isinstance(exc_info, ValueError)
    assert error_call.kwargs["error"] == "configuration path must be provided"


@pytest.mark.unit
def test_run_cli_command__run_id_determinism(tmp_path, monkeypatch):
    parser = argparse.ArgumentParser(prog="activity")
    config_path = tmp_path / "config.yaml"
    config_path.write_text("{}", encoding="utf-8")

    captured_cfgs: list[LoggerConfig] = []

    class DummyLogger:
        def info(self, *args, **kwargs):
            return None

        def error(self, *args, **kwargs):
            return None

        def warning(self, *args, **kwargs):
            return None

        def debug(self, *args, **kwargs):
            return None

    dummy_logger = DummyLogger()

    def fake_configure_logger(cfg: LoggerConfig) -> DummyLogger:
        captured_cfgs.append(cfg)
        return dummy_logger

    monkeypatch.setattr(cli_utils.cli, "configure_logger", fake_configure_logger)
    monkeypatch.setattr(
        cli_utils.cli,
        "apply_config_overrides",
        lambda *args, **kwargs: SimpleNamespace(io=SimpleNamespace()),
    )
    monkeypatch.setattr(cli_utils, "ensure_dirs", lambda cfg: None)

    def fake_run(cfg, args):
        return 0

    def make_namespace(
        *, input_name: str, output_name: str, limit: int
    ) -> argparse.Namespace:
        input_path = tmp_path / input_name
        input_path.write_text("activity_id\nACT1\n", encoding="utf-8")
        output_path = tmp_path / output_name
        namespace = argparse.Namespace(
            config=str(config_path),
            log_level="INFO",
            verbose=False,
            run_id=None,
            input_csv=input_path,
            final_out=output_path,
            output_csv=output_path,
            base_path=None,
            input_dir=None,
            output_dir=None,
            cache_dir=None,
            raw_out=None,
            date=None,
            sep=",",
            encoding="utf8",
            raw_format="csv",
            force=False,
            skip_existing=False,
            print_config=False,
            limit=limit,
            invocation=(
                parser.prog,
                "--input",
                str(input_path),
                "--final-out",
                str(output_path),
                "--limit",
                str(limit),
            ),
        )
        return namespace

    def run_and_get_id(namespace: argparse.Namespace) -> str:
        captured_cfgs.clear()
        log_cfg = LoggerConfig(level="INFO")
        exit_code = cli_utils.run_cli_command(
            args=namespace,
            parser=parser,
            base_parser=None,
            log_cfg=log_cfg,
            mapping={},
            run=fake_run,
            logger=dummy_logger,
        )
        assert exit_code == 0
        assert captured_cfgs, "configure_logger must be invoked"
        return captured_cfgs[-1].run_id

    first_args = make_namespace(input_name="input.csv", output_name="out.csv", limit=5)
    second_args = make_namespace(input_name="input.csv", output_name="out.csv", limit=5)
    third_args = make_namespace(input_name="input.csv", output_name="out.csv", limit=7)

    run_id_a = run_and_get_id(first_args)
    run_id_b = run_and_get_id(second_args)
    run_id_c = run_and_get_id(third_args)

    assert run_id_a == run_id_b
    assert run_id_a != run_id_c


@pytest.mark.unit
def test_resolve_standard_output_naming__temporary_working_path() -> None:
    output_path = Path(".output.targets_20240101.csv.tmp")

    table_name, date_tag = cli_utils._resolve_standard_output_naming(
        output_path,
        cfg=None,
        command="target",
        invocation=("target", "--final-out", str(output_path)),
    )

    assert table_name == "targets"
    assert date_tag == "20240101"
