"""Unit tests for :mod:`library.cli_utils`."""

from __future__ import annotations

import argparse
from collections.abc import Iterable, Sequence
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import Mock

import pandas as pd
import pandas.testing as tm
import pytest

from library import cli_utils
from library.cli.pipeline_definition import PipelineDefinition
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


@pytest.mark.unit
@pytest.mark.pipeline_scenario(
    "csv_loading",
    "normalization",
    "enrichment",
    "transformation_rules",
    "missing_data",
    "logging",
    "assembly",
    "export",
    "degradation",
    "idempotence",
)
def test_parquet_chunk_store__fallback_to_pickle_when_engine_missing(
    monkeypatch: pytest.MonkeyPatch, caplog: pytest.LogCaptureFixture
) -> None:
    """The chunk store must remain usable when parquet engines are absent."""

    def _raise_import_error(*_args, **_kwargs):
        raise ImportError("no parquet engine available")

    monkeypatch.setattr("pandas.io.parquet.get_engine", _raise_import_error)

    caplog.set_level("WARNING", logger="chembl")

    store = cli_utils._ParquetChunkStore()
    frame = pd.DataFrame({"value": [1, 2], "label": ["a", "b"]})

    store.append(frame)

    assert store.row_count == len(frame)
    assert store._backend == "pickle"
    assert store._paths[0].suffix == ".pkl"

    (restored,) = list(store.iter_frames())
    tm.assert_frame_equal(restored, frame.reset_index(drop=True))

    fallback_messages = [rec.message for rec in caplog.records]
    assert any(
        "chunk_store_backend_fallback" in message for message in fallback_messages
    )


@pytest.mark.unit
@pytest.mark.parametrize("strict_mode", [False, True])
def test_run_pipeline__metadata_hook_must_return_dataframe(
    strict_mode: bool, tmp_path: Path
) -> None:
    output_path = tmp_path / "output.csv"
    failure_path = tmp_path / "failures.csv"

    frame = pd.DataFrame({"value": [1]})

    writer_calls: list[list[pd.DataFrame]] = []

    def writer(
        chunks: Iterable[pd.DataFrame],
        destination: Path,
        column_order: Sequence[str] | None,
        resolved_keys: Sequence[str],
    ) -> Path:
        writer_calls.append(list(chunks))
        destination.write_text("should not be written", encoding="utf-8")
        return destination

    def broken_hook(df: pd.DataFrame) -> pd.DataFrame | None:
        return None

    definition = PipelineDefinition(
        schema=None,
        schema_name="test",
        writer=writer,
        metadata_hooks=[broken_hook],
        table_quality=lambda _path: None,
        strict_mode=strict_mode,
        command="unit-test",
    )

    logger = Mock()

    result = cli_utils.run_pipeline(
        definition=definition,
        fetcher=lambda: [frame],
        output_path=output_path,
        failure_path=failure_path,
        emit_standard_outputs=False,
        emit_legacy_artifacts=False,
        logger=logger,
    )

    assert result.exit_code == 1
    assert not writer_calls, "writer must not be invoked when hook output is invalid"

    invalid_calls = [
        call
        for call in logger.error.call_args_list
        if call.args and call.args[0] == "metadata_hook_invalid_return"
    ]
    assert invalid_calls, "metadata_hook_invalid_return should be logged"

    invalid_kwargs = invalid_calls[0].kwargs
    assert invalid_kwargs["hook"].endswith("broken_hook")
    assert "broken_hook" in invalid_kwargs["error"]
    assert invalid_kwargs["strict_mode"] is strict_mode
