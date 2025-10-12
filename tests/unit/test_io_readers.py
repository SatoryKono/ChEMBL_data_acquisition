"""Unit tests for :mod:`library.io.readers`."""

from __future__ import annotations

import argparse
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest

from library import io
from scripts import make_activity_postprocessing as activity_cli


@pytest.mark.unit
@pytest.mark.pipeline_scenario("csv_loading")
def test_read_csv__wraps_parser_errors(monkeypatch, tmp_path: Path) -> None:
    sample_path = tmp_path / "input.csv"
    sample_path.write_text("col\nvalue\n", encoding="utf-8")

    def _raise_parser_error(*_args: object, **_kwargs: object) -> pd.DataFrame:
        raise pd.errors.ParserError("boom")

    monkeypatch.setattr(pd, "read_csv", _raise_parser_error)

    cfg = SimpleNamespace(csv_sep=",", csv_encoding="utf-8")

    with pytest.raises(io.CsvReadError) as excinfo:
        io.read_csv(sample_path, cfg=cfg)

    error = excinfo.value
    assert error.path == sample_path
    assert isinstance(error.original_error, pd.errors.ParserError)


@pytest.mark.unit
def test_make_activity_postprocessing_run__returns_non_zero_on_csv_error(
    monkeypatch, tmp_path: Path
) -> None:
    input_path = tmp_path / "activities.csv"
    output_path = tmp_path / "output.csv"

    args = argparse.Namespace(
        input=str(input_path),
        output=str(output_path),
        config=None,
        log_level="INFO",
        run_id="run",
    )

    pipeline_config = SimpleNamespace(pipeline_version="1.0", params={})
    csv_cfg = SimpleNamespace(sep=",", encoding="utf-8", chunksize=1000)

    monkeypatch.setattr(
        activity_cli,
        "get_pipeline_config",
        lambda *_, **__: pipeline_config,
    )
    monkeypatch.setattr(
        activity_cli,
        "get_csv_runtime_config",
        lambda *_, **__: csv_cfg,
    )
    monkeypatch.setattr(
        activity_cli,
        "generate_metrics_report",
        lambda *_, **__: None,
    )

    def _raise_csv_error(*_args, **_kwargs):
        raise io.CsvReadError(input_path, ValueError("csv failed"))

    monkeypatch.setattr(
        activity_cli,
        "run_postprocessing_pipeline",
        _raise_csv_error,
    )

    exit_code = activity_cli.run(args)

    assert exit_code == 1
