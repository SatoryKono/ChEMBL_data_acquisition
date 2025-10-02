"""Tests for :mod:`library.cli.parser.prepare_io_paths`."""

from __future__ import annotations

import argparse
import warnings
from io import StringIO

import pytest

from library.cli import prepare_io_paths
from library.cli import parser as cli_parser
from library.logging_setup import Logger, LoggerConfig


def _build_namespace() -> argparse.Namespace:
    return argparse.Namespace(
        base_path=None,
        input_dir=None,
        output_dir=None,
        _deprecated_out=argparse.SUPPRESS,
        output_csv="output.csv",
        final_out=None,
        input_csv="input.csv",
        raw_format=None,
        raw_out=None,
        date=None,
        force=False,
        skip_existing=False,
    )


def test_prepare_io_paths_warns_when_logger_stream_closed(monkeypatch: pytest.MonkeyPatch) -> None:
    closed_stream = StringIO()
    closed_stream.close()
    closed_logger = Logger(LoggerConfig(stream=closed_stream))
    monkeypatch.setattr(cli_parser.log, "logger", closed_logger)

    args = _build_namespace()

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", DeprecationWarning)
        prepare_io_paths(args, output_stem="assays")

    assert any(
        "--final-out" in str(item.message)
        for item in caught
    ), "Deprecation warning should fall back to warnings when logger stream is closed"
