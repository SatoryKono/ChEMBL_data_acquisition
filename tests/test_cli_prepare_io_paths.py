"""Tests for :mod:`library.cli.parser.prepare_io_paths`."""

from __future__ import annotations

import argparse
import warnings
from io import StringIO
from pathlib import Path

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


def _namespace_with(**overrides: object) -> argparse.Namespace:
    namespace = _build_namespace()
    for key, value in overrides.items():
        setattr(namespace, key, value)
    return namespace


def test_prepare_io_paths_prefers_final_out_over_output(tmp_path: Path) -> None:
    args = _namespace_with(
        base_path=tmp_path,
        input_csv=tmp_path / "input.csv",
        output_csv="legacy.csv",
        final_out="custom.csv",
    )

    prepare_io_paths(args, output_stem="assays")

    expected = (tmp_path / "custom.csv").resolve()
    assert args.final_out == expected
    assert args.output_csv == expected


def test_prepare_io_paths_uses_deprecated_output_when_final_missing(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    warning_calls: list[None] = []

    def fake_warning() -> None:
        warning_calls.append(None)

    monkeypatch.setattr(cli_parser, "_emit_deprecated_output_warning", fake_warning)

    args = _namespace_with(
        base_path=tmp_path,
        input_csv=tmp_path / "input.csv",
        output_csv="legacy.csv",
        final_out=None,
    )

    prepare_io_paths(args, output_stem="assays")

    expected = (tmp_path / "legacy.csv").resolve()
    assert args.final_out == expected
    assert args.output_csv == expected
    assert warning_calls == [None]


def test_prepare_io_paths_handles_final_out_only(tmp_path: Path) -> None:
    args = _namespace_with(
        base_path=tmp_path,
        input_csv=tmp_path / "input.csv",
        output_csv=None,
        final_out="final.csv",
    )

    prepare_io_paths(args, output_stem="assays")

    expected = (tmp_path / "final.csv").resolve()
    assert args.final_out == expected
    assert args.output_csv == expected
