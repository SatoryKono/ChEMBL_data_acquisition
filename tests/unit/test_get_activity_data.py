"""Unit tests for :mod:`scripts.get_activity_data`."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterable

import pytest

from scripts import get_activity_data


def _make_args(tmp_path: Path) -> argparse.Namespace:
    """Return an :class:`argparse.Namespace` populated with sane defaults."""

    input_csv = tmp_path / "input.csv"
    input_csv.write_text("activity_id\nA1\n", encoding="utf-8")
    return argparse.Namespace(
        input_csv=input_csv,
        output_csv=tmp_path / "output.csv",
        offset=0,
        workers=None,
        skip_existing=False,
        force=False,
        dry_run=False,
        invocation=None,
    )


@pytest.mark.parametrize(
    "invocation, expected",
    [
        (None, (get_activity_data.PROGRAM_NAME,)),
        ((), ()),
        ([], ()),
        (["run"], ("run",)),
        (("fetch", 2), ("fetch", "2")),
        ([Path("/tmp/file.csv")], ("/tmp/file.csv",)),
        ([""], ("",)),
        (["--limit", 5], ("--limit", "5")),
        (["Σ", Path("relative")], ("Σ", "relative")),
        ([42.0, 0], ("42.0", "0")),
    ],
)
def test_args_invocation__cases(invocation: Iterable[object] | None, expected: tuple[str, ...]) -> None:
    namespace = argparse.Namespace(invocation=invocation)
    assert get_activity_data._args_invocation(namespace) == expected


def test_run_chembl__invalid_limit_logs_error(cfg, tmp_path, caplog) -> None:
    args = _make_args(tmp_path)
    cfg.activity.limit = -5

    caplog.set_level("ERROR")
    exit_code = get_activity_data.run_chembl(cfg, args)

    assert exit_code == 1


def test_run_chembl__dry_run_short_circuits(cfg, tmp_path, caplog) -> None:
    args = _make_args(tmp_path)
    cfg.activity.dry_run = True
    cfg.activity.limit = 7

    caplog.set_level("INFO")
    exit_code = get_activity_data.run_chembl(cfg, args)

    assert exit_code == 0


def test_run__skip_existing_respects_flags(cfg, tmp_path, monkeypatch) -> None:
    args = _make_args(tmp_path)
    args.output_csv = None
    args.skip_existing = True
    args.force = False

    calls: list[str] = []

    def fake_default_output_path(input_path: Path, _io_cfg) -> Path:  # pragma: no cover - exercised via tests
        assert input_path == args.input_csv
        destination = tmp_path / "output.csv"
        destination.write_text("existing", encoding="utf-8")
        return destination

    monkeypatch.setattr(get_activity_data.io, "default_output_path", fake_default_output_path)
    monkeypatch.setattr(get_activity_data, "run_chembl", lambda *_: calls.append("run") or 0)

    exit_code = get_activity_data.run(cfg, args)

    assert exit_code == 0
    assert calls == []
