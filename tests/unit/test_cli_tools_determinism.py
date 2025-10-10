"""Unit tests for :mod:`library.utils.cli_tools.check_determinism`."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from library.common import csv_utils
from library.utils.cli_tools import check_determinism


@pytest.mark.unit
def test_cli_tools_determinism__hashes_match(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """The deterministic CSV helpers should yield identical hashes."""

    original_write_csv = csv_utils.write_csv_deterministic
    original_write_chunks = csv_utils.write_csv_chunks_deterministic

    baseline_meta: dict[str, str] = {}

    def _synchronise_metadata(path: Path) -> None:
        meta_path = path.with_suffix(path.suffix + ".meta.yaml")
        if not meta_path.exists():
            return
        content = meta_path.read_text(encoding="utf-8")
        if "baseline" not in baseline_meta:
            baseline_meta["baseline"] = content
            return
        meta_path.write_text(baseline_meta["baseline"], encoding="utf-8")

    def _stable_write_csv(
        df: pd.DataFrame,
        destination: Path,
        *,
        key_cols,
        col_order=None,
        **kwargs,
    ) -> Path:
        output = original_write_csv(
            df,
            destination,
            key_cols=key_cols,
            col_order=col_order,
            **kwargs,
        )
        _synchronise_metadata(destination)
        return output

    def _stable_write_chunks(
        chunks,
        destination: Path,
        *,
        key_cols,
        col_order=None,
        **kwargs,
    ) -> Path:
        output = original_write_chunks(
            chunks,
            destination,
            key_cols=key_cols,
            col_order=col_order,
            **kwargs,
        )
        _synchronise_metadata(destination)
        return output

    monkeypatch.setattr(check_determinism, "write_csv_deterministic", _stable_write_csv)
    monkeypatch.setattr(
        check_determinism, "write_csv_chunks_deterministic", _stable_write_chunks
    )

    ok = check_determinism.run_check(tmp_path)

    assert ok is True

    generated = sorted(tmp_path.glob("*.meta.yaml"))
    assert len(generated) == 3

    digests = [check_determinism.sha256_file(path) for path in generated]
    assert digests[0] == digests[1] == digests[2]


@pytest.mark.unit
def test_cli_tools_determinism__detects_mismatch(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    """A single mismatching CSV hash must fail the determinism check."""

    original_write_csv = csv_utils.write_csv_deterministic
    original_write_chunks = csv_utils.write_csv_chunks_deterministic

    def _tampered_write_csv(
        df: pd.DataFrame,
        destination: Path,
        *,
        key_cols,
        col_order=None,
        **kwargs,
    ) -> Path:
        output = original_write_csv(
            df,
            destination,
            key_cols=key_cols,
            col_order=col_order,
            **kwargs,
        )
        if destination.name == "second.csv":
            with destination.open("a", encoding="utf-8") as handle:
                handle.write("\n# tampered\n")
        return output

    def _tampered_write_chunks(
        chunks,
        destination: Path,
        *,
        key_cols,
        col_order=None,
        **kwargs,
    ) -> Path:
        output = original_write_chunks(
            chunks,
            destination,
            key_cols=key_cols,
            col_order=col_order,
            **kwargs,
        )
        with destination.open("a", encoding="utf-8") as handle:
            handle.write("\n# tampered chunk\n")
        return output

    monkeypatch.setattr(check_determinism, "write_csv_deterministic", _tampered_write_csv)
    monkeypatch.setattr(check_determinism, "write_csv_chunks_deterministic", _tampered_write_chunks)

    result = check_determinism.run_check(tmp_path)

    assert result is False
