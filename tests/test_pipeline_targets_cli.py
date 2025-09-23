"""CLI tests for :mod:`scripts.pipeline_targets_main`."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd
import pytest

from library.config import Config
from library.pipeline_targets import PipelineResult
from scripts import pipeline_targets_main as cli


class _DummyLogger:
    def info(self, *args: Any, **kwargs: Any) -> None:  # pragma: no cover - trivial
        pass


def test_cli_forwards_batch_size(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    input_csv = tmp_path / "targets.csv"
    input_csv.write_text("target_chembl_id\nCHEMBL1\n", encoding="utf8")
    output_csv = tmp_path / "out.csv"

    captured: dict[str, Any] = {}

    def fake_run_pipeline(
        chunk_iterator: Any,
        cfg: Config,
        *,
        chembl_fetcher: Any,
        batch_size: int | None,
        **_: Any,
    ) -> PipelineResult:
        captured["chunks"] = [list(chunk) for chunk in chunk_iterator()]
        captured["batch_size"] = batch_size
        return PipelineResult(chembl=pd.DataFrame({"target_chembl_id": ["CHEMBL1"]}))

    def fake_write_csv(
        df: pd.DataFrame,
        path: Path | str,
        *,
        cfg: Config,
        sep: str | None = None,
        encoding: str | None = None,
    ) -> Path:
        captured["written_path"] = Path(path)
        captured["written_df"] = df.copy()
        return Path(path)

    monkeypatch.setattr(cli, "configure_logger", lambda cfg: _DummyLogger())
    monkeypatch.setattr(cli, "apply_config_overrides", lambda *a: Config())
    monkeypatch.setattr(cli, "ensure_dirs", lambda cfg: None)
    monkeypatch.setattr(cli, "print_config", lambda cfg: None)
    monkeypatch.setattr(cli, "run_pipeline", fake_run_pipeline)
    monkeypatch.setattr(cli, "write_csv", fake_write_csv)

    args = [
        "--input",
        str(input_csv),
        "--output",
        str(output_csv),
        "--batch-size",
        "25",
    ]
    exit_code = cli.main(args)
    assert exit_code == 0
    assert captured["batch_size"] == 25
    assert captured["chunks"] == [["CHEMBL1"]]
    assert captured["written_path"] == output_csv
    assert list(captured["written_df"]["target_chembl_id"]) == ["CHEMBL1"]
