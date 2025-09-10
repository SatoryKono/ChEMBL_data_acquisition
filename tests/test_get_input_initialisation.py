from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
import get_input_initialisation as cli
from library.config import Config


def test_run_creates_quality_reports(tmp_path: Path, monkeypatch):
    same_doc = tmp_path / "same.xlsx"
    all_doc = tmp_path / "all.xlsx"
    same_doc.write_text("dummy")
    all_doc.write_text("dummy")
    out_dir = tmp_path / "out"

    tables = {
        "assay": pd.DataFrame({"a": [1, 2, 3], "b": ["x", "y", "z"]}),
        "activity": pd.DataFrame({"x": [1, 2, 3], "y": [1, 2, 3]}),
        "pairs_same_document": pd.DataFrame({"id": [3, 4]}),
        "pairs_independent": pd.DataFrame({"id": [5]}),
        "pairs_non_independent": pd.DataFrame({"id": [6]}),
        "activity_independent_status": pd.DataFrame({"id": [7]}),
    }

    def fake_load_same_doc(path: Path):  # pragma: no cover - simple stub
        return {}

    def fake_load_all_doc(path: Path):  # pragma: no cover - simple stub
        return {}

    def fake_build_combined_tables(*_args, **_kwargs):
        return tables

    monkeypatch.setattr(cli.lib, "load_same_doc", fake_load_same_doc)
    monkeypatch.setattr(cli.lib, "load_all_doc", fake_load_all_doc)
    monkeypatch.setattr(cli.lib, "build_combined_tables", fake_build_combined_tables)

    args = argparse.Namespace(
        same_doc=same_doc,
        all_doc=all_doc,
        out_dir=out_dir,
        format="csv",
        dictionary_dir=tmp_path,
    )
    result = cli.run(args, Config())
    assert result == 0

    assert (
        out_dir / "status" / "independent" / "activity_independent_status.csv"
    ).exists()

    for name in tables:
        quality = out_dir / "data_validity_report" / f"{name}_quality_report_table.csv"
        corr = (
            out_dir
            / "data_validity_report"
            / f"{name}_data_correlation_report_table.csv"
        )
        assert quality.exists(), quality
        assert corr.exists(), corr
        df = pd.read_csv(quality)
        assert "column" in df.columns


def test_run_missing_activity_logs_error(tmp_path: Path, monkeypatch, caplog) -> None:
    """``run`` should report missing required tables."""
    same_doc = tmp_path / "same.xlsx"
    all_doc = tmp_path / "all.xlsx"
    same_doc.write_text("dummy")
    all_doc.write_text("dummy")

    def fake_load_same_doc(path: Path):  # pragma: no cover - simple stub
        return {}

    def fake_load_all_doc(path: Path):  # pragma: no cover - simple stub
        return {}

    def fake_build_combined_tables(*_args, **_kwargs):
        raise KeyError("activity")

    monkeypatch.setattr(cli.lib, "load_same_doc", fake_load_same_doc)
    monkeypatch.setattr(cli.lib, "load_all_doc", fake_load_all_doc)
    monkeypatch.setattr(cli.lib, "build_combined_tables", fake_build_combined_tables)

    args = argparse.Namespace(
        same_doc=same_doc,
        all_doc=all_doc,
        out_dir=tmp_path / "out",
        format="csv",
        dictionary_dir=tmp_path,
    )
    with caplog.at_level("ERROR"):
        result = cli.run(args, Config())
    assert result == 1
    assert "required table 'activity' missing" in caplog.text
