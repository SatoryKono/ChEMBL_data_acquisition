from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
import get_input_initialisation as cli


def test_run_creates_quality_reports(tmp_path: Path, monkeypatch):
    same_doc = tmp_path / "same.xlsx"
    all_doc = tmp_path / "all.xlsx"
    same_doc.write_text("dummy")
    all_doc.write_text("dummy")
    out_dir = tmp_path / "out"

    tables = {
        "assay": pd.DataFrame({"a": [1, 2, 3], "b": ["x", "y", "z"]}),
        "activity": pd.DataFrame({"x": [1, 2, 3], "y": [1, 2, 3]}),
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
    result = cli.run(args)
    assert result == 0

    for name in tables:
        quality = out_dir / f"{name}_quality_report_table.csv"
        corr = out_dir / f"{name}_data_correlation_report_table.csv"
        assert quality.exists(), quality
        assert corr.exists(), corr
        df = pd.read_csv(quality)
        assert "column" in df.columns
