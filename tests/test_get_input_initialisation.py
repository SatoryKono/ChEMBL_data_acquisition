from __future__ import annotations

import argparse
import json
import io
from pathlib import Path

import pandas as pd
from library.cli import LoggerConfig, configure_logger
import get_input_initialisation as cli
from library.config import Config


def test_run_creates_quality_reports(tmp_path: Path, monkeypatch) -> None:
    """``run`` should produce quality reports for all tables."""
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

    def fake_load_same_doc(
        path: Path,
    ) -> dict[str, object]:  # pragma: no cover - simple stub
        return {}

    def fake_load_all_doc(
        path: Path,
    ) -> dict[str, object]:  # pragma: no cover - simple stub
        return {}

    def fake_build_combined_tables(*_args, **_kwargs) -> dict[str, pd.DataFrame]:
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
    result = cli.run(Config(), args)
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


def test_run_missing_activity_logs_error(tmp_path: Path, monkeypatch) -> None:
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
    buf = io.StringIO()
    configure_logger(LoggerConfig(stream=buf))
    result = cli.run(Config(), args)
    assert result == 1
    lines = buf.getvalue().splitlines()
    assert lines
    record = json.loads(lines[-1])
    assert "required table 'activity' missing" in record.get("msg", "")


def test_run_uses_config_output_dir(tmp_path: Path, monkeypatch) -> None:
    same_doc = tmp_path / "same.xlsx"
    all_doc = tmp_path / "all.xlsx"
    same_doc.write_text("dummy")
    all_doc.write_text("dummy")

    cfg = Config()
    cfg.init.output_dir = tmp_path / "default"

    tables = {"assay": pd.DataFrame({"id": [1]})}

    monkeypatch.setattr(cli.lib, "load_same_doc", lambda _p: {})
    monkeypatch.setattr(cli.lib, "load_all_doc", lambda _p: {})
    monkeypatch.setattr(cli.lib, "build_combined_tables", lambda *_a, **_k: tables)

    args = argparse.Namespace(
        same_doc=same_doc,
        all_doc=all_doc,
        out_dir=None,
        format="csv",
        dictionary_dir=tmp_path,
    )

    assert cli.run(cfg, args) == 0
    assert (cfg.init.output_dir / "assay.csv").exists()
