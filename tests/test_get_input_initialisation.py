from __future__ import annotations

import argparse
import io
import json
from pathlib import Path

import pandas as pd

from library.cli import LoggerConfig, configure_logger
from library.config import Config
from scripts import get_input_initialisation as cli


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
    cfg = Config(api={"user_agent": "test/1.0 (mailto:test@example.org)"})
    result = cli.run(cfg, args)
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
    cfg = Config(api={"user_agent": "test/1.0 (mailto:test@example.org)"})
    result = cli.run(cfg, args)
    assert result == 1
    lines = buf.getvalue().splitlines()
    assert lines
    record = json.loads(lines[-1])
    assert "required table 'activity' missing" in record.get("msg", "")


def test_run_missing_columns_logs_error(tmp_path: Path, monkeypatch) -> None:
    """``run`` should report missing required columns."""
    same_doc = tmp_path / "same.xlsx"
    all_doc = tmp_path / "all.xlsx"
    same_doc.write_text("dummy")
    all_doc.write_text("dummy")

    def fake_load_same_doc(path: Path):  # pragma: no cover - simple stub
        return {}

    def fake_load_all_doc(path: Path):  # pragma: no cover - simple stub
        return {}

    def fake_build_combined_tables(*_args, **_kwargs):
        raise KeyError("['target_chembl_id', 'gene_index'] not in index")

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
    cfg = Config(api={"user_agent": "test/1.0 (mailto:test@example.org)"})
    result = cli.run(cfg, args)
    assert result == 1
    lines = buf.getvalue().splitlines()
    assert lines
    record = json.loads(lines[-1])
    msg = record.get("msg", "")
    assert "required column(s) missing" in msg
    assert "target_chembl_id" in msg


def test_run_uses_config_output_dir(tmp_path: Path, monkeypatch) -> None:
    same_doc = tmp_path / "same.xlsx"
    all_doc = tmp_path / "all.xlsx"
    same_doc.write_text("dummy")
    all_doc.write_text("dummy")

    cfg = Config(api={"user_agent": "test/1.0 (mailto:test@example.org)"})
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


def test_main_missing_dictionary_dir(tmp_path: Path, monkeypatch) -> None:
    """``main`` should handle absent dictionary directory gracefully."""

    cfg_path = tmp_path / "cfg.yaml"
    cfg_path.write_text('log:\n  level: "ERROR"\n')
    same_doc = tmp_path / "same.xlsx"
    all_doc = tmp_path / "all.xlsx"
    same_doc.write_text("dummy")
    all_doc.write_text("dummy")
    out_dir = tmp_path / "out"

    # Skip heavy processing by stubbing ``run``
    monkeypatch.setattr(cli, "run", lambda _cfg, _args: 0)

    result = cli.main(
        [
            "--config",
            str(cfg_path),
            "--same-doc",
            str(same_doc),
            "--all-doc",
            str(all_doc),
            "--out-dir",
            str(out_dir),
        ]
    )
    assert result == 0
