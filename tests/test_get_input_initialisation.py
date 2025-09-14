from __future__ import annotations

import argparse
import io
import json
from pathlib import Path

import pandas as pd
import pytest

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
        "activity_independent_status": pd.DataFrame(
            {"Filtered.new": ["good", "bad"], "independent_IC50": [1, 0]}
        ),
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
    cfg = Config.model_validate({"api": {"user_agent": "test@example.org"}})
    result = cli.run(cfg, args)
    assert result == 0

    assert (out_dir / "independent" / "activity_independent.csv").exists()

    assert (
        out_dir
        / "status"
        / "independent"
        / "activity_independent_status_statistics.csv"
    ).exists()

    expected = set(tables)
    expected.remove("activity_independent_status")
    expected.update({"activity_independent", "activity_independent_status_statistics"})

    for name in expected:
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
    cfg = Config.model_validate({"api": {"user_agent": "test@example.org"}})
    result = cli.run(cfg, args)
    assert result == 1
    lines = buf.getvalue().splitlines()
    assert lines
    record = json.loads(lines[-1])
    assert record["event"] == "missing_table"
    assert record["table"] == "activity"


@pytest.mark.parametrize(
    "error_msg",
    [
        "activity table missing expected columns: foo, bar",
        "required column(s) baz missing",
    ],
)
def test_run_missing_columns_logs_specific_error(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, error_msg: str
) -> None:
    """``run`` should report missing required columns with a specific message."""
    same_doc = tmp_path / "same.xlsx"
    all_doc = tmp_path / "all.xlsx"
    same_doc.write_text("dummy")
    all_doc.write_text("dummy")

    def fake_load_same_doc(path: Path) -> dict[str, object]:  # pragma: no cover
        return {}

    def fake_load_all_doc(path: Path) -> dict[str, object]:  # pragma: no cover
        return {}

    def fake_build_combined_tables(*_args, **_kwargs):
        raise KeyError(error_msg)

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
    cfg = Config.model_validate({"api": {"user_agent": "test@example.org"}})
    result = cli.run(cfg, args)
    assert result == 1
    lines = buf.getvalue().splitlines()
    assert lines
    record = json.loads(lines[-1])
    assert record["event"] == "missing_columns"
    assert record["details"] == error_msg


def test_run_uses_config_output_dir(tmp_path: Path, monkeypatch) -> None:
    same_doc = tmp_path / "same.xlsx"
    all_doc = tmp_path / "all.xlsx"
    same_doc.write_text("dummy")
    all_doc.write_text("dummy")

    cfg = Config.model_validate({"api": {"user_agent": "test@example.org"}})
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
    monkeypatch.setenv("CHEMBL_DA__API__USER_AGENT", "test@example.org")

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


def test_status_merge_preserves_columns(tmp_path: Path, monkeypatch) -> None:
    """Status merging should retain original entity columns."""

    same_doc = tmp_path / "same.xlsx"
    all_doc = tmp_path / "all.xlsx"
    same_doc.write_text("dummy")
    all_doc.write_text("dummy")
    out_dir = tmp_path / "out"

    id_cols = {
        "activity": "activity_chembl_id",
        "assay": "assay_chembl_id",
        "document": "document_chembl_id",
        "target": "target_chembl_id",
        "testitem": "molecule_chembl_id",
    }
    tables: dict[str, pd.DataFrame] = {}
    for entity, id_col in id_cols.items():
        base = pd.DataFrame({id_col: [f"{entity}1"], "orig": [1]})
        status = pd.DataFrame(
            {id_col: [f"{entity}1"], "Filtered.new": ["good"], "independent_IC50": [1]}
        )
        tables[f"{entity}_independent"] = base
        tables[f"{entity}_independent_status"] = status

    monkeypatch.setattr(cli.lib, "load_same_doc", lambda _p: {})
    monkeypatch.setattr(cli.lib, "load_all_doc", lambda _p: {})
    monkeypatch.setattr(cli.lib, "build_combined_tables", lambda *_a, **_k: tables)
    monkeypatch.setattr(cli.lib, "generate_pair_entity_tables", lambda t, _m: t)
    monkeypatch.setattr(cli, "analyze_table_quality", lambda *_a, **_k: None)

    args = argparse.Namespace(
        same_doc=same_doc,
        all_doc=all_doc,
        out_dir=out_dir,
        format="csv",
        dictionary_dir=tmp_path,
    )
    cfg = Config.model_validate({"api": {"user_agent": "test@example.org"}})
    assert cli.run(cfg, args) == 0

    for entity, id_col in id_cols.items():
        df = pd.read_csv(out_dir / "independent" / f"{entity}_independent.csv")
        assert id_col in df.columns
        assert "orig" in df.columns


def test_run_handles_system_status_table(tmp_path: Path, monkeypatch) -> None:
    """System status tables should be merged or emitted with logging."""

    same_doc = tmp_path / "same.xlsx"
    all_doc = tmp_path / "all.xlsx"
    same_doc.write_text("dummy")
    all_doc.write_text("dummy")
    out_dir = tmp_path / "out"

    tables = {
        "system_independent_status": pd.DataFrame(
            {
                "system_id": ["T1_TG1_IC50"],
                "Filtered.new": ["good"],
                "independent_IC50": [1],
            }
        )
    }

    monkeypatch.setattr(cli.lib, "load_same_doc", lambda _p: {})
    monkeypatch.setattr(cli.lib, "load_all_doc", lambda _p: {})
    monkeypatch.setattr(cli.lib, "build_combined_tables", lambda *_a, **_k: tables)
    monkeypatch.setattr(cli.lib, "generate_pair_entity_tables", lambda t, _m: t)
    monkeypatch.setattr(cli, "analyze_table_quality", lambda *_a, **_k: None)

    buf = io.StringIO()
    configure_logger(LoggerConfig(stream=buf))
    args = argparse.Namespace(
        same_doc=same_doc,
        all_doc=all_doc,
        out_dir=out_dir,
        format="csv",
        dictionary_dir=tmp_path,
    )
    cfg = Config.model_validate({"api": {"user_agent": "test@example.org"}})
    assert cli.run(cfg, args) == 0

    out_file = out_dir / "independent" / "system_independent.csv"
    assert out_file.exists()
    df = pd.read_csv(out_file)
    assert "system_id" in df.columns
    assert "Filtered" in df.columns

    logs = [json.loads(line) for line in buf.getvalue().splitlines() if line]
    assert any(
        "base table 'system_independent' missing" in rec.get("msg", "") for rec in logs
    )
