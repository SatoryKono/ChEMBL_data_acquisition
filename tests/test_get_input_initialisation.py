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
        "assay": pd.DataFrame({"assay_chembl_id": ["a1", "a2"], "b": [1, 2]}),
        "activity": pd.DataFrame(
            {
                "activity_chembl_id": [1, 2, 3],
                "assay_chembl_id": ["a1", "a2", "a1"],
                "document_chembl_id": ["d1", "d2", "d3"],
                "target_chembl_id": ["t1", "t2", "t3"],
                "molecule_chembl_id": ["m1", "m2", "m3"],
            }
        ),
        "document": pd.DataFrame({"document_chembl_id": ["d1", "d2", "d3"]}),
        "target": pd.DataFrame({"target_chembl_id": ["t1", "t2", "t3"]}),
        "testitem": pd.DataFrame({"molecule_chembl_id": ["m1", "m2", "m3"]}),
        "pairs_same_document": pd.DataFrame(
            {"activity_chembl_id1": [1], "activity_chembl_id2": [2]}
        ),
        "pairs_independent": pd.DataFrame(
            {"activity_chembl_id1": [1], "activity_chembl_id2": [3]}
        ),
        "pairs_non_independent": pd.DataFrame(
            {"activity_chembl_id1": [2], "activity_chembl_id2": [3]}
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
    cfg = Config()
    cfg.api.user_agent = "test@example.org"
    result = cli.run(cfg, args)
    assert result == 0

    assert (out_dir / "independent" / "activity_independent.csv").exists()

    data_files = [
        p for p in out_dir.rglob("*.csv") if "data_validity_report" not in p.parts
    ]
    for path in data_files:
        report = (
            out_dir / "data_validity_report" / f"{path.stem}_quality_report_table.csv"
        )

        corr = (
            out_dir
            / "data_validity_report"
            / f"{path.stem}_data_correlation_report_table.csv"
        )
        assert report.exists(), report
        assert corr.exists(), corr
        df = pd.read_csv(report)
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
    cfg = Config()
    cfg.api.user_agent = "test@example.org"
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
    cfg = Config()
    cfg.api.user_agent = "test@example.org"
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

    cfg = Config()
    cfg.api.user_agent = "test@example.org"
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
    cfg_path.write_text('system:\n  log:\n    level: "ERROR"\n')
    same_doc = tmp_path / "same.xlsx"
    all_doc = tmp_path / "all.xlsx"
    same_doc.write_text("dummy")
    all_doc.write_text("dummy")
    out_dir = tmp_path / "out"

    # Skip heavy processing by stubbing ``run``
    monkeypatch.setattr(cli, "run", lambda _cfg, _args: 0)
    monkeypatch.setenv(
        "CHEMBL_DA__SOURCES__CHEMBL__API__USER_AGENT", "test@example.org"
    )

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
