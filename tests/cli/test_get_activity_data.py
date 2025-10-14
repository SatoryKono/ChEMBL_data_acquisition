"""Smoke tests for the simplified activity CLI in :mod:`scripts.get_activity_data`."""

from __future__ import annotations

from types import SimpleNamespace

import pandas as pd
import pytest

from library.io.output_writer import StandardOutputArtifacts
from scripts import get_activity_data


def test_parse_activity_args_defaults() -> None:
    namespace, args = get_activity_data.parse_activity_args([])

    assert namespace.limit == 1000
    assert args.limit == 1000
    assert args.output_dir == get_activity_data.Path("data/output")
    assert namespace.date_tag == args.date_tag


def test_run_activity_cli_smoke(monkeypatch: pytest.MonkeyPatch, tmp_path) -> None:
    dataset = pd.DataFrame(
        {
            "activity_id": pd.Series([1], dtype="Int64"),
            "assay_chembl_id": pd.Series(["A1"], dtype=pd.StringDtype()),
            "target_chembl_id": pd.Series(["T1"], dtype=pd.StringDtype()),
            "standard_type": pd.Series(["IC50"], dtype=pd.StringDtype()),
            "standard_relation": pd.Series(["="], dtype=pd.StringDtype()),
            "standard_value": pd.Series([1.0], dtype="Float64"),
            "standard_units": pd.Series(["NM"], dtype=pd.StringDtype()),
        }
    )

    saved_outputs: dict[str, object] = {}

    def _fake_fetch(limit: int) -> pd.DataFrame:
        saved_outputs["fetch_limit"] = limit
        return dataset

    def _fake_generate(frame: pd.DataFrame) -> SimpleNamespace:
        assert frame is dataset
        return SimpleNamespace(
            dataset=frame,
            correlation_report=pd.DataFrame(),
            quality_report=pd.DataFrame(),
            qc_summary={"row_count": 1, "column_count": 7, "non_null_ratio": 1.0},
        )

    def _fake_save_outputs(*args, **kwargs) -> StandardOutputArtifacts:
        saved_outputs["save_args"] = (args, kwargs)
        base = tmp_path / "output.activity_20251019"
        return StandardOutputArtifacts(
            dataset=base.with_suffix(".csv"),
            correlation_report=base.with_name(base.name + "_data_correlation_report_table.csv"),
            quality_report=base.with_name(base.name + "_quality_report_table.csv"),
        )

    def _fake_save_metadata(*args, **kwargs) -> None:
        saved_outputs["metadata"] = (args, kwargs)

    def _fake_load_config(path=None):
        return SimpleNamespace(path=get_activity_data.Path(path or "config/config.yaml"))

    monkeypatch.setattr(get_activity_data, "fetch_normalize_activity", _fake_fetch)
    monkeypatch.setattr(get_activity_data, "generate_activity_reports", _fake_generate)
    monkeypatch.setattr(get_activity_data, "save_standard_outputs", _fake_save_outputs)
    monkeypatch.setattr(get_activity_data, "save_metadata", _fake_save_metadata)
    monkeypatch.setattr(get_activity_data, "load_config", _fake_load_config)

    exit_code = get_activity_data.run_activity_cli(
        [
            "--limit",
            "5",
            "--date-tag",
            "20251019",
            "--output-dir",
            str(tmp_path),
            "--config",
            str(tmp_path / "config.yaml"),
        ]
    )

    assert exit_code == 0
    assert saved_outputs["fetch_limit"] == 5
    args, kwargs = saved_outputs["save_args"]
    assert args[0].equals(dataset)
    assert kwargs["key_columns"] == ["activity_id"]
    assert "metadata" in saved_outputs

