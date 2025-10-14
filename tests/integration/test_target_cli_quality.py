from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest
from tests.helpers import ensure_normalizer_stub

from library.config import Config

ensure_normalizer_stub()

from scripts import get_target_data


def _install_target_stubs(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    monkeypatch.setattr(get_target_data, "normalize_targets", lambda df: df)
    monkeypatch.setattr(get_target_data, "add_pipeline_metadata", lambda df: df)
    monkeypatch.setattr(
        get_target_data,
        "_prepare_targets_for_schema",
        lambda df: (df, set(), set()),
    )

    class _ValidationResult:
        def __init__(self, data: pd.DataFrame) -> None:
            self.data = data
            self.failure_cases = pd.DataFrame()

    monkeypatch.setattr(
        get_target_data,
        "validate_targets",
        lambda df, return_result=True: _ValidationResult(df),
    )
    monkeypatch.setattr(
        get_target_data,
        "run_target_postprocess_if_requested",
        lambda *_, **__: None,
    )
    monkeypatch.setattr(
        get_target_data,
        "build_table_quality_hook",
        lambda *_, **__: lambda _df: None,
    )

    def _fake_save_standard_outputs(
        frame: pd.DataFrame,
        correlation: pd.DataFrame,
        quality: pd.DataFrame,
        *,
        table_name: str,
        date_tag: str,
        output_path: Path,
        key_columns,
    ) -> SimpleNamespace:
        del correlation, quality, table_name, date_tag, key_columns
        output_path.parent.mkdir(parents=True, exist_ok=True)
        frame.to_csv(output_path, index=False)
        quality_path = output_path.parent / f"{output_path.stem}_quality.csv"
        correlation_path = (
            output_path.parent / f"{output_path.stem}_correlation.csv"
        )
        quality_path.write_text("quality", encoding="utf-8")
        correlation_path.write_text("correlation", encoding="utf-8")
        return SimpleNamespace(
            dataset=output_path,
            quality_report=quality_path,
            correlation_report=correlation_path,
        )

    monkeypatch.setattr(
        get_target_data.io,
        "save_standard_outputs",
        _fake_save_standard_outputs,
    )
    monkeypatch.setattr(get_target_data.io, "save_metadata", lambda *_, **__: None)


@pytest.mark.integration
def test_get_target_cli__doc_quality_parameters(
    tmp_path: Path, cfg: Config, monkeypatch: pytest.MonkeyPatch
) -> None:
    cfg.io.output_dir = tmp_path
    _install_target_stubs(monkeypatch, tmp_path)

    doc_cfg = cfg.system.doc_quality
    doc_cfg.include_columns = ("sequence_length", "complexity")
    doc_cfg.exclude_columns = ("pref_name",)
    doc_cfg.sample_rows = 11
    doc_cfg.correlation_max_columns = 4

    frame = pd.DataFrame(
        {
            "target_chembl_id": ["CHEMBL1", "CHEMBL2"],
            "pref_name": ["Alpha", "Beta"],
            "sequence_length": [123, 456],
            "complexity": [1.5, 2.5],
        }
    )

    captured: dict[str, dict[str, object]] = {}

    def _fake_generate_qc_report(
        data: pd.DataFrame,
        *,
        table_name: str,
        include_columns,
        exclude_columns,
        sample_rows,
        method: str = "pearson",
        profiler=None,
    ) -> pd.DataFrame:
        del data, table_name, method, profiler
        captured["qc"] = {
            "include_columns": include_columns,
            "exclude_columns": exclude_columns,
            "sample_rows": sample_rows,
        }
        return pd.DataFrame({"metric": [1.0]})

    def _fake_generate_correlation_report(
        data: pd.DataFrame,
        *,
        table_name: str,
        include_columns,
        exclude_columns,
        sample_rows,
        method: str = "pearson",
        profiler=None,
    ) -> pd.DataFrame:
        del data, table_name, method, profiler
        captured["correlation"] = {
            "include_columns": include_columns,
            "exclude_columns": exclude_columns,
            "sample_rows": sample_rows,
        }
        return pd.DataFrame({"metric": [0.0]})

    monkeypatch.setattr(
        get_target_data, "generate_qc_report", _fake_generate_qc_report
    )
    monkeypatch.setattr(
        get_target_data,
        "generate_correlation_report",
        _fake_generate_correlation_report,
    )

    output_path = tmp_path / "output.targets_20240101.csv"

    exit_code = get_target_data.validate_and_write(
        frame,
        output_path,
        cfg,
        emit_legacy_artifacts=False,
    )

    assert exit_code == 0
    assert "qc" in captured and "correlation" in captured
    qc_call = captured["qc"]
    corr_call = captured["correlation"]
    assert qc_call["include_columns"] == doc_cfg.include_columns
    assert qc_call["exclude_columns"] == doc_cfg.exclude_columns
    assert qc_call["sample_rows"] == doc_cfg.sample_rows
    assert corr_call["include_columns"] == doc_cfg.include_columns
    assert corr_call["exclude_columns"] == doc_cfg.exclude_columns
    assert corr_call["sample_rows"] == doc_cfg.sample_rows


@pytest.mark.integration
def test_get_target_cli__correlation_column_sampling(
    tmp_path: Path, cfg: Config, monkeypatch: pytest.MonkeyPatch
) -> None:
    cfg.io.output_dir = tmp_path
    _install_target_stubs(monkeypatch, tmp_path)

    doc_cfg = cfg.system.doc_quality
    doc_cfg.include_columns = None
    doc_cfg.exclude_columns = ("pref_name",)
    doc_cfg.sample_rows = None
    doc_cfg.correlation_max_columns = 2

    frame = pd.DataFrame(
        {
            "target_chembl_id": ["CHEMBL1", "CHEMBL2"],
            "pref_name": ["Alpha", "Beta"],
            "sequence_length": [123, 456],
            "complexity": [1.5, 2.5],
            "confidence_score": [0.8, 0.9],
        }
    )

    captured: dict[str, dict[str, object]] = {}

    def _fake_generate_qc_report(
        data: pd.DataFrame,
        *,
        table_name: str,
        include_columns,
        exclude_columns,
        sample_rows,
        method: str = "pearson",
        profiler=None,
    ) -> pd.DataFrame:
        del data, table_name, method, profiler
        captured["qc"] = {
            "include_columns": include_columns,
            "exclude_columns": exclude_columns,
            "sample_rows": sample_rows,
        }
        return pd.DataFrame({"metric": [1.0]})

    def _fake_generate_correlation_report(
        data: pd.DataFrame,
        *,
        table_name: str,
        include_columns,
        exclude_columns,
        sample_rows,
        method: str = "pearson",
        profiler=None,
    ) -> pd.DataFrame:
        del data, table_name, method, profiler
        captured["correlation"] = {
            "include_columns": include_columns,
            "exclude_columns": exclude_columns,
            "sample_rows": sample_rows,
        }
        return pd.DataFrame({"metric": [0.0]})

    monkeypatch.setattr(
        get_target_data, "generate_qc_report", _fake_generate_qc_report
    )
    monkeypatch.setattr(
        get_target_data,
        "generate_correlation_report",
        _fake_generate_correlation_report,
    )

    output_path = tmp_path / "output.targets_20240202.csv"

    exit_code = get_target_data.validate_and_write(
        frame,
        output_path,
        cfg,
        emit_legacy_artifacts=False,
    )

    assert exit_code == 0
    assert "correlation" in captured
    corr_call = captured["correlation"]
    assert corr_call["include_columns"] == ("sequence_length", "complexity")
    assert corr_call["exclude_columns"] == doc_cfg.exclude_columns
    assert corr_call["sample_rows"] is None
    qc_call = captured["qc"]
    assert qc_call["include_columns"] is None
    assert qc_call["exclude_columns"] == doc_cfg.exclude_columns
    assert qc_call["sample_rows"] is None
