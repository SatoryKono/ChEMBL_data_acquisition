from __future__ import annotations

import math
from pathlib import Path

import pandas as pd
import pytest
from library import io
from library.common.csv_utils import sha256_file
from library.pipelines.testitem import cli
from library.pipelines.testitem.catalog import ParentLookupStats


class _StatsSupplier:
    """Callable wrapper returning a predetermined :class:`ParentLookupStats`."""

    def __init__(self, stats: ParentLookupStats) -> None:
        self._stats = stats
        self.calls = 0

    def __call__(self) -> ParentLookupStats:
        self.calls += 1
        return self._stats


def _base_stats() -> ParentLookupStats:
    return ParentLookupStats(
        source="lookup",
        missing=0,
        unique=2,
        attached=2,
        uncovered=0,
        failed_ids=(),
        hierarchy_attached=1,
        fallback_attached=1,
        no_parent=0,
    )


def _unwrap_finalize_result(
    result: int | tuple[int, io.StandardOutputArtifacts | None]
) -> tuple[int, io.StandardOutputArtifacts | None]:
    if isinstance(result, tuple):
        return result
    return result, None


@pytest.mark.integration
def test_finalize_output__writes_csv_and_metadata(
    tmp_path: Path, sample_input_csv: Path, cfg
) -> None:
    cfg.system.doc_quality.enable = False

    chunk = pd.DataFrame(
        {
            "molecule_chembl_id": pd.Series([" CHEMBL1 ", "CHEMBL2"], dtype="string"),
            "parent_molecule_chembl_id": pd.Series(["CHEMBL10", pd.NA], dtype="string"),
            "natural_product": pd.Series([True, False], dtype="boolean"),
            "salt_chembl_id": pd.Series(["CHEMBL1", pd.NA], dtype="string"),
            "pubchem_cid": pd.Series([123, 456], dtype="Int64"),
        }
    )
    output_path = tmp_path / "final.csv"
    stats_supplier = _StatsSupplier(_base_stats())

    result = cli.finalize_output(
        [chunk],
        cfg=cfg,
        output=output_path,
        parent_stats_supplier=stats_supplier,
        input_csv=sample_input_csv,
        missing_ids=["CHEMBL999"],
    )

    exit_code, artifacts = _unwrap_finalize_result(result)
    assert exit_code == 0
    dataset_path = artifacts.dataset if artifacts is not None else output_path
    assert dataset_path.exists()

    final = pd.read_csv(dataset_path)
    assert list(final.columns[:3]) == [
        "molecule_chembl_id",
        "parent_molecule_chembl_id",
        "salt_chembl_id",
    ]
    assert list(final["molecule_chembl_id"]) == ["CHEMBL1", "CHEMBL2"]
    assert "pipeline_version" in final.columns
    assert "timestamp_utc" in final.columns
    assert sha256_file(dataset_path)
    assert stats_supplier.calls == 1


@pytest.mark.integration
def test_finalize_output__normalises_hidden_output_path(
    tmp_path: Path, sample_input_csv: Path, cfg
) -> None:
    cfg.system.doc_quality.enable = False

    chunk = pd.DataFrame(
        {
            "molecule_chembl_id": pd.Series(["CHEMBL1"], dtype="string"),
            "natural_product": pd.Series([True], dtype="boolean"),
        }
    )
    working_output = tmp_path / ".output.testitems_20240101.csv.tmp"
    stats_supplier = _StatsSupplier(_base_stats())

    exit_code, artifacts = cli.finalize_output(
        [chunk],
        cfg=cfg,
        output=working_output,
        parent_stats_supplier=stats_supplier,
        input_csv=sample_input_csv,
    )

    assert exit_code == 0
    assert artifacts is not None
    assert artifacts.dataset.name == "output.testitem_20240101.csv"
    assert artifacts.dataset.parent == working_output.parent
    assert artifacts.dataset.exists()
    assert not working_output.exists()


@pytest.mark.integration
def test_finalize_output__missing_required_columns_fails(
    tmp_path: Path, sample_input_csv: Path, cfg, monkeypatch: pytest.MonkeyPatch
) -> None:
    cfg.system.doc_quality.enable = False

    chunk = pd.DataFrame({"other": [1]})
    output_path = tmp_path / "invalid.csv"
    warnings: list[tuple[str, dict[str, object]]] = []

    def capture_warning(event: str, **fields: object) -> None:
        warnings.append((event, fields))

    monkeypatch.setattr(cli.logger, "warning", capture_warning)
    stats_supplier = _StatsSupplier(_base_stats())

    result = cli.finalize_output(
        [chunk],
        cfg=cfg,
        output=output_path,
        parent_stats_supplier=stats_supplier,
        input_csv=sample_input_csv,
        emit_legacy_artifacts=True,
    )

    exit_code, _ = _unwrap_finalize_result(result)
    assert exit_code == 1
    assert (
        "validation_skipped",
        {"missing_columns": ["molecule_chembl_id"]},
    ) in warnings


@pytest.mark.integration
def test_finalize_output__optional_columns_missing_warns(
    tmp_path: Path, sample_input_csv: Path, cfg, monkeypatch: pytest.MonkeyPatch
) -> None:
    cfg.system.doc_quality.enable = False

    chunk = pd.DataFrame({"molecule_chembl_id": ["CHEMBL1"]})
    output_path = tmp_path / "optional.csv"
    warnings: list[tuple[str, dict[str, object]]] = []

    def capture_warning(event: str, **fields: object) -> None:
        warnings.append((event, fields))

    monkeypatch.setattr(cli.logger, "warning", capture_warning)
    stats_supplier = _StatsSupplier(_base_stats())

    result = cli.finalize_output(
        [chunk],
        cfg=cfg,
        output=output_path,
        parent_stats_supplier=stats_supplier,
        input_csv=sample_input_csv,
        emit_legacy_artifacts=True,
    )

    exit_code, _ = _unwrap_finalize_result(result)
    assert exit_code == 0
    assert any(event == "optional_columns_missing" for event, _ in warnings)


@pytest.mark.integration
def test_finalize_output__omits_pubchem_columns_when_disabled(
    tmp_path: Path, sample_input_csv: Path, cfg
) -> None:
    cfg.system.doc_quality.enable = False
    cfg.pubchem.enable = False

    chunk = pd.DataFrame({"molecule_chembl_id": ["CHEMBL1"]})
    output_path = tmp_path / "no_pubchem.csv"
    stats_supplier = _StatsSupplier(_base_stats())

    result = cli.finalize_output(
        [chunk],
        cfg=cfg,
        output=output_path,
        parent_stats_supplier=stats_supplier,
        input_csv=sample_input_csv,
        emit_legacy_artifacts=True,
    )

    exit_code, artifacts = _unwrap_finalize_result(result)
    assert exit_code == 0

    dataset_path = artifacts.dataset if artifacts is not None else output_path
    final = pd.read_csv(dataset_path)
    for column in cli._PUBCHEM_OPTIONAL_COLUMNS:
        assert column not in final.columns


@pytest.mark.integration
def test_finalize_output__omits_salt_column_when_enrichment_disabled(
    tmp_path: Path, sample_input_csv: Path, cfg
) -> None:
    cfg.system.doc_quality.enable = False
    cfg.testitem_molecule_enrichment.enable = False

    chunk = pd.DataFrame({"molecule_chembl_id": ["CHEMBL1"]})
    output_path = tmp_path / "no_salt.csv"
    stats_supplier = _StatsSupplier(_base_stats())

    result = cli.finalize_output(
        [chunk],
        cfg=cfg,
        output=output_path,
        parent_stats_supplier=stats_supplier,
        input_csv=sample_input_csv,
        emit_legacy_artifacts=True,
    )

    exit_code, artifacts = _unwrap_finalize_result(result)
    assert exit_code == 0

    dataset_path = artifacts.dataset if artifacts is not None else output_path
    final = pd.read_csv(dataset_path)
    assert cli._SALT_OPTIONAL_COLUMN not in final.columns


@pytest.mark.integration
def test_finalize_output__aligns_nullable_numeric_columns(
    tmp_path: Path, sample_input_csv: Path, cfg
) -> None:
    cfg.system.doc_quality.enable = False

    first_chunk = pd.DataFrame(
        {
            "molecule_chembl_id": pd.Series(["CHEMBL1"], dtype="string"),
            "first_approval": pd.Series([1999.0], dtype="Float64"),
        }
    )
    second_chunk = pd.DataFrame(
        {
            "molecule_chembl_id": pd.Series(["CHEMBL2"], dtype="string"),
        }
    )
    output_path = tmp_path / "numeric_alignment.csv"
    stats_supplier = _StatsSupplier(_base_stats())

    result = cli.finalize_output(
        [first_chunk, second_chunk],
        cfg=cfg,
        output=output_path,
        parent_stats_supplier=stats_supplier,
        input_csv=sample_input_csv,
    )

    exit_code, artifacts = _unwrap_finalize_result(result)
    assert exit_code == 0

    dataset_path = artifacts.dataset if artifacts is not None else output_path
    final = pd.read_csv(dataset_path)
    assert "first_approval" in final.columns
    assert final["first_approval"].iloc[0] == pytest.approx(1999.0)
    assert math.isnan(final["first_approval"].iloc[1])


@pytest.mark.integration
def test_finalize_output__writes_failure_cases(
    tmp_path: Path, sample_input_csv: Path, cfg
) -> None:
    cfg.system.doc_quality.enable = False

    chunk = pd.DataFrame(
        {
            "molecule_chembl_id": ["CHEMBL1"],
            "natural_product": ["unexpected"],
        }
    )
    output_path = tmp_path / "failures.csv"
    stats_supplier = _StatsSupplier(_base_stats())

    result = cli.finalize_output(
        [chunk],
        cfg=cfg,
        output=output_path,
        parent_stats_supplier=stats_supplier,
        input_csv=sample_input_csv,
        emit_legacy_artifacts=True,
    )

    exit_code, _ = _unwrap_finalize_result(result)
    assert exit_code == 1
    failure_path = tmp_path / "failures_failure_cases.csv"
    assert failure_path.exists()
    failure_contents = pd.read_csv(failure_path)
    assert not failure_contents.empty
    assert "natural_product" in failure_contents["column"].tolist()


@pytest.mark.integration
def test_finalize_output__quality_report_failure_non_fatal(
    tmp_path: Path,
    sample_input_csv: Path,
    cfg,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    cfg.system.doc_quality.enable = True
    cfg.system.doc_quality.fatal_on_error = False

    chunk = pd.DataFrame({"molecule_chembl_id": ["CHEMBL1"]})
    output_path = tmp_path / "quality.csv"
    stats_supplier = _StatsSupplier(_base_stats())
    recorded: list[dict[str, object]] = []

    def fake_build_quality(*_args: object, **_kwargs: object):
        def _hook(*__args: object, **__kwargs: object) -> None:
            raise RuntimeError("quality failed")

        return _hook

    def capture_failure(*_args: object, **kwargs: object) -> None:
        recorded.append(kwargs)

    monkeypatch.setattr(cli, "build_table_quality_hook", fake_build_quality)
    monkeypatch.setattr(cli, "record_quality_failure", capture_failure)

    result = cli.finalize_output(
        [chunk],
        cfg=cfg,
        output=output_path,
        parent_stats_supplier=stats_supplier,
        input_csv=sample_input_csv,
        emit_legacy_artifacts=True,
    )

    exit_code, _ = _unwrap_finalize_result(result)
    assert exit_code == 0
    assert recorded and recorded[0]["error"] == "quality failed"


@pytest.mark.integration
def test_finalize_output__quality_report_failure_fatal(
    tmp_path: Path,
    sample_input_csv: Path,
    cfg,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    cfg.system.doc_quality.enable = True
    cfg.system.doc_quality.fatal_on_error = True

    chunk = pd.DataFrame({"molecule_chembl_id": ["CHEMBL1"]})
    output_path = tmp_path / "quality_fatal.csv"
    stats_supplier = _StatsSupplier(_base_stats())

    def fake_build_quality(*_args: object, **_kwargs: object):
        def _hook(*__args: object, **__kwargs: object) -> None:
            raise RuntimeError("quality failed")

        return _hook

    monkeypatch.setattr(cli, "build_table_quality_hook", fake_build_quality)

    result = cli.finalize_output(
        [chunk],
        cfg=cfg,
        output=output_path,
        parent_stats_supplier=stats_supplier,
        input_csv=sample_input_csv,
        emit_legacy_artifacts=True,
    )

    exit_code, _ = _unwrap_finalize_result(result)
    assert exit_code == 1


@pytest.mark.integration
def test_finalize_output__empty_input_produces_placeholder(
    tmp_path: Path, sample_input_csv: Path, cfg
) -> None:
    cfg.system.doc_quality.enable = False

    output_path = tmp_path / "empty.csv"
    stats_supplier = _StatsSupplier(_base_stats())

    result = cli.finalize_output(
        [],
        cfg=cfg,
        output=output_path,
        parent_stats_supplier=stats_supplier,
        input_csv=sample_input_csv,
    )

    exit_code, artifacts = _unwrap_finalize_result(result)
    assert exit_code == 0
    dataset_path = artifacts.dataset if artifacts is not None else output_path
    frame = pd.read_csv(dataset_path)
    assert list(frame.columns[:1]) == ["molecule_chembl_id"]
    assert frame.empty


@pytest.mark.integration
def test_finalize_output__idempotent_results(
    tmp_path: Path, sample_input_csv: Path, cfg
) -> None:
    cfg.system.doc_quality.enable = False

    chunk = pd.DataFrame({"molecule_chembl_id": ["CHEMBL1", "CHEMBL2"]})
    output_path = tmp_path / "stable.csv"
    stats_supplier = _StatsSupplier(_base_stats())

    first_result = cli.finalize_output(
        [chunk],
        cfg=cfg,
        output=output_path,
        parent_stats_supplier=stats_supplier,
        input_csv=sample_input_csv,
    )
    first_exit, first_artifacts = _unwrap_finalize_result(first_result)
    first_dataset = (
        first_artifacts.dataset if first_artifacts is not None else output_path
    )
    first_hash = sha256_file(first_dataset)

    second_result = cli.finalize_output(
        [chunk.copy()],
        cfg=cfg,
        output=output_path,
        parent_stats_supplier=stats_supplier,
        input_csv=sample_input_csv,
    )
    second_exit, second_artifacts = _unwrap_finalize_result(second_result)
    second_dataset = (
        second_artifacts.dataset if second_artifacts is not None else output_path
    )
    second_hash = sha256_file(second_dataset)

    assert first_exit == 0 == second_exit
    assert first_dataset == second_dataset
    assert first_hash == second_hash
