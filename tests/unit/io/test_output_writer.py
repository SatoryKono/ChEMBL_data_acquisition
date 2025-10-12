"""Unit tests for :mod:`library.io.output_writer`."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from library.io import output_writer


@pytest.fixture()
def sample_frames() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    dataset = pd.DataFrame(
        {
            "identifier": ["row-1", "row-2"],
            "value": [1, 2],
        }
    )
    correlation = pd.DataFrame(
        {
            "identifier": ["row-1", "row-2"],
            "value": [0.75, 1.0],
        }
    )
    quality = pd.DataFrame(
        {
            "column": ["identifier", "value"],
            "non_null": [2, 2],
        }
    )
    return dataset, correlation, quality


@pytest.mark.unit
def test_save_standard_outputs__writes_expected_csvs(
    tmp_path: Path, sample_frames, monkeypatch: pytest.MonkeyPatch
) -> None:
    dataset, correlation, quality = sample_frames

    monkeypatch.setattr(output_writer, "OUTPUT_DIR", tmp_path)

    artifacts = output_writer.save_standard_outputs(
        dataset,
        correlation,
        quality,
        table_name="documents",
        date_tag="20240101",
    )

    expected_names = {
        "dataset": "output.documents_20240101.csv",
        "correlation": "output.documents_20240101_data_correlation_report_table.csv",
        "quality": "output.documents_20240101_quality_report_table.csv",
    }

    assert artifacts.dataset == tmp_path / expected_names["dataset"]
    assert artifacts.correlation_report == tmp_path / expected_names["correlation"]
    assert artifacts.quality_report == tmp_path / expected_names["quality"]

    for path in expected_names.values():
        assert (tmp_path / path).exists(), f"missing {path}"

    pd.testing.assert_frame_equal(pd.read_csv(artifacts.dataset), dataset)
    pd.testing.assert_frame_equal(pd.read_csv(artifacts.correlation_report), correlation)
    pd.testing.assert_frame_equal(pd.read_csv(artifacts.quality_report), quality)


@pytest.mark.unit
def test_save_standard_outputs__creates_output_directory(tmp_path: Path, sample_frames) -> None:
    dataset, correlation, quality = sample_frames

    destination = tmp_path / "nested" / "dir"
    artifacts = output_writer.save_standard_outputs(
        dataset,
        correlation,
        quality,
        table_name="assays",
        date_tag="20240101",
        output_dir=destination,
    )

    assert destination.exists() and destination.is_dir()
    assert artifacts.dataset.parent == destination


@pytest.mark.unit
def test_save_standard_outputs__supports_string_output_dir(tmp_path: Path, sample_frames) -> None:
    dataset, correlation, quality = sample_frames
    destination = tmp_path / "string-dir"

    artifacts = output_writer.save_standard_outputs(
        dataset,
        correlation,
        quality,
        table_name="targets",
        date_tag="20240101",
        output_dir=str(destination),
    )

    assert artifacts.dataset.parent == destination
    assert artifacts.dataset.exists()


@pytest.mark.unit
def test_save_standard_outputs__uses_canonical_naming_and_cleans_source(
    tmp_path: Path, sample_frames
) -> None:
    dataset, correlation, quality = sample_frames

    legacy_path = tmp_path / "legacy" / "result.tmp.csv"
    legacy_path.parent.mkdir(parents=True, exist_ok=True)
    legacy_path.write_text("legacy")
    (legacy_path.parent / "result.tmp.csv.meta.yaml").write_text("meta")

    artifacts = output_writer.save_standard_outputs(
        dataset,
        correlation,
        quality,
        table_name="testitem",
        date_tag="20240101",
        output_path=legacy_path,
    )

    expected_dataset = legacy_path.parent / "output.testitem_20240101.csv"
    assert artifacts.dataset == expected_dataset
    assert artifacts.correlation_report == legacy_path.parent / (
        "output.testitem_20240101_data_correlation_report_table.csv"
    )
    assert artifacts.quality_report == legacy_path.parent / (
        "output.testitem_20240101_quality_report_table.csv"
    )

    assert not legacy_path.exists()
    assert not (legacy_path.parent / "result.tmp.csv.meta.yaml").exists()


@pytest.mark.unit
def test_save_standard_outputs__idempotent_and_without_extraneous_files(
    tmp_path: Path, sample_frames, monkeypatch: pytest.MonkeyPatch
) -> None:
    dataset, correlation, quality = sample_frames

    monkeypatch.setattr(output_writer, "OUTPUT_DIR", tmp_path)

    artifacts_first = output_writer.save_standard_outputs(
        dataset,
        correlation,
        quality,
        table_name="targets",
        date_tag="20240101",
    )

    first_paths = {
        "dataset": artifacts_first.dataset,
        "correlation": artifacts_first.correlation_report,
        "quality": artifacts_first.quality_report,
    }
    first_snapshots = {
        name: path.read_bytes() for name, path in first_paths.items()
    }

    artifacts_second = output_writer.save_standard_outputs(
        dataset.copy(),
        correlation.copy(),
        quality.copy(),
        table_name="targets",
        date_tag="20240101",
    )

    assert artifacts_second == artifacts_first

    observed = {child.name for child in tmp_path.iterdir() if child.is_file()}
    expected = {path.name for path in first_paths.values()}
    assert observed == expected, "unexpected artefacts detected"

    for name, path in {
        "dataset": artifacts_second.dataset,
        "correlation": artifacts_second.correlation_report,
        "quality": artifacts_second.quality_report,
    }.items():
        assert path.exists(), f"missing {name}"
        assert path.read_bytes() == first_snapshots[name]


@pytest.mark.unit
def test_set_output_directory__updates_default(tmp_path: Path, sample_frames) -> None:
    dataset, correlation, quality = sample_frames

    original_dir = output_writer.OUTPUT_DIR
    try:
        output_writer.set_output_directory(tmp_path)

        artifacts = output_writer.save_standard_outputs(
            dataset,
            correlation,
            quality,
            table_name="documents",
            date_tag="20240202",
        )

        assert artifacts.dataset.parent == tmp_path
        assert artifacts.correlation_report.parent == tmp_path
        assert artifacts.quality_report.parent == tmp_path
    finally:
        output_writer.set_output_directory(original_dir)
