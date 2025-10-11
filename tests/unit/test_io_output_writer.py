"""Unit tests for :mod:`library.io.output_writer`."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from library.io import output_writer


@pytest.fixture()
def _sample_frames() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    dataset = pd.DataFrame(
        {
            "identifier": ["row-1", "row-2"],
            "value": [1, 2],
        }
    )
    correlation = pd.DataFrame(
        {
            "identifier": [1.0, 0.0],
            "value": [0.0, 1.0],
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
    tmp_path: Path, _sample_frames
) -> None:
    dataset, correlation, quality = _sample_frames

    captured: dict[str, pd.DataFrame] = {}

    def _fake_write_csv(
        frame: pd.DataFrame,
        destination: Path,
        **_: object,
    ) -> Path:
        path = Path(destination)
        path.parent.mkdir(parents=True, exist_ok=True)
        frame.to_csv(path, index=False)
        captured[path.name] = frame.copy()
        return path

    monkeypatch = pytest.MonkeyPatch()
    monkeypatch.setattr(output_writer, "write_csv_deterministic", _fake_write_csv)

    artefacts = output_writer.save_standard_outputs(
        dataset,
        correlation,
        quality,
        table_name="documents",
        date_tag="20240101",
        output_dir=tmp_path,
    )
    monkeypatch.undo()

    expected_names = {
        "dataset": "output.documents_20240101.csv",
        "correlation": "output.documents_20240101_data_correlation_report_table.csv",
        "quality": "output.documents_20240101_quality_report_table.csv",
    }

    assert set(captured) == set(expected_names.values())
    assert artefacts.dataset.name == expected_names["dataset"]
    assert artefacts.correlation_report.name == expected_names["correlation"]
    assert artefacts.quality_report.name == expected_names["quality"]
    assert artefacts.dataset.parent == tmp_path
    assert artefacts.correlation_report.parent == tmp_path
    assert artefacts.quality_report.parent == tmp_path

    output_files = sorted(p.name for p in tmp_path.glob("*.csv"))
    assert set(output_files) == set(expected_names.values())

    pd.testing.assert_frame_equal(
        pd.read_csv(artefacts.dataset), dataset
    )
    pd.testing.assert_frame_equal(
        pd.read_csv(artefacts.correlation_report), correlation
    )
    pd.testing.assert_frame_equal(
        pd.read_csv(artefacts.quality_report), quality
    )


@pytest.mark.unit
def test_save_standard_outputs__creates_output_directory(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    dataset = pd.DataFrame({"id": [1]})
    correlation = pd.DataFrame({"id": [1.0]})
    quality = pd.DataFrame({"column": ["id"]})

    output_dir = tmp_path / "missing" / "nested"
    output_writer.save_standard_outputs(
        dataset,
        correlation,
        quality,
        table_name="demo",
        date_tag="20240101",
        output_dir=output_dir,
    )

    assert output_dir.exists() and output_dir.is_dir()


@pytest.mark.unit
@pytest.mark.unit
def test_save_standard_outputs__supports_string_base_directory(
    tmp_path: Path, _sample_frames
) -> None:
    dataset, correlation, quality = _sample_frames
    base = tmp_path / "custom"

    artifacts = output_writer.save_standard_outputs(
        dataset,
        correlation,
        quality,
        table_name="documents",
        date_tag="20240101",
        output_dir=str(base),
    )

    assert artifacts.dataset.parent == base
    assert artifacts.dataset.exists()


@pytest.mark.unit
def test_save_standard_outputs__requires_output_directory(
    _sample_frames: tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]
) -> None:
    dataset, correlation, quality = _sample_frames

    with pytest.raises(ValueError, match="output_dir must be provided"):
        output_writer.save_standard_outputs(
            dataset,
            correlation,
            quality,
            table_name="documents",
            date_tag="20240101",
        )


@pytest.mark.unit
def test_save_standard_outputs__uses_explicit_output_path(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, _sample_frames
) -> None:
    dataset, correlation, quality = _sample_frames

    # Ensure the fallback directory would differ from our explicit path.
    fallback_dir = tmp_path / "fallback"
    monkeypatch.setattr(output_writer, "OUTPUT_DIR", fallback_dir)

    output_path = tmp_path / "final" / ".output.targets_20240101.csv.tmp"

    captured: dict[Path, pd.DataFrame] = {}

    def _fake_write_csv(
        frame: pd.DataFrame,
        destination: Path,
        **_: object,
    ) -> Path:
        path = Path(destination)
        path.parent.mkdir(parents=True, exist_ok=True)
        captured[path] = frame.copy()
        return path

    monkeypatch.setattr(output_writer, "write_csv_deterministic", _fake_write_csv)

    artefacts = output_writer.save_standard_outputs(
        dataset,
        correlation,
        quality,
        table_name="targets",
        date_tag="20240101",
        output_path=output_path,
    )

    assert artefacts.dataset == output_path
    expected_stem = output_path.with_suffix("").name
    expected_correlation = output_path.parent / (
        f"{expected_stem}_data_correlation_report_table.csv"
    )
    expected_quality = output_path.parent / (
        f"{expected_stem}_quality_report_table.csv"
    )

    assert artefacts.correlation_report == expected_correlation
    assert artefacts.quality_report == expected_quality

    assert set(captured) == {
        artefacts.dataset,
        artefacts.correlation_report,
        artefacts.quality_report,
    }

    # Fallback directory should remain unused.
    assert not fallback_dir.exists()
