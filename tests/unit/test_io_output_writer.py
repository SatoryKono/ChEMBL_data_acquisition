"""Unit tests for :mod:`library.io.output_writer`."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from library.config import IoCfg
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
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, _sample_frames
) -> None:
    dataset, correlation, quality = _sample_frames

    captured: dict[str, pd.DataFrame] = {}

    def _fake_write_csv(
        frame: pd.DataFrame,
        destination: Path,
        *,
        cfg: IoCfg,
        **_: object,
    ) -> Path:
        assert cfg.output_dir == tmp_path
        path = Path(destination)
        path.parent.mkdir(parents=True, exist_ok=True)
        frame.to_csv(path, index=False, encoding=cfg.csv_encoding, sep=cfg.csv_sep)
        captured[path.name] = frame.copy()
        return path

    monkeypatch.setattr(output_writer, "write_csv_deterministic", _fake_write_csv)

    cfg = IoCfg(output_dir=tmp_path, csv_sep=",", csv_encoding="utf-8", exist_ok=True)

    artefacts = output_writer.save_standard_outputs(
        dataset,
        correlation,
        quality,
        table_name="documents",
        date_tag="20240101",
        cfg=cfg,
    )

    expected_names = {
        "dataset": "output.documents_20240101.csv",
        "correlation": "_output.documents_20240101_data_correlation_report.csv",
        "quality": "_output.documents_20240101_quality_report.csv",
    }

    assert set(captured) == set(expected_names.values())
    assert artefacts.dataset.name == expected_names["dataset"]
    assert artefacts.correlation_report.name == expected_names["correlation"]
    assert artefacts.quality_report.name == expected_names["quality"]

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
def test_save_standard_outputs__fails_when_directory_missing(tmp_path: Path) -> None:
    cfg = IoCfg(output_dir=tmp_path / "missing", exist_ok=False)

    dataset = pd.DataFrame({"id": [1]})
    correlation = pd.DataFrame({"id": [1.0]})
    quality = pd.DataFrame({"column": ["id"]})

    with pytest.raises(FileNotFoundError):
        output_writer.save_standard_outputs(
            dataset,
            correlation,
            quality,
            table_name="demo",
            date_tag="20240101",
            cfg=cfg,
        )
