from __future__ import annotations

import argparse
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path

import yaml

from library.common.run_context import RunContext
from library.io import metadata_writer


@dataclass
class _SampleParams:
    identifier: str
    location: Path


def test_save_metadata_serialises_complex_arguments(tmp_path: Path) -> None:
    output_dir = tmp_path / "meta"
    args = argparse.Namespace(
        params=_SampleParams(identifier="step-1", location=tmp_path / "data"),
        options={"threshold": 0.5, "flags": [True, False]},
        executed_at=datetime(2024, 1, 5, 12, 30, tzinfo=UTC),
    )

    artifacts = [
        tmp_path / "output.documents_20240105.csv",
        tmp_path / "output.documents_20240105_quality_report_table.csv",
    ]

    meta_path = metadata_writer.save_metadata(
        table_name="documents",
        date_tag="20240105",
        args=args,
        qc_summary={"rows": 2, "columns": 3},
        artifacts=artifacts,
        output_dir=output_dir,
        sources=["chembl", "pubmed"],
        stats_extra={"records": 2},
        run_context=RunContext(run_id="test", generated_at="2024-01-05T12:34:56Z"),
    )

    data = yaml.safe_load(meta_path.read_text(encoding="utf-8"))

    assert data["table"] == "documents"
    assert data["outputs"] == [path.name for path in artifacts]
    assert data["qc_summary"] == {"rows": 2, "columns": 3}
    assert data["stats"] == {"records": 2}
    assert data["parameters"]["params"]["identifier"] == "step-1"
    assert data["parameters"]["params"]["location"] == str(tmp_path / "data")
    assert data["parameters"]["options"]["flags"] == [True, False]
    assert data["parameters"]["executed_at"] == "2024-01-05T12:30:00Z"
    assert data["generated_at"] == "2024-01-05T12:34:56Z"


def test_save_metadata_uses_default_outputs_when_missing(tmp_path: Path) -> None:
    meta_path = metadata_writer.save_metadata(
        table_name="targets",
        date_tag="20240102",
        args=None,
        output_dir=tmp_path,
        artifacts=None,
        sources=None,
    )

    data = yaml.safe_load(meta_path.read_text(encoding="utf-8"))
    assert data["outputs"] == [
        "output.targets_20240102.csv",
        "output.targets_20240102_quality_report_table.csv",
        "output.targets_20240102_data_correlation_report_table.csv",
    ]
