"""Integration tests for the tissue data pipeline."""

from __future__ import annotations

import datetime as dt
import json
from pathlib import Path

import pandas as pd
import pytest
import yaml

from library.io import default_output_path
import library.io.metadata as io_metadata
from library.pipelines.tissue import TISSUE_COLUMN_ORDER
from library.pipelines.tissue.pipeline import (
    TissuePipelineOptions,
    run_tissue_pipeline,
)
from library.pipelines.common import metadata as pipeline_metadata_module
from library.common import run_context as run_context_module
from library.common.run_context import RunContext


RUN_GENERATED_AT = "2020-01-01T00:00:00+00:00"


class _DummyClient:
    def close(self) -> None:  # pragma: no cover - interface compatibility
        return None


def test_run_tissue_pipeline__writes_normalised_output(tmp_path: Path, cfg, monkeypatch) -> None:
    """The pipeline writes a deterministic CSV and reports missing identifiers."""

    input_csv = tmp_path / "tissue.csv"
    input_csv.write_text("tissue_chembl_id\nCHEMBLT1\nCHEMBLT2\nCHEMBLT3\n", encoding="utf-8")
    output_csv = tmp_path / "output.csv"

    expected_calls: list[dict[str, object]] = []

    def fake_get_tissues(
        ids,
        *,
        cfg,
        client,
        chunk_size: int,
        timeout: float | None,
    ) -> pd.DataFrame:
        expected_calls.append(
            {
                "ids": list(ids),
                "chunk_size": chunk_size,
                "timeout": timeout,
            }
        )
        return pd.DataFrame(
            {
                "tissue_chembl_id": ["CHEMBLT2", "CHEMBLT1"],
                "pref_name": ["Beta", "Alpha"],
                "uberon_id": [pd.NA, "UBERON:0001"],
                "efo_id": ["EFO:0002", pd.NA],
                "bto_id": [pd.NA, pd.NA],
                "caloha_id": [pd.NA, pd.NA],
            }
        )

    monkeypatch.setattr(
        "library.pipelines.tissue.pipeline.get_tissues",
        fake_get_tissues,
    )

    options = TissuePipelineOptions(
        input_csv=input_csv,
        output_csv=output_csv,
        column="tissue_chembl_id",
        batch_size=2,
        limit=None,
        offset=0,
        timeout=None,
    )

    monkeypatch.setattr(
        run_context_module,
        "_CURRENT",
        RunContext(run_id="tissue-test", generated_at=RUN_GENERATED_AT),
    )

    result = run_tissue_pipeline(cfg, options, client=_DummyClient())

    assert result.exit_code == 0
    assert result.records == 3
    assert result.missing_ids == ("CHEMBLT3",)
    assert result.output_path == output_csv
    assert output_csv.exists()
    assert expected_calls == [
        {
            "ids": ["CHEMBLT1", "CHEMBLT2", "CHEMBLT3"],
            "chunk_size": 2,
            "timeout": cfg.tissue.timeout,
        }
    ]

    output_df = pd.read_csv(output_csv, sep=cfg.io.csv_sep, dtype="string")
    assert list(output_df.columns) == TISSUE_COLUMN_ORDER
    assert output_df["tissue_chembl_id"].tolist() == [
        "CHEMBLT1",
        "CHEMBLT2",
        "CHEMBLT3",
    ]
    assert pd.isna(output_df.loc[2, "pref_name"])
    assert (tmp_path / "output_validation_failures.csv").exists() is False
    meta_path = Path(f"{output_csv}.meta.yaml")
    assert meta_path.exists()
    metadata = yaml.safe_load(meta_path.read_text(encoding="utf-8"))
    assert metadata["generated_at"] == RUN_GENERATED_AT


@pytest.mark.integration
def test_run_tissue_pipeline__deterministic_output_from_fixtures(
    tmp_path: Path,
    cfg,
    snapshot_resource: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Pipeline output matches recorded fixtures and respects deterministic naming."""

    resource_dir = snapshot_resource / "tissue_pipeline"
    payload_path = resource_dir / "integration_payload.json"
    expected_csv = resource_dir / "integration_expected_output.csv"
    payload = json.loads(payload_path.read_text(encoding="utf-8"))

    input_dir = tmp_path / "inputs"
    input_dir.mkdir(parents=True, exist_ok=True)
    input_csv = input_dir / "tissue_ids.csv"
    with input_csv.open("w", encoding="utf-8", newline="") as handle:
        handle.write("tissue_chembl_id\n")
        for identifier in payload["input_ids"]:
            handle.write(f"{identifier}\n")

    cfg.io.output_dir = tmp_path / "exports"
    cfg.io.csv_sep = ","
    cfg.io.csv_encoding = "utf-8"
    cfg.tissue.timeout = 9.25

    output_path = default_output_path(input_csv, cfg.io)

    expected_calls: list[dict[str, object]] = []

    def fake_get_tissues(ids, *, cfg, client, chunk_size: int, timeout: float | None):
        expected_calls.append(
            {
                "ids": list(ids),
                "chunk_size": chunk_size,
                "timeout": timeout,
            }
        )
        frame = pd.DataFrame(payload["records"])
        return frame.iloc[[1, 0, 2]].reset_index(drop=True)

    monkeypatch.setattr(
        "library.pipelines.tissue.pipeline.get_tissues",
        fake_get_tissues,
    )

    options = TissuePipelineOptions(
        input_csv=input_csv,
        output_csv=output_path,
        column="tissue_chembl_id",
        batch_size=2,
        limit=None,
        offset=0,
        timeout=None,
    )

    monkeypatch.setattr(
        pipeline_metadata_module,
        "datetime",
        dt.datetime,
        raising=False,
    )
    monkeypatch.setattr(
        io_metadata,
        "datetime",
        dt.datetime,
        raising=False,
    )
    pipeline_metadata_module.get_timestamp_utc.cache_clear()
    pipeline_metadata_module.pipeline_metadata.cache_clear()

    monkeypatch.setattr(
        run_context_module,
        "_CURRENT",
        RunContext(run_id="tissue-test", generated_at=RUN_GENERATED_AT),
    )

    result = run_tissue_pipeline(cfg, options, client=_DummyClient())

    assert result.exit_code == 0
    assert result.missing_ids == ()
    assert result.records == 3
    assert result.output_path == output_path
    assert result.output_path.name == output_path.name
    assert expected_calls == [
        {
            "ids": payload["input_ids"],
            "chunk_size": 2,
            "timeout": cfg.tissue.timeout,
        }
    ]

    output_df = pd.read_csv(result.output_path, sep=cfg.io.csv_sep, dtype="string")
    expected_df = pd.read_csv(expected_csv, sep=cfg.io.csv_sep, dtype="string")

    pd.testing.assert_frame_equal(output_df, expected_df)
    assert list(output_df.columns) == TISSUE_COLUMN_ORDER
    assert output_df["tissue_chembl_id"].tolist() == sorted(payload["input_ids"])

    meta_path = Path(f"{result.output_path}.meta.yaml")
    assert meta_path.exists()
    metadata = yaml.safe_load(meta_path.read_text(encoding="utf-8"))
    assert metadata["columns"] == TISSUE_COLUMN_ORDER
    assert metadata["generated_at"] == RUN_GENERATED_AT
    assert metadata["command"]

