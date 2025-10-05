"""Integration tests for the tissue data pipeline."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from library.pipelines.tissue import TISSUE_COLUMN_ORDER
from library.pipelines.tissue.pipeline import (
    TissuePipelineOptions,
    run_tissue_pipeline,
)


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

