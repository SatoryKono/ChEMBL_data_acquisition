"""Smoke tests for the document CLI entry point."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd
import pytest

from library.cli.commands import get_document_data as document_cli


@pytest.mark.integration
def test_document_cli_writes_standard_outputs(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """The CLI should emit the canonical dataset, QC reports and metadata."""

    input_csv = tmp_path / "documents.csv"
    input_csv.write_text("document_chembl_id\nDOC0001\nDOC0002\n", encoding="utf-8")

    class _StubContext:
        def __init__(self, cfg: Any) -> None:
            self.cfg = cfg
            self.chembl_client = object()

        def __enter__(self) -> "_StubContext":
            return self

        def __exit__(self, exc_type, exc, tb) -> None:
            return None

    def _fake_get_documents(
        ids: Any,
        *,
        cfg: Any,
        client: Any,
        chunk_size: int = 5,
        timeout: float | None = None,
    ) -> pd.DataFrame:
        identifiers = list(ids)
        rows: list[dict[str, object]] = []
        for index, identifier in enumerate(identifiers):
            rows.append(
                {
                    "document_chembl_id": identifier,
                    "title": f"Synthetic Title {index}",
                    "doc_type": "journal",
                    "doi": f"10.5555/{index}",
                    "year": "2024",
                    "journal": "Synthetic Journal",
                }
            )
        return pd.DataFrame(rows)

    def _fake_qc_report(df: pd.DataFrame, table_name: str, profiler: Any) -> pd.DataFrame:
        return pd.DataFrame({"metric": ["rows"], "value": [len(df)]})

    def _fake_corr_report(df: pd.DataFrame, table_name: str, profiler: Any) -> pd.DataFrame:
        return pd.DataFrame({"metric": [], "value": []})

    monkeypatch.setattr(document_cli, "ETLContext", lambda cfg: _StubContext(cfg))
    monkeypatch.setattr(document_cli.cl, "get_documents", _fake_get_documents)
    monkeypatch.setattr(document_cli, "generate_qc_report", _fake_qc_report)
    monkeypatch.setattr(document_cli, "generate_correlation_report", _fake_corr_report)

    output_target = tmp_path / "output.csv"
    args = [
        "--config",
        str(Path("config/config.yaml")),
        "--mode",
        "chembl",
        "--input",
        str(input_csv),
        "--final-out",
        str(output_target),
    ]

    exit_code = document_cli.main(args)
    assert exit_code == 0

    dataset_files = list(tmp_path.glob("output.documents_*.csv"))
    assert len(dataset_files) == 1
    dataset_path = dataset_files[0]
    correlation_path = dataset_path.with_name(f"{dataset_path.stem}_data_correlation_report_table.csv")
    quality_path = dataset_path.with_name(f"{dataset_path.stem}_quality_report_table.csv")
    meta_path = dataset_path.with_name(f"{dataset_path.name}.meta.yaml")

    assert dataset_path.is_file()
    assert correlation_path.is_file()
    assert quality_path.is_file()
    assert meta_path.is_file()

    dataset = pd.read_csv(dataset_path)
    assert set(["document_chembl_id", "title", "doc_type"]).issubset(dataset.columns)
    assert len(dataset) == 2
