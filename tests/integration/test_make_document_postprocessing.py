from __future__ import annotations

import json
from datetime import UTC, datetime
from textwrap import dedent

import pandas as pd
from scripts import make_document_postprocessing as cli

from library.postprocessing.documents.schema import DOCUMENT_SCHEMA


def test_make_document_postprocessing__end_to_end(tmp_path, monkeypatch):
    log_dir = tmp_path / "logs"
    monkeypatch.setenv("CHEMBL_POSTPROCESS_LOG_DIR", str(log_dir))

    input_path = tmp_path / "documents_raw.csv"
    df = pd.DataFrame(
        {
            "document_chembl_id": ["CHEMBL_DOC2", "CHEMBL_DOC1"],
            "title": ["  First Document  ", "Second"],
            "doc_type": ["article", "journal"],
            "year": [2020, 2019],
            "journal": ["  Nature ", "Science"],
        }
    )
    df.to_csv(input_path, index=False)

    output_path = tmp_path / "documents_processed.csv"
    exit_code = cli.main(["--input", str(input_path), "--output", str(output_path)])
    assert exit_code == 0

    result = pd.read_csv(output_path)
    expected_columns = [
        col for col in DOCUMENT_SCHEMA.column_order if col in result.columns
    ]
    assert list(result.columns) == expected_columns
    assert result["document_chembl_id"].tolist() == ["CHEMBL_DOC1", "CHEMBL_DOC2"]
    assert result["title"].tolist() == ["Second", "First Document"]
    assert result["doc_type"].tolist() == ["journal", "article"]
    assert result["publication_year"].tolist() == [2019, 2020]
    assert result["year"].tolist() == [2019, 2020]
    assert result["journal"].tolist() == ["Science", "Nature"]
    assert result["pipeline_version"].notna().all()

    report_path = output_path.parent / "documents.postprocess.report.json"
    payload = json.loads(report_path.read_text(encoding="utf-8"))
    assert payload["table"] == "documents"
    assert payload["metrics"]["output"]["rows"] == len(result)
    extras = payload.get("extras")
    assert extras and extras["output_postprocessed"] == str(output_path)

    date_str = datetime.now(UTC).strftime("%Y%m%d")
    log_path = log_dir / f"make_documents_postprocessing_{date_str}.log"
    assert log_path.exists()

    snapshot = output_path.read_bytes()
    second_exit = cli.main(["--input", str(input_path), "--output", str(output_path)])
    assert second_exit == 0
    assert output_path.read_bytes() == snapshot


def test_make_document_postprocessing__fills_missing_title_from_known_sources(
    tmp_path, monkeypatch
):
    log_dir = tmp_path / "logs"
    monkeypatch.setenv("CHEMBL_POSTPROCESS_LOG_DIR", str(log_dir))

    input_path = tmp_path / "documents_raw.csv"
    df = pd.DataFrame(
        {
            "document_chembl_id": ["CHEMBL_DOC3"],
            "ChEMBL.title": ["  Derived Title  "],
            "doc_type": ["article"],
            "year": [2022],
            "journal": [" Example Journal "],
        }
    )
    df.to_csv(input_path, index=False)

    output_path = tmp_path / "documents_processed.csv"
    exit_code = cli.main(["--input", str(input_path), "--output", str(output_path)])
    assert exit_code == 0

    result = pd.read_csv(output_path)
    assert "title" in result.columns
    assert result.loc[0, "title"] == "Derived Title"
    assert result.loc[0, "doc_type"] == "article"
    assert result.loc[0, "journal"] == "Example Journal"


def test_make_document_postprocessing__uses_cli_pipeline_override(
    tmp_path, monkeypatch
):
    log_dir = tmp_path / "logs"
    monkeypatch.setenv("CHEMBL_POSTPROCESS_LOG_DIR", str(log_dir))

    pipeline_path = tmp_path / "documents_override.yaml"
    pipeline_path.write_text(
        dedent(
            """
            pipeline_version: "override-2024"
            enabled_steps:
              - name: normalize_document_fields
                callable: "library.postprocessing.documents.steps:normalize_document_fields"
                params:
                  trim_whitespace: false
                  normalise_unicode: false
              - name: finalize_document_records
                callable: "library.postprocessing.documents.steps:finalize_document_records"
                params:
                  enforce_schema: true
                  ensure_unique_ids: true
            """
        ).strip()
    )

    input_path = tmp_path / "documents_raw.csv"
    df = pd.DataFrame(
        {
            "document_chembl_id": ["CHEMBL_DOC2"],
            "title": ["  First Document  "],
            "doc_type": ["article"],
            "year": [2020],
            "journal": ["  Nature "],
        }
    )
    df.to_csv(input_path, index=False)

    output_path = tmp_path / "documents_processed.csv"
    exit_code = cli.main(
        [
            "--input",
            str(input_path),
            "--output",
            str(output_path),
            "--config",
            str(pipeline_path),
        ]
    )
    assert exit_code == 0

    result = pd.read_csv(output_path)
    assert result.loc[0, "title"] == "  First Document  "
    assert result.loc[0, "journal"] == "  Nature "
    assert result.loc[0, "pipeline_version"] == "override-2024"

    report_path = output_path.parent / "documents.postprocess.report.json"
    payload = json.loads(report_path.read_text(encoding="utf-8"))
    assert payload["metrics"]["pipeline_version"] == "override-2024"
    step_names = [step["name"] for step in payload["metrics"]["steps"]]
    assert step_names == [
        "normalize_document_fields",
        "finalize_document_records",
    ]


def test_make_document_postprocessing__resolves_pipeline_version_from_defaults(
    tmp_path, monkeypatch
):
    log_dir = tmp_path / "logs"
    monkeypatch.setenv("CHEMBL_POSTPROCESS_LOG_DIR", str(log_dir))

    pipeline_path = tmp_path / "documents_defaults.yaml"
    pipeline_path.write_text(
        dedent(
            """
            enabled_steps:
              - name: normalize_document_fields
                callable: "library.postprocessing.documents.steps:normalize_document_fields"
                params:
                  trim_whitespace: false
                  normalise_unicode: false
              - name: finalize_document_records
                callable: "library.postprocessing.documents.steps:finalize_document_records"
                params:
                  enforce_schema: true
                  ensure_unique_ids: true
            params:
              defaults:
                pipeline_version: defaults-override
            """
        ).strip()
    )

    input_path = tmp_path / "documents_raw.csv"
    df = pd.DataFrame(
        {
            "document_chembl_id": ["CHEMBL_DOC2"],
            "title": ["  First Document  "],
            "doc_type": ["article"],
            "year": [2020],
            "journal": ["  Nature "],
        }
    )
    df.to_csv(input_path, index=False)

    output_path = tmp_path / "documents_processed.csv"
    exit_code = cli.main(
        [
            "--input",
            str(input_path),
            "--output",
            str(output_path),
            "--config",
            str(pipeline_path),
        ]
    )
    assert exit_code == 0

    result = pd.read_csv(output_path)
    assert result.loc[0, "pipeline_version"] == "defaults-override"

    report_path = output_path.parent / "documents.postprocess.report.json"
    payload = json.loads(report_path.read_text(encoding="utf-8"))
    assert payload["metrics"]["pipeline_version"] == "defaults-override"
