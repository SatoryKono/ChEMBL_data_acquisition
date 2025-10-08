from __future__ import annotations

import importlib
import logging
import types
from pathlib import Path

import pandas as pd
import pytest

from library.postprocessing import document as document_module


def _install_stub_qa(monkeypatch: pytest.MonkeyPatch, metrics: dict[str, object]) -> None:
    """Patch :mod:`importlib` to return a stub QA module with *metrics*."""

    class _StubCrosswalk:
        @staticmethod
        def load(path: Path) -> object:  # pragma: no cover - behaviour is trivial
            return object()

    stub_module = types.SimpleNamespace(
        Crosswalk=_StubCrosswalk,
        CROSSWALK_PATH=str(Path("stub_crosswalk.yaml")),
        DEFAULT_DIFF_LIMIT=100,
        DEFAULT_REPORT_DIR=Path("output") / "document",
        run_document_postprocessing_check=lambda **kwargs: metrics,
    )

    original_import = importlib.import_module

    def _fake_import(name: str, package: str | None = None):
        if name in {
            "qa.check_document_postprocessing",
            "library.qa.check_document_postprocessing",
        }:
            return stub_module
        return original_import(name, package)

    monkeypatch.setattr(importlib, "import_module", _fake_import)


def _setup_frames() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    reference = pd.DataFrame(
        {
            "document_chembl_id": ["DOC1", "DOC2"],
            "primary_pubmed_id": ["PM1", "PM2"],
        }
    )
    output = pd.DataFrame(
        {
            "document_chembl_id": ["DOC1"],
            "primary_pubmed_id": ["PM1"],
        }
    )
    harmonised = pd.DataFrame(
        {
            "document_chembl_id": ["DOC1"],
            "primary_pubmed_id": ["PM1"],
            "invalid": [False],
            "invalid.doi": [False],
            "invalid.PMID": [False],
            "invalid.reference": [False],
        }
    )
    return reference, output, harmonised


@pytest.mark.unit
def test_preprocess_documents_csv__qa_subset_missing_rows(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, caplog: pytest.LogCaptureFixture
) -> None:
    data_dir = tmp_path / "data"
    (data_dir / "input" / "full").mkdir(parents=True)
    (data_dir / "output" / "document").mkdir(parents=True)
    (data_dir / "input" / "full" / "document.csv").write_text("reference", encoding="utf-8")
    (data_dir / "output" / "document" / "output.document_20250101.csv").write_text(
        "output",
        encoding="utf-8",
    )

    reference, output, harmonised = _setup_frames()

    monkeypatch.setattr(document_module, "_load_reference_document", lambda path: reference)
    monkeypatch.setattr(document_module, "_load_output_document", lambda path: output)
    monkeypatch.setattr(document_module, "_harmonise_documents", lambda *_: harmonised)

    metrics = {
        "status": "FAIL",
        "report_json": str(data_dir / "output" / "document" / "report.json"),
        "diff_path": None,
        "issues": ["1 M-output-only keys detected"],
        "missing_rows": {"python_only": 0, "m_only": 1},
        "diff_rows": 0,
    }
    _install_stub_qa(monkeypatch, metrics)

    caplog.set_level(logging.INFO)

    result = document_module.preprocess_documents_csv(
        base_path=str(data_dir),
        ref_document_rel="input\\full\\document.csv",
        out_document_rel="output\\document\\output.document_20250101.csv",
    )

    expected_path = data_dir / "output" / "document" / "preprocessed_output.document_20250101.csv"
    assert Path(result) == expected_path
    assert expected_path.exists()

    assert any(
        "document_postprocess_qa_skipped_subset" in record.message for record in caplog.records
    )


@pytest.mark.unit
def test_preprocess_documents_csv__qa_failure_propagates(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    data_dir = tmp_path / "data"
    (data_dir / "input" / "full").mkdir(parents=True)
    (data_dir / "output" / "document").mkdir(parents=True)
    (data_dir / "input" / "full" / "document.csv").write_text("reference", encoding="utf-8")
    (data_dir / "output" / "document" / "output.document_20250101.csv").write_text(
        "output",
        encoding="utf-8",
    )

    reference, output, harmonised = _setup_frames()

    monkeypatch.setattr(document_module, "_load_reference_document", lambda path: reference)
    monkeypatch.setattr(document_module, "_load_output_document", lambda path: output)
    monkeypatch.setattr(document_module, "_harmonise_documents", lambda *_: harmonised)

    metrics = {
        "status": "FAIL",
        "report_json": str(data_dir / "output" / "document" / "report.json"),
        "diff_path": str(data_dir / "output" / "document" / "diff.csv"),
        "issues": ["Match rate for preferred below threshold (0.0000% < 99.50%)"],
        "missing_rows": {"python_only": 0, "m_only": 0},
        "diff_rows": 1,
    }
    _install_stub_qa(monkeypatch, metrics)

    with pytest.raises(RuntimeError):
        document_module.preprocess_documents_csv(
            base_path=str(data_dir),
            ref_document_rel="input\\full\\document.csv",
            out_document_rel="output\\document\\output.document_20250101.csv",
        )


@pytest.mark.unit
def test_load_output_document__preserves_page_ranges(tmp_path: Path) -> None:
    csv_path = tmp_path / "output.document_20251010.csv"

    row: dict[str, object] = {}
    for column, dtype in document_module.OUTPUT_DTYPE.items():
        row[column] = 0 if dtype == "Int64" else ""

    row["PubMed.StartPage"] = "11-12"
    row["PubMed.EndPage"] = "15"
    row["PubMed.Volume"] = "IV"
    row["PubMed.Issue"] = "S1"
    row["ChEMBL.volume"] = "VII"
    row["ChEMBL.issue"] = "Supplement"
    row["ChEMBL.first_page"] = "11-12"
    row["ChEMBL.last_page"] = "15A"

    pd.DataFrame([row]).to_csv(csv_path, index=False)

    frame = document_module._load_output_document(csv_path)

    assert frame.loc[0, "ChEMBL.first_page"] == "11-12"
    assert frame.loc[0, "ChEMBL.last_page"] == "15A"
    assert frame.loc[0, "ChEMBL.volume"] == "VII"
    assert frame.loc[0, "ChEMBL.issue"] == "Supplement"
    assert frame.loc[0, "PubMed.StartPage"] == "11-12"
    assert frame.loc[0, "PubMed.Volume"] == "IV"

    assert frame["ChEMBL.first_page"].dtype == pd.StringDtype()
    assert frame["PubMed.StartPage"].dtype == pd.StringDtype()
