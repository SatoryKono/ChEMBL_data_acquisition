from __future__ import annotations

import shutil
from pathlib import Path

import pandas as pd
import pytest

from library.pipelines.document import postprocessing as stage_dp
from library.postprocessing import document as doc


FIXTURE_DIR = Path("tests/data/postprocessing_document")
BOOL_COLUMNS = {
    "review",
    "experimental",
    "document_contains_external_links",
    "invalid",
    "invalid.doi",
    "invalid.PMID",
    "invalid.reference",
}


def _normalise_strings(df: pd.DataFrame) -> pd.DataFrame:
    normalised = df.apply(
        lambda col: col.map(lambda value: "" if pd.isna(value) else str(value))
    )
    for column in BOOL_COLUMNS.intersection(normalised.columns):
        normalised[column] = normalised[column].str.lower()
    return normalised


def test_normalize_journal_and_padding_functions() -> None:
    assert doc.normalize_journal(" J. Clin. Chem. ") == "j clin chem"
    assert doc.normalize_journal(None) == ""

    assert doc.pad4(5) == "0005"
    assert doc.pad4(None) == "0000"
    assert doc.pad2("9") == "09"
    assert doc.pad2(0) == "00"
    assert doc.pad_pmid8("123") == "00000123"
    assert doc.pad_pmid8(0) == ""


def test_eq_ne_text() -> None:
    assert doc.eq_text("10.1000/test", "10.1000/test")
    assert not doc.eq_text("", "10.1000/test")
    assert doc.ne_text("10.2000/test", "10.1000/test")
    assert not doc.ne_text("", "10.1000/test")


def test_invalid_flag_detection() -> None:
    ref_frame = doc._load_reference_document(FIXTURE_DIR / "document.csv")
    out_frame = doc._load_output_document(FIXTURE_DIR / "output.document_20230101.csv")

    harmonised = doc._harmonise_documents(out_frame, ref_frame)
    indexed = harmonised.set_index("document_chembl_id")

    assert bool(indexed.loc["DOC2", "invalid.doi"]) is True
    assert bool(indexed.loc["DOC2", "invalid.reference"]) is True
    assert bool(indexed.loc["DOC2", "invalid.PMID"]) is False

    assert bool(indexed.loc["DOC3", "invalid.PMID"]) is True
    assert bool(indexed.loc["DOC3", "invalid.doi"]) is False
    assert bool(indexed.loc["DOC3", "invalid.reference"]) is False

    aggregated = indexed["invalid"]
    recomputed = (
        indexed["invalid.doi"]
        | indexed["invalid.PMID"]
        | indexed["invalid.reference"]
    )
    pd.testing.assert_series_equal(
        aggregated.sort_index(),
        recomputed.sort_index(),
        check_names=False,
    )

    assert "unknown" in indexed.loc["DOC2", "reference"]
    assert indexed.loc["DOC1", "reference"].startswith("journal of testing")
    assert bool(indexed.loc["DOC1", "review"]) is True
    assert bool(indexed.loc["DOC1", "experimental"]) is False
    assert bool(indexed.loc["DOC3", "review"]) is False
    assert bool(indexed.loc["DOC3", "experimental"]) is True

    mesh_lists = indexed.loc[:, [
        "PubMed.mesh.descriptors",
        "PubMed.mesh.qualifiers",
        "PubMed.chemical.list",
        "OpenAlex.mesh.descriptors",
    ]]
    for column in mesh_lists.columns:
        assert all(value == value.lower() for value in mesh_lists[column].fillna(""))


def test_preprocess_documents_csv_integration(tmp_path: Path) -> None:
    base_dir = tmp_path / "data"
    (base_dir / "input" / "full").mkdir(parents=True)
    (base_dir / "output" / "document").mkdir(parents=True)

    shutil.copy(FIXTURE_DIR / "document.csv", base_dir / "input" / "full" / "document.csv")
    shutil.copy(
        FIXTURE_DIR / "output.document_20230101.csv",
        base_dir / "output" / "document" / "output.document_20230101.csv",
    )
    shutil.copy(
        FIXTURE_DIR / "ref_document.csv",
        base_dir / "input" / "full" / "ref_document.csv",
    )

    result_path = doc.preprocess_documents_csv(
        base_path=str(base_dir),
        ref_document_rel="input\\full\\document.csv",
        out_document_rel="output\\document\\output.document_20230101.csv",
        qa_reference_rel=None,
    )

    result = Path(result_path)
    expected = FIXTURE_DIR / "preprocessed_output.document_20230101.csv"

    produced_df = pd.read_csv(result, dtype=str)
    expected_df = pd.read_csv(expected, dtype=str)

    assert result.name == "preprocessed_output.document_20230101.csv"
    assert list(produced_df.columns) == list(doc.FINAL_COLUMN_ORDER)
    pd.testing.assert_frame_equal(
        _normalise_strings(produced_df),
        _normalise_strings(expected_df),
        check_dtype=False,
    )

    completed = produced_df["completed"].tolist()
    assert completed == sorted(completed)

    qa_json = result.parent / "qa_document_postprocessing_report_20230101.json"
    qa_md = result.parent / "qa_document_postprocessing_report_20230101.md"
    qa_diff = result.parent / "qa_document_postprocessing_diff_20230101.csv"

    assert not qa_json.exists()
    assert not qa_md.exists()
    assert not qa_diff.exists()


def test_stage_postprocess_documents_missing_reference(tmp_path: Path) -> None:
    out_frame = pd.read_csv(
        FIXTURE_DIR / "output.document_20230101.csv",
        dtype=str,
    )
    missing_path = tmp_path / "missing.csv"

    with pytest.raises(FileNotFoundError, match="Reference document CSV not found"):
        stage_dp.postprocess_documents(out_frame, ref_document_path=missing_path)


def test_preprocess_document_export_missing_reference(tmp_path: Path) -> None:
    out_frame = pd.read_csv(
        FIXTURE_DIR / "output.document_20230101.csv",
        dtype=str,
    )
    missing_path = tmp_path / "missing.csv"

    with pytest.raises(FileNotFoundError, match="Reference document CSV not found"):
        doc.preprocess_document_export(out_frame, ref_document_path=missing_path)
