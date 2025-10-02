"""Tests for the QA checker comparing Power Query and Python outputs."""

from __future__ import annotations

import json
import shutil
from pathlib import Path

import pandas as pd
import pytest


from library import document_postprocessing as dp

from qa.check_document_postprocessing import MAX_DIFF_KEY_EXPORT

from qa.check_document_postprocessing import main as qa_main
from qa.check_document_postprocessing import run_document_postprocessing_check


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


def _prepare_environment(tmp_path: Path) -> tuple[Path, Path, Path]:
    base_dir = tmp_path / "data"
    reference_dir = base_dir / "input" / "full"
    output_dir = base_dir / "output" / "document"
    reference_dir.mkdir(parents=True)
    output_dir.mkdir(parents=True)

    reference_path = reference_dir / "document.csv"

    legacy_reference_path = reference_dir / "ref_document.csv"
    candidate_path = output_dir / "output.document_20230101.csv"

    shutil.copy(FIXTURE_DIR / "document.csv", reference_path)
    shutil.copy(FIXTURE_DIR / "document.csv", legacy_reference_path)
    shutil.copy(FIXTURE_DIR / "output.document_20230101.csv", candidate_path)


    return base_dir, reference_path, candidate_path


def test_run_document_postprocessing_check_pass(tmp_path: Path) -> None:
    base_dir, reference_path, candidate_path = _prepare_environment(tmp_path)

    result = run_document_postprocessing_check(
        base_path=base_dir,
        reference_path=reference_path,
        candidate_path=candidate_path,
    )

    assert result.passed is True
    assert result.diff_csv is None

    produced_path = candidate_path.with_name(
        f"preprocessed_{candidate_path.name}"
    )
    assert produced_path.exists()

    produced_df = pd.read_csv(produced_path, dtype=str)
    expected_df = pd.read_csv(
        FIXTURE_DIR / "preprocessed_output.document_20230101.csv", dtype=str
    )
    for column in BOOL_COLUMNS:
        if column in produced_df.columns:
            produced_df[column] = produced_df[column].str.lower()
    for column in BOOL_COLUMNS:
        if column in expected_df.columns:
            expected_df[column] = expected_df[column].str.lower()

    pd.testing.assert_frame_equal(
        produced_df.fillna(""), expected_df.fillna(""), check_dtype=False
    )

    with result.report_json.open("r", encoding="utf-8") as handle:
        payload = json.load(handle)
    assert payload["status"] == "passed"
    assert payload["structure"]["columns_equal"] is True
    assert payload["structure"]["column_order_equal"] is True
    assert payload["reference"]["review"]["share_true"] == pytest.approx(2 / 3)
    assert payload["candidate"]["invariants"]["invalid_formula"]["passed"] is True

    markdown = result.report_markdown.read_text(encoding="utf-8")
    assert "- Status: **passed**" in markdown
    assert "- Column sets identical: ✅ yes" in markdown
    assert "- Column order identical: ✅ yes" in markdown
    assert "- Diff excerpt: not generated" in markdown


def test_run_document_postprocessing_check_failure(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    base_dir, reference_path, candidate_path = _prepare_environment(tmp_path)

    original = dp.postprocess_documents

    def mutate(
        frame, *, required_columns=None, ref_document=None, ref_document_path=None
    ):
        processed = original(
            frame,
            required_columns=required_columns,
            ref_document=ref_document,
            ref_document_path=ref_document_path,
        )
        processed.loc[0, "doi"] = "10.9999/incorrect"
        return processed

    monkeypatch.setattr(dp, "postprocess_documents", mutate)

    result = run_document_postprocessing_check(
        base_path=base_dir,
        reference_path=reference_path,
        candidate_path=candidate_path,
        max_diff_rows=5,
    )

    assert result.passed is False
    assert result.diff_csv is not None

    diff_df = pd.read_csv(result.diff_csv, dtype=str)
    assert len(diff_df) <= 5
    assert "doi" in diff_df["column"].unique()

    with result.report_json.open("r", encoding="utf-8") as handle:
        payload = json.load(handle)
    assert payload["status"] == "failed"
    assert payload["differences"]["rows_with_differences"] >= 1

    markdown = result.report_markdown.read_text(encoding="utf-8")
    assert "- Status: **failed**" in markdown
    assert "- Column sets identical:" in markdown
    assert "- Column order identical:" in markdown
    assert result.diff_csv.name in markdown


def test_diff_extract_limited_by_keys(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    base_dir, reference_path, candidate_path = _prepare_environment(tmp_path)

    candidate_template = pd.read_csv(candidate_path, dtype=str)
    candidate_rows: list[pd.Series] = []
    for idx in range(150):
        record = candidate_template.iloc[0].copy()
        record["ChEMBL.document_chembl_id"] = f"DOC{idx:05d}"
        record["ChEMBL.doi"] = f"10.1000/{idx:03d}"
        record["ChEMBL.pubmed_id"] = f"{100000 + idx}"
        record["PubMed.PMID"] = f"{100000 + idx}"
        record["PubMed.YearCompleted"] = "2020"
        record["PubMed.MonthCompleted"] = "01"
        record["PubMed.DayCompleted"] = "01"
        record["ChEMBL.year"] = "2020"
        record["ChEMBL.volume"] = "1"
        record["ChEMBL.issue"] = "1"
        record["ChEMBL.first_page"] = "1"
        record["ChEMBL.last_page"] = "2"
        candidate_rows.append(record)
    candidate_raw = pd.DataFrame(candidate_rows)
    candidate_raw.to_csv(candidate_path, index=False)

    original = dp.postprocess_documents

    def mutate(
        frame,
        *,
        required_columns=None,
        ref_document=None,
        ref_document_path=None,
    ):
        processed = original(
            frame,
            required_columns=required_columns,
            ref_document=ref_document,
            ref_document_path=ref_document_path,
        )
        processed["doi"] = processed["doi"].astype(str) + ".mismatch"
        return processed

    monkeypatch.setattr(dp, "postprocess_documents", mutate)

    result = run_document_postprocessing_check(
        base_path=base_dir,
        reference_path=reference_path,
        candidate_path=candidate_path,
    )

    assert result.passed is False
    assert result.diff_csv is not None

    diff_df = pd.read_csv(result.diff_csv, dtype=str)
    assert diff_df["PMID"].nunique() == MAX_DIFF_KEY_EXPORT
    assert set(diff_df["column"]) == {"doi"}


@pytest.mark.parametrize(
    "mutate",
    [False, True],
)
def test_cli_exit_codes(
    tmp_path: Path, mutate: bool, monkeypatch: pytest.MonkeyPatch
) -> None:
    base_dir, reference_path, candidate_path = _prepare_environment(tmp_path)
    if mutate:
        original = dp.postprocess_documents

        def mutate_fn(
            frame,
            *,
            required_columns=None,
            ref_document=None,
            ref_document_path=None,
        ):
            processed = original(
                frame,
                required_columns=required_columns,
                ref_document=ref_document,
                ref_document_path=ref_document_path,
            )
            processed.loc[0, "doi"] = "bad-doi"
            return processed

        monkeypatch.setattr(dp, "postprocess_documents", mutate_fn)

    exit_code = qa_main(
        [
            "--base-path",
            str(base_dir),
            "--out",
            "output\\document\\output.document_20230101.csv",
        ]
    )

    if mutate:
        assert exit_code == 1
    else:
        assert exit_code == 0
