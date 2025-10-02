"""Tests for the QA checker comparing Power Query and Python outputs."""

from __future__ import annotations

import json
import shutil
from pathlib import Path

import pandas as pd
import pytest


from library import document_postprocessing as dp
from library.config import IoCfg

import qa.check_document_postprocessing as qa_module

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

def test_structure_metrics_detects_missing_candidate_columns(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    base_dir, reference_path, candidate_path = _prepare_environment(tmp_path)

    original_postprocess = dp.postprocess_file

    def drop_column_postprocess(
        source: Path,
        target_dir: Path,
        *,
        cfg: IoCfg,
        ref_document_path: Path,
    ) -> Path:
        processed_path = original_postprocess(
            source,
            target_dir,
            cfg=cfg,
            ref_document_path=ref_document_path,
        )
        processed_df = pd.read_csv(processed_path, dtype=str)
        if "doi" in processed_df.columns:
            processed_df = processed_df.drop(columns=["doi"])
        processed_df.to_csv(processed_path, index=False)
        return processed_path

    monkeypatch.setattr(dp, "postprocess_file", drop_column_postprocess)

    result = run_document_postprocessing_check(
        base_path=base_dir,
        reference_path=reference_path,
        candidate_path=candidate_path,
    )

    assert result.passed is False

    with result.report_json.open("r", encoding="utf-8") as handle:
        payload = json.load(handle)

    assert payload["structure"]["columns_equal"] is False
    assert any(
        issue.startswith("Candidate missing columns") and "doi" in issue
        for issue in payload["issues"]
    )

def test_diff_extract_limited_by_keys(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    base_dir, reference_path, candidate_path = _prepare_environment(tmp_path)

    rows: list[dict[str, str]] = []
    for idx in range(150):
        row = {column: "" for column in dp.FINAL_COLUMN_ORDER}
        row["PMID"] = f"{100000 + idx}"
        row["document_chembl_id"] = f"DOC{idx:05d}"
        row["completed"] = "2020-01-01"
        row["doi"] = f"10.1000/{idx:03d}"
        row["reference"] = "Journal, 1(1), p.1-2"
        row["sortorder.document"] = f"ISSN:{idx:08d}"
        row["review"] = "false"
        row["experimental"] = "true"
        row["document_contains_external_links"] = "false"
        row["invalid"] = "false"
        row["invalid.doi"] = "false"
        row["invalid.PMID"] = "false"
        row["invalid.reference"] = "false"
        row["year"] = "2020"
        row["month"] = "01"
        row["day"] = "01"
        rows.append(row)

    expected_df = pd.DataFrame(rows)
    actual_df = expected_df.copy()
    actual_df["doi"] = actual_df["doi"] + ".mismatch"

    expected_df.to_csv(candidate_path, index=False)

    monkeypatch.setattr(
        qa_module,
        "_power_query_expected",
        lambda ref_document, out_document: expected_df.copy(),
    )

    def fake_postprocess_file(
        source: Path,
        target_dir: Path,
        *,
        cfg: IoCfg,
        ref_document_path: Path,
    ) -> Path:
        processed_path = source.with_name(f"preprocessed_{source.name}")
        actual_df.to_csv(processed_path, index=False)
        return processed_path

    monkeypatch.setattr(dp, "postprocess_file", fake_postprocess_file)


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
            str(candidate_path.relative_to(base_dir)),

        ]
    )

    if mutate:
        assert exit_code == 1
    else:
        assert exit_code == 0
