"""Tests for the QA checker comparing Power Query and Python outputs."""

from __future__ import annotations

import json
import shutil
from pathlib import Path

import pandas as pd
import pytest

from qa.check_document_postprocessing import MAX_DIFF_KEY_EXPORT
from qa.check_document_postprocessing import main as qa_main
from qa.check_document_postprocessing import run_document_postprocessing_check


FIXTURE_DIR = Path("tests/data/postprocessing_document")


def _prepare_environment(tmp_path: Path) -> tuple[Path, Path, Path]:
    base_dir = tmp_path / "data"
    reference_dir = base_dir / "input" / "full"
    output_dir = base_dir / "output"
    reference_dir.mkdir(parents=True)
    output_dir.mkdir(parents=True)

    reference_path = reference_dir / "document.csv"
    candidate_path = output_dir / "output.document_20230101.csv"

    shutil.copy(FIXTURE_DIR / "ref_document.csv", reference_path)
    shutil.copy(
        FIXTURE_DIR / "preprocessed_output.document_20230101.csv",
        candidate_path,
    )

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

    with result.report_json.open("r", encoding="utf-8") as handle:
        payload = json.load(handle)
    assert payload["status"] == "passed"
    assert payload["structure"]["columns_equal"] is True
    assert payload["structure"]["column_order_equal"] is True
    assert payload["reference"]["review"]["share_true"] == pytest.approx(2 / 3)
    assert payload["candidate"]["invariants"]["invalid_formula"]["passed"] is True

    markdown = result.report_markdown.read_text(encoding="utf-8")
    assert "- Status: **passed**" in markdown
    assert "- Diff excerpt: not generated" in markdown


def test_run_document_postprocessing_check_failure(tmp_path: Path) -> None:
    base_dir, reference_path, candidate_path = _prepare_environment(tmp_path)

    mutated = pd.read_csv(candidate_path, dtype=str)
    mutated.loc[0, "doi"] = "10.9999/incorrect"
    mutated.to_csv(candidate_path, index=False)

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
    assert result.diff_csv.name in markdown


def test_diff_extract_limited_by_keys(tmp_path: Path) -> None:
    base_dir, reference_path, candidate_path = _prepare_environment(tmp_path)

    template = pd.read_csv(reference_path, dtype=str)
    rows: list[pd.Series] = []
    for idx in range(150):
        record = template.iloc[0].copy()
        record["PMID"] = f"{100000 + idx}"
        record["document_chembl_id"] = f"DOC{idx:05d}"
        record["completed"] = "2020-01-01"
        record["doi"] = f"10.1000/{idx:03d}"
        rows.append(record)
    reference_df = pd.DataFrame(rows)
    reference_df.to_csv(reference_path, index=False)

    candidate_df = reference_df.copy()
    candidate_df["doi"] = candidate_df["doi"] + ".mismatch"
    candidate_df.to_csv(candidate_path, index=False)

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
def test_cli_exit_codes(tmp_path: Path, mutate: bool) -> None:
    base_dir, reference_path, candidate_path = _prepare_environment(tmp_path)
    if mutate:
        mutated = pd.read_csv(candidate_path, dtype=str)
        mutated.loc[0, "doi"] = "bad-doi"
        mutated.to_csv(candidate_path, index=False)

    exit_code = qa_main(
        [
            "--base-path",
            str(base_dir),
            "--out",
            "output\\output.document_20230101.csv",
        ]
    )

    if mutate:
        assert exit_code == 1
    else:
        assert exit_code == 0
