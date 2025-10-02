"""Tests for the QA checker comparing Power Query and Python outputs."""

from __future__ import annotations

import json
import shutil
from pathlib import Path

import pandas as pd
import pytest

from library import document_postprocessing as dp
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
            "--ref",
            "input\\full\\document.csv",
            "--actual",
            "output\\document\\output.document_20230101.csv",
        ]
    )

    if mutate:
        assert exit_code == 1
    else:
        assert exit_code == 0
