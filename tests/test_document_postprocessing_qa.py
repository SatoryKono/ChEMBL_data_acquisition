"""Tests for the QA checker comparing Power Query and Python outputs."""

from __future__ import annotations

import json
import shutil
from pathlib import Path

import pandas as pd
import pytest

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
