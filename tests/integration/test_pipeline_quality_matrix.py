"""Regression tests covering the QA checklist for the CSV pipeline."""

from __future__ import annotations

import hashlib
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import pandas as pd
import pandas.testing as pdt
import pytest

REQUIRED_COLUMNS = [
    "molecule_chembl_id",
    "name",
    "smiles",
    "category",
    "score",
]
LOGGER_NAME = "tests.pipeline_quality"
EXPORT_COLUMNS = [
    "molecule_chembl_id",
    "name",
    "preferred_name",
    "smiles",
    "category_normalized",
    "score",
    "has_smiles",
    "is_high_score",
    "source",
]


@dataclass
class PipelineResult:
    output_path: Path
    frame: pd.DataFrame


def _normalise_booleans(frame: pd.DataFrame, columns: Iterable[str]) -> pd.DataFrame:
    normalised = frame.copy()
    for column in columns:
        normalised[column] = normalised[column].map({True: "true", False: "false"})
    return normalised


def run_pipeline(input_path: Path, dictionary_path: Path, destination: Path) -> PipelineResult:
    logger = logging.getLogger(LOGGER_NAME)
    logger.info(
        "pipeline_start input=%s dictionary=%s",
        str(input_path),
        str(dictionary_path),
    )

    frame = pd.read_csv(input_path, keep_default_na=False)
    missing = [column for column in REQUIRED_COLUMNS if column not in frame.columns]
    if missing:
        logger.error(
            "schema_missing_columns missing=%s",
            ",".join(sorted(str(column) for column in missing)),
        )
        raise ValueError(f"Input schema missing columns: {missing}")

    if frame.empty:
        logger.error("input_empty path=%s", str(input_path))
        raise ValueError("Input payload does not contain any rows")

    frame = frame.copy()
    frame["molecule_chembl_id"] = frame["molecule_chembl_id"].astype("string").str.strip()
    frame["name"] = frame["name"].astype("string").str.strip()
    frame["smiles"] = frame["smiles"].astype("string").str.strip()
    frame["category"] = frame["category"].astype("string").str.strip()
    frame["score"] = pd.to_numeric(frame["score"], errors="coerce").fillna(0.0)

    frame["smiles"] = frame["smiles"].replace({"N/A": "", "n/a": ""})
    frame["category"] = frame["category"].str.lower()

    duplicates = frame["molecule_chembl_id"].duplicated(keep=False)
    if duplicates.any():
        dup_ids = sorted(frame.loc[duplicates, "molecule_chembl_id"].unique().tolist())
        logger.warning(
            "duplicate_rows_dropped ids=%s count=%d",
            ",".join(dup_ids),
            int(duplicates.sum()),
        )
        frame = frame.sort_values(["molecule_chembl_id", "score"], ascending=[True, False])
        frame = frame.drop_duplicates("molecule_chembl_id", keep="first")

    dictionary = pd.read_csv(dictionary_path, keep_default_na=False)
    dictionary["molecule_chembl_id"] = dictionary["molecule_chembl_id"].astype("string").str.strip()
    dictionary["preferred_name"] = dictionary["preferred_name"].astype("string").str.strip()
    dictionary["category_override"] = (
        dictionary["category_override"].astype("string").str.strip().str.lower()
    )

    merged = frame.merge(dictionary, on="molecule_chembl_id", how="left")

    merged["preferred_name"] = merged["preferred_name"].astype("string").fillna("").str.strip()
    merged["category_override"] = (
        merged["category_override"].astype("string").fillna("").str.strip().str.lower()
    )

    missing_lookup = merged["preferred_name"].eq("")
    if missing_lookup.any():
        missing_ids = sorted(
            merged.loc[missing_lookup, "molecule_chembl_id"].unique().tolist()
        )
        logger.warning(
            "dictionary_missing_entries ids=%s count=%d",
            ",".join(missing_ids),
            len(missing_ids),
        )

    merged.loc[merged["name"].eq(""), "name"] = merged.loc[merged["name"].eq(""), "preferred_name"]
    merged.loc[merged["name"].eq(""), "name"] = "Undefined"
    merged.loc[merged["preferred_name"].eq(""), "preferred_name"] = pd.NA

    category_override = merged["category_override"].where(~merged["category_override"].eq(""))
    merged["category_normalized"] = category_override.fillna(merged["category"].where(~merged["category"].eq(""), "unknown"))

    merged["has_smiles"] = merged["smiles"].ne("")
    merged.loc[merged["smiles"].eq(""), "smiles"] = pd.NA
    merged["is_high_score"] = merged["score"] >= 1.0
    merged["source"] = "test-fixture"

    merged = merged[EXPORT_COLUMNS]
    merged = merged.sort_values("molecule_chembl_id").reset_index(drop=True)
    merged = _normalise_booleans(merged, ["has_smiles", "is_high_score"])

    destination.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(destination, index=False)
    logger.info(
        "pipeline_success output=%s rows=%d",
        str(destination),
        len(merged),
    )

    return PipelineResult(destination, merged)


@pytest.fixture()
def pipeline_resources(snapshot_resource: Path) -> Path:
    return snapshot_resource / "pipeline_quality"


@pytest.mark.integration
def test_pipeline_end_to_end__produces_golden_output(
    tmp_path: Path,
    pipeline_resources: Path,
    caplog: pytest.LogCaptureFixture,
) -> None:
    caplog.set_level(logging.INFO, LOGGER_NAME)
    input_path = pipeline_resources / "input_valid.csv"
    dictionary_path = pipeline_resources / "dictionary.csv"
    output_path = tmp_path / "final.csv"

    result = run_pipeline(input_path, dictionary_path, output_path)

    expected = pd.read_csv(
        pipeline_resources / "expected_output.csv",
        dtype={"has_smiles": "string", "is_high_score": "string"},
    )
    expected = expected.convert_dtypes(dtype_backend="numpy_nullable")
    actual = result.frame.convert_dtypes(dtype_backend="numpy_nullable")
    pdt.assert_frame_equal(actual, expected, check_like=True)
    assert output_path.exists()
    assert "pipeline_success" in caplog.text
    assert "duplicate_rows_dropped" in caplog.text
    assert "dictionary_missing_entries" in caplog.text


@pytest.mark.integration
def test_pipeline_schema_validation__missing_column(
    tmp_path: Path,
    pipeline_resources: Path,
    caplog: pytest.LogCaptureFixture,
) -> None:
    caplog.set_level(logging.INFO, LOGGER_NAME)
    input_path = pipeline_resources / "input_missing_column.csv"
    dictionary_path = pipeline_resources / "dictionary.csv"
    output_path = tmp_path / "final.csv"

    with pytest.raises(ValueError) as excinfo:
        run_pipeline(input_path, dictionary_path, output_path)

    assert "missing columns" in str(excinfo.value)
    assert "schema_missing_columns" in caplog.text
    assert not output_path.exists()


@pytest.mark.integration
def test_pipeline_degradation_cases__empty_file(
    tmp_path: Path,
    pipeline_resources: Path,
    caplog: pytest.LogCaptureFixture,
) -> None:
    caplog.set_level(logging.INFO, LOGGER_NAME)
    input_path = pipeline_resources / "input_empty.csv"
    dictionary_path = pipeline_resources / "dictionary.csv"
    output_path = tmp_path / "final.csv"

    with pytest.raises(ValueError) as excinfo:
        run_pipeline(input_path, dictionary_path, output_path)

    assert "does not contain any rows" in str(excinfo.value)
    assert "input_empty" in caplog.text
    assert not output_path.exists()


@pytest.mark.integration
def test_pipeline_idempotence__repeated_runs_are_identical(
    tmp_path: Path,
    pipeline_resources: Path,
) -> None:
    input_path = pipeline_resources / "input_valid.csv"
    dictionary_path = pipeline_resources / "dictionary.csv"

    first_output = tmp_path / "run-1" / "final.csv"
    second_output = tmp_path / "run-2" / "final.csv"

    first = run_pipeline(input_path, dictionary_path, first_output)
    second = run_pipeline(input_path, dictionary_path, second_output)

    assert first.frame.equals(second.frame)

    def _digest(path: Path) -> str:
        data = path.read_bytes()
        return hashlib.sha256(data).hexdigest()

    assert _digest(first_output) == _digest(second_output)
