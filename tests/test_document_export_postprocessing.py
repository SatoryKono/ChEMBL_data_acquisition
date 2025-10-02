from pathlib import Path

import pandas as pd

from library.pipelines.document import postprocessing as dp
from library.config import IoCfg
from library.postprocessing import document as document_export_postprocessing


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


def _frame_from_csv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, dtype=str)


def _normalise_strings(df: pd.DataFrame) -> pd.DataFrame:
    normalised = df.apply(
        lambda col: col.map(lambda value: "" if pd.isna(value) else str(value))
    )
    for column in BOOL_COLUMNS.intersection(normalised.columns):
        normalised[column] = normalised[column].str.lower()
    return normalised


def test_preprocess_document_export_matches_stage_projection() -> None:
    """The export wrapper mirrors the Stage post-processing output."""

    input_frame = _frame_from_csv(FIXTURE_DIR / "output.document_20230101.csv")
    reference_frame = _frame_from_csv(FIXTURE_DIR / "document.csv")
    expected_frame = _frame_from_csv(
        FIXTURE_DIR / "preprocessed_output.document_20230101.csv"
    )

    result = document_export_postprocessing.preprocess_document_export(
        input_frame,
        ref_document=reference_frame,
    )

    assert list(result.columns) == list(dp.FINAL_COLUMN_ORDER)
    assert "preferred_title" not in result.columns
    assert "metadata_sources" not in result.columns

    pd.testing.assert_frame_equal(
        _normalise_strings(result),
        _normalise_strings(expected_frame),
        check_dtype=False,
    )


def test_postprocess_export_file_writes_stage_projection(tmp_path: Path) -> None:
    """The CSV artefact follows the Stage column ordering."""

    raw_frame = _frame_from_csv(FIXTURE_DIR / "output.document_20230101.csv")
    input_path = tmp_path / "output.document_20240101.csv"
    raw_frame.to_csv(input_path, index=False)

    cfg = IoCfg(csv_sep=",", csv_encoding="utf-8")
    output_path = document_export_postprocessing.postprocess_export_file(
        input_path,
        cfg=cfg,
        ref_document_path=FIXTURE_DIR / "document.csv",
    )

    assert output_path.name == f"preprocessed_{input_path.name}"

    produced = pd.read_csv(output_path, dtype=str)
    expected = _frame_from_csv(FIXTURE_DIR / "preprocessed_output.document_20230101.csv")

    assert list(produced.columns) == list(dp.FINAL_COLUMN_ORDER)
    pd.testing.assert_frame_equal(
        _normalise_strings(produced),
        _normalise_strings(expected),
        check_dtype=False,
    )
