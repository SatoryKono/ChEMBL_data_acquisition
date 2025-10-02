"""Tests for :mod:`library.document_postprocessing`."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from library import document_postprocessing as dp
from library.config import IoCfg


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


def test_helper_functions_match_m_script_behaviour() -> None:
    """Utility helpers mimic the Power Query guards."""

    assert dp.null_or_empty("")
    assert dp.null_or_empty(0)
    assert not dp.null_or_empty("1")

    assert dp.normalize_journal(" Journal. Of Testing. ") == "journal of testing"
    assert dp.pad4("5") == "0005"
    assert dp.pad2(9) == "09"
    assert dp.pad_pmid8(123) == "00000123"

    assert dp.eq_text("abc", "abc") is True
    assert dp.eq_text("", "abc") is False
    assert dp.ne_text("abc", "def") is True
    assert dp.ne_text("", "def") is False


def test_postprocess_documents_matches_stage_pipeline() -> None:
    """The Python port produces the same result as the Stage pipeline."""

    ref_frame = _frame_from_csv(FIXTURE_DIR / "document.csv")
    out_frame = _frame_from_csv(FIXTURE_DIR / "output.document_20230101.csv")
    expected = _frame_from_csv(
        FIXTURE_DIR / "preprocessed_output.document_20230101.csv"
    )

    result = dp.postprocess_documents(out_frame, ref_document=ref_frame)

    assert list(result.columns) == list(dp.FINAL_COLUMN_ORDER)
    pd.testing.assert_frame_equal(
        _normalise_strings(result),
        _normalise_strings(expected),
        check_dtype=False,
    )

    indexed = result.set_index("document_chembl_id")
    assert bool(indexed.loc["DOC1", "review"]) is True
    assert bool(indexed.loc["DOC1", "experimental"]) is False
    assert bool(indexed.loc["DOC2", "invalid.doi"]) is True
    assert bool(indexed.loc["DOC2", "invalid.reference"]) is True
    assert bool(indexed.loc["DOC3", "document_contains_external_links"]) is True


def test_postprocess_file_creates_prefixed_output(tmp_path: Path) -> None:
    """``postprocess_file`` writes UTF-8 CSV files using the expected prefix."""

    input_path = FIXTURE_DIR / "output.document_20230101.csv"
    ref_path = FIXTURE_DIR / "document.csv"
    cfg = IoCfg(csv_sep=",", csv_encoding="utf-8")

    destination = dp.postprocess_file(
        input_path,
        tmp_path,
        cfg=cfg,
        ref_document_path=ref_path,
    )

    assert destination.name.startswith(dp.OUTPUT_PREFIX)

    produced = pd.read_csv(destination, dtype=str)
    expected = _frame_from_csv(
        FIXTURE_DIR / "preprocessed_output.document_20230101.csv"
    )
    pd.testing.assert_frame_equal(
        _normalise_strings(produced),
        _normalise_strings(expected),
        check_dtype=False,
    )

