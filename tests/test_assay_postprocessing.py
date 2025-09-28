"""Tests for :mod:`library.assay_postprocessing`."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pandera.errors as pa_errors
import pytest

from library import assay_postprocessing as ap
from library.config import IoCfg
from library.constants import AssayPostprocessSchema


def test_postprocess_assays_counts() -> None:
    """``postprocess_assays`` adds per-target counts."""

    df = pd.DataFrame(
        {
            "document_chembl_id": ["doc1", "doc1", "doc1", "doc2"],
            "target_chembl_id": ["t1", "t1", "t2", "t1"],
            "assay_chembl_id": ["a1", "a2", "a3", "a4"],
        }
    )
    out = ap.postprocess_assays(df)

    assert out.loc[0, "assay_with_same_target"] == 2
    assert out.loc[2, "assay_with_same_target"] == 1
    assert out.loc[3, "assay_with_same_target"] == 1


def test_postprocess_assays_invalid_schema() -> None:
    """Missing required columns trigger schema errors."""

    df = pd.DataFrame({"document_chembl_id": ["doc1"], "assay_chembl_id": ["a1"]})
    with pytest.raises(pa_errors.SchemaError):
        ap.postprocess_assays(df)


def test_postprocess_file_roundtrip(tmp_path: Path) -> None:
    """``postprocess_file`` honours ``IoCfg`` defaults for CSV handling."""

    df = pd.DataFrame(
        {
            "document_chembl_id": ["doc1", "doc1"],
            "target_chembl_id": ["t1", "t1"],
            "assay_chembl_id": ["a1", "a2"],
        }
    )
    cfg = IoCfg(csv_sep=";", csv_encoding="utf8")
    input_path = tmp_path / "in.csv"
    df.to_csv(input_path, index=False, sep=cfg.csv_sep, encoding=cfg.csv_encoding)
    output_path = tmp_path / "out.csv"

    ap.postprocess_file(input_path, output_path, cfg=cfg)

    result = pd.read_csv(output_path, sep=cfg.csv_sep, encoding=cfg.csv_encoding)
    assert list(result["assay_with_same_target"]) == [2, 2]


def test_assay_postprocess_schema() -> None:
    """Direct validation against the schema works as expected."""

    df_valid = pd.DataFrame(
        {"document_chembl_id": ["doc1"], "target_chembl_id": ["t1"]}
    )
    AssayPostprocessSchema.validate(df_valid)

    df_invalid = pd.DataFrame({"document_chembl_id": ["doc1"]})
    with pytest.raises(pa_errors.SchemaError):
        AssayPostprocessSchema.validate(df_invalid)
