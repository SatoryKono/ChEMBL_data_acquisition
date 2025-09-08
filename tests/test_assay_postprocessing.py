"""Tests for :mod:`library.assay_postprocessing`."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from library import assay_postprocessing as ap


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


def test_postprocess_file_roundtrip(tmp_path: Path) -> None:
    """``postprocess_file`` reads, processes and writes data."""

    df = pd.DataFrame(
        {
            "document_chembl_id": ["doc1", "doc1"],
            "target_chembl_id": ["t1", "t1"],
            "assay_chembl_id": ["a1", "a2"],
        }
    )
    input_path = tmp_path / "in.csv"
    df.to_csv(input_path, index=False)
    output_path = tmp_path / "out.csv"

    ap.postprocess_file(input_path, output_path)

    result = pd.read_csv(output_path)
    assert list(result["assay_with_same_target"]) == [2, 2]
