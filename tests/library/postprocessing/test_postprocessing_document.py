"""Tests for document post-processing helpers.

Changelog
~~~~~~~~
- 2025-02-??: Added regression test for ``postprocess_export_file`` basename
  normalisation when dealing with temporary export artefacts.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from library.config import IoCfg
from library.postprocessing import document as document_postprocessing


def test_postprocess_export_file_strips_tmp_suffix(monkeypatch, tmp_path: Path) -> None:
    """Temporary exports lose the leading dot and ``.tmp`` suffix when saved."""

    source_path = tmp_path / ".output.documents_test.csv.tmp"
    source_path.write_text("document_chembl_id\nCHEMBL:TEST\n", encoding="utf-8")

    processed_frame = pd.DataFrame(
        {
            "document_chembl_id": ["CHEMBL:TEST"],
            "completed": ["true"],
        }
    )
    monkeypatch.setattr(
        document_postprocessing,
        "preprocess_document_export",
        lambda *_, **__: processed_frame,
    )

    cfg = IoCfg(csv_sep=",", csv_encoding="utf-8")

    destination = document_postprocessing.postprocess_export_file(
        source_path,
        cfg=cfg,
    )

    assert destination.name == "preprocessed_output.documents_test.csv"
    assert destination.exists()
