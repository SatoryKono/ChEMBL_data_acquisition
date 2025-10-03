from pathlib import Path

import pandas as pd

from library.config import IoCfg
from library.postprocessing import document as document_export_postprocessing


FIXTURE_DIR = Path("tests/data/postprocessing_document")


def test_postprocess_export_file_strips_tmp_suffix(tmp_path: Path) -> None:
    """Temporary suffixes are removed when deriving the output path."""

    source = FIXTURE_DIR / "output.document_20230101.csv"
    input_path = tmp_path / ".output.documents_test.csv.tmp"
    input_path.write_bytes(source.read_bytes())

    cfg = IoCfg(csv_sep=",", csv_encoding="utf-8")
    output_path = document_export_postprocessing.postprocess_export_file(
        input_path,
        cfg=cfg,
        ref_document_path=FIXTURE_DIR / "document.csv",
    )

    assert output_path.name == "preprocessed_output.documents_test.csv"
    assert output_path.exists()

    produced = pd.read_csv(output_path, dtype=str)
    assert not produced.empty
