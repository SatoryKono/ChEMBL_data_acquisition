from pathlib import Path

import pandas as pd
import pytest

from library.config import Config
from library.pipelines.cellline import (
    CELL_LINE_BASE_COLUMNS,
    prepare_cellline_dataframe,
    read_cellline_identifiers,
)


@pytest.fixture()
def cellline_input_csv(tmp_path: Path) -> Path:
    path = tmp_path / "cellline.csv"
    path.write_text(
        "cell_chembl_id\nCHEMBL3307636\n CHEMBL3307790 \n#N/A\n",
        encoding="utf-8",
    )
    return path


def test_read_cellline_identifiers__applies_limit_and_offset(
    cfg: Config, cellline_input_csv: Path
) -> None:
    identifiers = read_cellline_identifiers(
        cellline_input_csv,
        column="cell_chembl_id",
        io_cfg=cfg.io,
        limit=1,
        offset=1,
    )

    assert identifiers == ["CHEMBL3307790"]


def test_prepare_cellline_dataframe__adds_missing_identifiers() -> None:
    raw = pd.DataFrame(
        {
            "cell_chembl_id": ["CHEMBL3307636"],
            "cell_name": ["22Rv1"],
            "cell_id": [1234],
            "cell_source_tax_id": [9606],
        }
    )
    prepared, missing = prepare_cellline_dataframe(
        raw,
        ["CHEMBL3307636", "CHEMBL3307790"],
    )

    assert tuple(missing) == ("CHEMBL3307790",)
    assert set(prepared["cell_chembl_id"].dropna()) == {
        "CHEMBL3307636",
        "CHEMBL3307790",
    }
    assert list(prepared.columns) == CELL_LINE_BASE_COLUMNS


def test_prepare_cellline_dataframe__coerces_nullable_types() -> None:
    raw = pd.DataFrame(
        {
            "cell_chembl_id": ["CHEMBL1"],
            "cell_id": [pd.NA],
            "cell_source_tax_id": [9606],
            "cell_name": ["Example"],
        }
    )
    prepared, _ = prepare_cellline_dataframe(raw, ["CHEMBL1"])

    assert prepared["cell_id"].dtype == "Int64"
    assert prepared["cell_source_tax_id"].dtype == "Int64"
    assert all(
        prepared[col].dtype == "string"
        for col in [
            "cell_chembl_id",
            "cell_name",
            "cell_description",
            "cell_source_organism",
            "cell_source_tissue",
            "cellosaurus_id",
            "cl_lincs_id",
            "clo_id",
            "efo_id",
        ]
    )
