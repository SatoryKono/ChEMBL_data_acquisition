from pathlib import Path

import pandas as pd
import pytest

from library.pipelines.testitem import cli


@pytest.mark.unit
def test_read_input_ids__load_and_validate(sample_input_csv: Path, cfg) -> None:
    frame = pd.read_csv(sample_input_csv)
    frame["value"] = [10, 20]
    frame.loc[0, "molecule_chembl_id"] = " CHEMBL1 "
    frame.to_csv(sample_input_csv, index=False)

    rc, result = cli.read_input_ids(
        sample_input_csv,
        column="molecule_chembl_id",
        io_cfg=cfg.io,
        limit=None,
        offset=0,
    )

    assert rc == 0
    assert result is not None
    ids = list(result.ids_iter)
    assert ids == ["CHEMBL1", "CHEMBL2"]
    assert result.sample_ids == ("CHEMBL1", "CHEMBL2")


@pytest.mark.unit
def test_read_input_ids__missing_file(cfg) -> None:
    rc, result = cli.read_input_ids(
        Path("does-not-exist.csv"),
        column="molecule_chembl_id",
        io_cfg=cfg.io,
        limit=None,
        offset=0,
    )

    assert rc == 1
    assert result is None


@pytest.mark.unit
def test_read_input_ids__limit_and_offset(sample_input_csv: Path, cfg) -> None:
    frame = pd.DataFrame(
        {
            "molecule_chembl_id": [f"CHEMBL{i}" for i in range(6)],
            "value": range(6),
        }
    )
    frame.to_csv(sample_input_csv, index=False)

    rc, result = cli.read_input_ids(
        sample_input_csv,
        column="molecule_chembl_id",
        io_cfg=cfg.io,
        limit=2,
        offset=3,
    )

    assert rc == 0
    assert result is not None
    ids = list(result.ids_iter)
    assert ids == ["CHEMBL3", "CHEMBL4"]
    assert result.sample_ids == ("CHEMBL3", "CHEMBL4")


@pytest.mark.unit
def test_read_input_ids__invalid_column(sample_input_csv: Path, cfg) -> None:
    frame = pd.DataFrame({"other_column": ["CHEMBL1"]})
    frame.to_csv(sample_input_csv, index=False)

    rc, result = cli.read_input_ids(
        sample_input_csv,
        column="molecule_chembl_id",
        io_cfg=cfg.io,
        limit=None,
        offset=0,
    )

    assert rc == 1
    assert result is None
