import csv
from pathlib import Path

import pandas as pd
import pytest
from library.pipelines.testitem import cli


@pytest.mark.unit
def test_read_input_ids__load_and_validate(tmp_path: Path, cfg) -> None:
    input_path = tmp_path / "input.csv"
    with input_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["molecule_chembl_id", "value"])
        writer.writerow([" CHEMBL1 ", "10"])
        writer.writerow(["CHEMBL2", "20"])

    rc, result = cli.read_input_ids(
        input_path,
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
def test_read_input_ids__limit_and_offset(tmp_path: Path, cfg) -> None:
    frame = pd.DataFrame(
        {
            "molecule_chembl_id": [f"CHEMBL{i}" for i in range(6)],
            "value": range(6),
        }
    )
    input_path = tmp_path / "input.csv"
    frame.to_csv(input_path, index=False)

    rc, result = cli.read_input_ids(
        input_path,
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
def test_read_input_ids__invalid_column(tmp_path: Path, cfg) -> None:
    input_path = tmp_path / "input.csv"
    with input_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["other_column"])
        writer.writerow(["CHEMBL1"])

    rc, result = cli.read_input_ids(
        input_path,
        column="molecule_chembl_id",
        io_cfg=cfg.io,
        limit=None,
        offset=0,
    )

    assert rc == 1
    assert result is None
