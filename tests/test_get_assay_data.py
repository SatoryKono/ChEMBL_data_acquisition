from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
import pytest

import library.cli_utils as cli_utils
from library import chembl_library as cl
from library import io
from library.config import Config
from scripts import get_assay_data as gas


def test_run_chembl_orders_columns(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, cfg: Config
) -> None:
    """Output columns should follow `AssaysSchema` order with extras alphabetic."""
    input_csv = tmp_path / "assays.csv"
    input_csv.write_text("assay_chembl_id\nA1\n")
    args = argparse.Namespace(input_csv=input_csv, output_csv=tmp_path / "out.csv")

    monkeypatch.setattr(io, "read_ids", lambda *_, **__: iter(["A1"]))

    df = pd.DataFrame(
        {
            "zzz": ["z"],
            "document_chembl_id": ["D1"],
            "aaa": ["a"],
            "assay_chembl_id": ["A1"],
            "target_chembl_id": ["T1"],
        }
    )
    monkeypatch.setattr(cl, "get_assays", lambda *_, **__: df)
    monkeypatch.setattr(gas.ap, "postprocess_assays", lambda df: df)
    monkeypatch.setattr(gas, "analyze_table_quality", lambda df, table_name: None)
    monkeypatch.setattr(cli_utils, "write_meta_yaml", lambda **kwargs: None)
    monkeypatch.setattr(cli_utils, "file_sha256", lambda p: "deadbeef")

    rc = gas.run_chembl(cfg, args)
    assert rc == 0
    out_df = pd.read_csv(args.output_csv)
    assert list(out_df.columns) == [
        "assay_chembl_id",
        "document_chembl_id",
        "target_chembl_id",
        "pipeline_version",
        "timestamp_utc",
        "aaa",
        "zzz",
    ]
