"""Unit tests for UniProt-specific helpers in ``scripts.get_target_data``."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
import pytest
from scripts import get_target_data


@pytest.mark.unit
def test_fetch_uniprot__respects_final_out_path(
    tmp_path: Path, cfg: get_target_data.Config, monkeypatch: pytest.MonkeyPatch
) -> None:
    """``fetch_uniprot`` must honour the provided ``output_csv`` path."""

    output_csv = tmp_path / "exports" / "custom_uniprot.csv"
    chembl_df = pd.DataFrame(
        {
            "target_chembl_id": ["CHEMBL1"],
            "uniprot_id": ["P12345"],
        }
    )

    cfg.io.output_dir = tmp_path / "out"
    cfg.target.uniprot.data_dir = tmp_path / "uniprot_data"
    cfg.target.all.data_dir = tmp_path / "uniprot_data"
    cfg.target.uniprot.data_dir.mkdir(parents=True, exist_ok=True)
    cfg.target.all.data_dir.mkdir(parents=True, exist_ok=True)

    expected_frame = pd.DataFrame(
        {
            "uniprot_id": ["P12345"],
            "original_id": ["P12345"],
            "source_column": ["uniprot_id"],
            "mapping_uniprot_id": [""],
        }
    )

    calls: dict[str, Path] = {}

    def fake_run_uniprot(
        local_cfg: get_target_data.Config, args: argparse.Namespace
    ) -> int:
        final_out = Path(args.final_out)
        calls["final_out"] = final_out
        calls["input_csv"] = Path(args.input_csv)

        final_out.parent.mkdir(parents=True, exist_ok=True)
        expected_frame.to_csv(
            final_out,
            index=False,
            sep=local_cfg.io.csv_sep,
            encoding=local_cfg.io.csv_encoding,
        )
        return 0

    monkeypatch.setattr(get_target_data, "run_uniprot", fake_run_uniprot)
    result = get_target_data.fetch_uniprot(cfg, chembl_df, output_csv)

    assert calls["final_out"] == output_csv
    assert calls["input_csv"].is_absolute()
    pd.testing.assert_frame_equal(result, expected_frame[result.columns])
