"""Integration tests for target pipeline IUPHAR integration."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest
from scripts import get_target_data

from library.config import Config


@pytest.mark.integration
def test_target_pipeline__iuphar_columns_populated(
    tmp_path: Path, cfg: Config, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Ensure the pipeline writes the requested IUPHAR output path."""

    cfg.io.csv_sep = ","
    cfg.io.csv_encoding = "utf-8"
    cfg.target.all.uniprot_column = "uniprot_id"
    cfg.target.iuphar.target_csv = tmp_path / "iuphar_targets.csv"
    cfg.target.iuphar.family_csv = tmp_path / "iuphar_families.csv"
    cfg.target.all.target_csv = cfg.target.iuphar.target_csv
    cfg.target.all.family_csv = cfg.target.iuphar.family_csv

    chembl_df = pd.DataFrame(
        {
            "target_chembl_id": ["CHEMBL1"],
            "uniprot_id": ["P12345"],
            "mapping_uniprot_id": ["P12345"],
            "pref_name": ["Alpha"],
        }
    )
    uniprot_df = pd.DataFrame(
        {
            "uniprot_id": ["P12345"],
            "mapping_uniprot_id": ["P12345"],
            "molecular_function": ["binding"],
        }
    )

    output_csv = tmp_path / "output.targets_iuphar.csv"
    expected_payload = {
        "uniprot_id": ["P12345"],
        "iuphar_target_id": ["TGT-001"],
        "iuphar_family_id": ["FAM-123"],
        "iuphar_type": ["Type"],
        "iuphar_class": ["Class"],
        "iuphar_subclass": ["Subclass"],
        "iuphar_chain": ["Chain"],
        "iuphar_name": ["Example"],
        "iuphar_full_id_path": ["TGT-001#FAM-123"],
        "iuphar_full_name_path": ["Example#Family"],
    }
    captured_paths: list[Path] = []

    class _StubIUPHARData:
        def map_uniprot_file(
            self,
            input_path: str | Path,
            output_path: str | Path,
            *,
            encoding: str = "utf-8",
            sep: str = ",",
        ) -> pd.DataFrame:
            captured_paths.append(Path(output_path))
            frame = pd.DataFrame(expected_payload)
            Path(output_path).parent.mkdir(parents=True, exist_ok=True)
            frame.to_csv(output_path, index=False, encoding=encoding, sep=sep)
            return frame

    class _StubIUPHARFactory:
        @classmethod
        def from_files(
            cls,
            target_path: str | Path,
            family_path: str | Path,
            *,
            encoding: str = "utf-8",
        ) -> _StubIUPHARData:
            return _StubIUPHARData()

    monkeypatch.setattr(get_target_data.ii, "IUPHARData", _StubIUPHARFactory)

    combined_df, iuphar_df = get_target_data.fetch_iuphar(
        cfg, chembl_df, uniprot_df, output_csv
    )

    assert captured_paths == [output_csv]
    assert output_csv.exists()

    written = pd.read_csv(output_csv, dtype=str)
    for column, values in expected_payload.items():
        assert column in written.columns
        assert list(written[column]) == values
        assert column in iuphar_df.columns
        assert list(iuphar_df[column]) == values

    assert "uniprot_id" in combined_df.columns
    assert list(combined_df["uniprot_id"]) == ["P12345"]
