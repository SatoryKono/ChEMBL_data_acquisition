from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
import pytest

from library.config import Config
from scripts import get_assay_data
from tests.helpers import ASSAY_ENRICHMENT_MIN_RATIO


def _ensure_parent(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)


@pytest.mark.integration
def test_get_assay_cli__enrichment_quality(
    tmp_path: Path, cfg: Config, monkeypatch: pytest.MonkeyPatch
) -> None:
    data_dir = Path(__file__).resolve().parents[1] / "data"
    input_csv = tmp_path / "assay.csv"
    input_csv.write_text((data_dir / "assay.csv").read_text(encoding="utf-8"), encoding="utf-8")
    output_csv = tmp_path / "out" / "assays.csv"
    dictionary_path = data_dir / "assay_dictionary.csv"

    def _stub_run_chembl(config: Config, args: argparse.Namespace) -> int:
        frame = pd.read_csv(args.input_csv)
        dictionary = pd.read_csv(dictionary_path)
        dictionary["assay_chembl_id"] = dictionary["assay_chembl_id"].astype("string")
        enriched = frame.merge(dictionary, on="assay_chembl_id", how="left")
        enriched["description"] = enriched["description"].astype("string").str.strip()
        enriched["description_length"] = enriched["description"].str.len().astype("Int64")
        enriched["year"] = pd.to_numeric(enriched["year"], errors="coerce").astype("Int64")
        quality_columns = ["assay_strain", "assay_group", "year", "accession"]
        completeness = 1.0 - enriched[quality_columns].isna().mean()
        if float(completeness.min()) < ASSAY_ENRICHMENT_MIN_RATIO:
            raise AssertionError(
                "assay enrichment below threshold "
                f"(threshold={ASSAY_ENRICHMENT_MIN_RATIO}, completeness={completeness.to_dict()})"
            )
        _ensure_parent(args.final_out)
        columns = [
            "assay_chembl_id",
            "target_chembl_id",
            "document_chembl_id",
            "description",
            "description_length",
            "assay_strain",
            "assay_group",
            "year",
            "accession",
        ]
        enriched.to_csv(args.final_out, index=False, columns=columns)
        return 0

    monkeypatch.setattr(get_assay_data, "run_chembl", _stub_run_chembl)

    args = argparse.Namespace(
        input_csv=input_csv,
        final_out=output_csv,
        skip_existing=False,
        force=False,
    )

    exit_code = get_assay_data.run(cfg, args)

    assert exit_code == 0
    result = pd.read_csv(output_csv)
    expected_columns = [
        "assay_chembl_id",
        "target_chembl_id",
        "document_chembl_id",
        "description",
        "description_length",
        "assay_strain",
        "assay_group",
        "year",
        "accession",
    ]
    assert list(result.columns) == expected_columns
    assert (result["description_length"] == result["description"].str.len()).all()
    quality_columns = ["assay_strain", "assay_group", "year", "accession"]
    completeness = 1.0 - result[quality_columns].isna().mean()
    assert (completeness >= ASSAY_ENRICHMENT_MIN_RATIO).all(), completeness.to_dict()
