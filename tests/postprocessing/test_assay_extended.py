"""Tests for :mod:`library.postprocessing.assay_extended`."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from library.postprocessing import enrich_assay_metadata


@pytest.mark.postprocessing
@pytest.mark.usefixtures("deterministic_env")
def test_enrich_assay_metadata__fills_missing_fields(tmp_path: Path) -> None:
    dictionary_root = tmp_path
    assay_dir = dictionary_root / "_assay"
    assay_dir.mkdir(parents=True)
    taxonomy_dir = dictionary_root / "_taxonomy"
    taxonomy_dir.mkdir(parents=True)
    document_dir = dictionary_root / "_document"
    document_dir.mkdir(parents=True)
    target_dir = dictionary_root / "_target"
    target_dir.mkdir(parents=True)

    pd.DataFrame(
        {
            "assay_chembl_id": ["CHEMBL1", "CHEMBL2"],
            "assay_tax_id": ["9606", "10090"],
            "assay_group": ["FallbackGroup", "Rodent"],
            "assay_strain": ["FallbackStrain", "Mouse"],
            "document_chembl_id": ["DOC-1", "DOC-2"],
            "target_chembl_id": ["CHEMBLT1", "CHEMBLT2"],
        }
    ).to_csv(assay_dir / "assay.csv", index=False)

    pd.DataFrame(
        {
            "tax_id": ["9606", "10090"],
            "assay_group": ["Human", "Rodent"],
            "assay_strain": ["HEK293", "Balb/c"],
        }
    ).to_csv(taxonomy_dir / "taxonomy.csv", index=False)

    pd.DataFrame(
        {
            "document_chembl_id": ["DOC-1", "DOC-2"],
            "year": [1999, 2010],
        }
    ).to_csv(document_dir / "document.csv", index=False)

    pd.DataFrame(
        {
            "target_chembl_id": ["CHEMBLT1", "CHEMBLT2"],
            "target_components": [
                "[{\"accession\": \"Q99999\"}, {\"accession\": \"P12345\"}]",
                "[]",
            ],
        }
    ).to_csv(target_dir / "output.target_20240101.csv", index=False)

    frame = pd.DataFrame(
        {
            "assay_chembl_id": ["CHEMBL1", "CHEMBL2"],
            "assay_tax_id": [pd.NA, "10090"],
            "assay_group": [pd.NA, "CustomGroup"],
            "assay_strain": [pd.NA, pd.NA],
            "document_chembl_id": ["DOC-1", pd.NA],
            "target_chembl_id": [pd.NA, "CHEMBLT2"],
            "accession": [pd.NA, "ExistingAcc"],
        }
    )

    result = enrich_assay_metadata(frame, dictionary_dir=dictionary_root)

    assert result.loc[0, "assay_group"] == "Human"
    assert result.loc[0, "assay_strain"] == "HEK293"
    assert result.loc[0, "assay_tax_id"] == "9606"
    assert result.loc[0, "target_chembl_id"] == "CHEMBLT1"
    assert result.loc[0, "accession"] == "P12345|Q99999"
    assert result.loc[0, "year"] == 1999

    assert result.loc[1, "assay_group"] == "CustomGroup"
    assert result.loc[1, "assay_strain"] == "Balb/c"
    assert result.loc[1, "document_chembl_id"] == "DOC-2"
    assert result.loc[1, "accession"] == "ExistingAcc"
    assert result.loc[1, "year"] == 2010

    assert str(result["year"].dtype) == "Int64"


@pytest.mark.postprocessing
def test_enrich_assay_metadata__missing_taxonomy_dir_logs_warning(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    dictionary_root = tmp_path
    (dictionary_root / "_assay").mkdir()
    (dictionary_root / "_document").mkdir()
    (dictionary_root / "_target").mkdir()

    pd.DataFrame({"assay_chembl_id": ["CHEMBL1"]}).to_csv(
        dictionary_root / "_assay" / "assay.csv", index=False
    )
    pd.DataFrame({"document_chembl_id": ["DOC-1"], "year": [2000]}).to_csv(
        dictionary_root / "_document" / "document.csv", index=False
    )
    pd.DataFrame(
        {"target_chembl_id": ["CHEMBLT1"], "target_components": ["[]"]}
    ).to_csv(dictionary_root / "_target" / "output.target_20240101.csv", index=False)

    frame = pd.DataFrame({"assay_chembl_id": ["CHEMBL1"]})

    with caplog.at_level("WARNING"):
        result = enrich_assay_metadata(frame, dictionary_dir=dictionary_root)

    assert "assay_extended_missing_taxonomy_dictionary" in caplog.text
    assert list(result.columns) == [
        "assay_chembl_id",
        "assay_tax_id",
        "assay_group",
        "assay_strain",
        "target_chembl_id",
        "document_chembl_id",
        "year",
        "accession",
    ]
    assert result.loc[0, "assay_chembl_id"] == "CHEMBL1"
    assert result["assay_tax_id"].isna().all()
    assert result["assay_group"].isna().all()
    assert result["assay_strain"].isna().all()
    assert result["target_chembl_id"].isna().all()
    assert result["document_chembl_id"].isna().all()
    assert result["accession"].isna().all()
    assert result["year"].isna().all()
    assert str(result["assay_tax_id"].dtype) == "string"
    assert str(result["assay_group"].dtype) == "string"
    assert str(result["assay_strain"].dtype) == "string"
    assert str(result["accession"].dtype) == "string"
    assert str(result["year"].dtype) == "Int64"
