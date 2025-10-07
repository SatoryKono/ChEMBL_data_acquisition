"""Tests for :mod:`library.postprocessing.assay_extended`."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from library.postprocessing import AssayExtendedError, enrich_assay_metadata


pytestmark = pytest.mark.integration


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

def test_enrich_assay_metadata__missing_taxonomy_dir_raises(tmp_path: Path) -> None:
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

    with pytest.raises(AssayExtendedError):
        enrich_assay_metadata(frame, dictionary_dir=dictionary_root)
