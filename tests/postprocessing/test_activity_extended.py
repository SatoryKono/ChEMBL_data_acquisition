"""Tests for :mod:`library.postprocessing.activity_extended`."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from library.postprocessing import ActivityExtendedError, process_activity_extended
from library.postprocessing.activity_extended import _FINAL_COLUMN_ORDER


def _write_activity_export(path: Path, rows: list[dict[str, object]]) -> None:
    frame = pd.DataFrame(rows)
    frame.to_csv(path, index=False, encoding="utf-8")


def _write_citation_fraction(path: Path) -> None:
    data = (
        "N,K_min_significant,test_used_at_threshold,p_value_at_threshold\n"
        "2,1,Fisher,0.05\n"
        "3,2,Fisher,0.05\n"
    )
    path.write_text(data, encoding="utf-8")


def _write_targets_dictionary(path: Path) -> None:
    rows = [
        {
            "target_chembl_id": "TAR1",
            "target_sort_order": "100",
            "multifunctional_enzyme": "FALSE",
            "IUPHAR_class": "ClassA",
            "IUPHAR_subclass": "SubclassA",
            "genus": "Homo",
            "superkingdom": "Eukaryota",
            "phylum": "Chordata",
            "taxon_id": "9606",
            "gene_index": "1",
            "taxon_index": "10",
        },
        {
            "target_chembl_id": "TAR2",
            "target_sort_order": "200",
            "multifunctional_enzyme": "TRUE",
            "IUPHAR_class": "ClassB",
            "IUPHAR_subclass": "SubclassB",
            "genus": "Escherichia",
            "superkingdom": "Bacteria",
            "phylum": "Proteobacteria",
            "taxon_id": "562",
            "gene_index": "2",
            "taxon_index": "20",
        },
    ]
    pd.DataFrame(rows).to_csv(path, index=False, encoding="utf-8")


@pytest.mark.postprocessing
def test_process_activity_extended__happy_path(tmp_path: Path) -> None:
    search_dir = tmp_path / "exports"
    search_dir.mkdir()
    dictionary_dir = tmp_path / "dictionary"
    (dictionary_dir / "_activity").mkdir(parents=True)
    (dictionary_dir / "_target").mkdir(parents=True)

    export_path = search_dir / "output.activity_20250101.csv"
    _write_activity_export(
        export_path,
        [
            {
                "activity_chembl_id": "ACT1",
                "salt_chembl_id": "SALT1",
                "molecule_chembl_id": "MOL1",
                "target_chembl_id": "TAR1",
                "assay_chembl_id": "ASSAY1",
                "document_chembl_id": "DOC1",
                "bao_endpoint": "EP1",
                "standard_type": "IC50",
                "standard_value": 10.0,
                "log_value": 5.0,
                "bao_format": "BAO_0001",
                "compound_key": "cmpd1",
                "compound_name": "Compound 1",
                "multmol_assay": pd.NA,
                "approx_cited_activity": 0,
                "shuffled_cit": 0,
                "exact_cited_activity": 1,
                "higly_correlated_cit": 0,
                "review_doc": 0,
                "rounded_data_citation": 0,
                "original_activity_approx": "approx",
                "original_activity_exact": "exact",
                "nstereo": 1,
            },
            {
                "activity_chembl_id": "ACT2",
                "salt_chembl_id": "SALT1",
                "molecule_chembl_id": "MOL1",
                "target_chembl_id": "TAR1",
                "assay_chembl_id": "ASSAY1",
                "document_chembl_id": "DOC1",
                "bao_endpoint": "EP1",
                "standard_type": "IC50",
                "standard_value": 11.0,
                "log_value": 5.1,
                "bao_format": "BAO_0001",
                "compound_key": "cmpd1",
                "compound_name": "Compound 1",
                "multmol_assay": pd.NA,
                "approx_cited_activity": 0,
                "shuffled_cit": 0,
                "exact_cited_activity": 0,
                "higly_correlated_cit": 0,
                "review_doc": 0,
                "rounded_data_citation": 0,
                "original_activity_approx": "approx",
                "original_activity_exact": "exact",
                "nstereo": 1,
            },
            {
                "activity_chembl_id": "ACT3",
                "salt_chembl_id": "SALT2",
                "molecule_chembl_id": "MOL2",
                "target_chembl_id": "TAR2",
                "assay_chembl_id": "ASSAY2",
                "document_chembl_id": "DOC2",
                "bao_endpoint": "EP2",
                "standard_type": "Ki",
                "standard_value": 7.0,
                "log_value": 4.9,
                "bao_format": "BAO_0002",
                "compound_key": "cmpd2",
                "compound_name": "Compound 2",
                "multmol_assay": "TRUE",
                "approx_cited_activity": 0,
                "shuffled_cit": 0,
                "exact_cited_activity": 0,
                "higly_correlated_cit": 0,
                "review_doc": 0,
                "rounded_data_citation": 0,
                "original_activity_approx": "approx",
                "original_activity_exact": "exact",
                "nstereo": 2,
            },
        ],
    )

    _write_citation_fraction(dictionary_dir / "_activity" / "citation_fraction.csv")
    _write_targets_dictionary(dictionary_dir / "_target" / "targets_type.csv")

    output_path = process_activity_extended(
        search_dir=search_dir,
        dictionary_dir=dictionary_dir,
    )

    assert output_path == search_dir / "extended.output.activity_20250101.csv"
    assert output_path.exists()

    result = pd.read_csv(output_path)
    expected_columns = list(_FINAL_COLUMN_ORDER)
    assert result.columns.tolist() == expected_columns

    assert result["saltform_id"].tolist() == ["SALT1", "SALT1", "SALT2"]
    assert result["unknown_chirality"].tolist() == [False, False, True]
    assert result["multmol_assay"].tolist() == [True, True, True]
    assert result["is_citation"].tolist() == [True, False, False]
    assert result["high_citation_rate"].tolist() == [True, True, False]
    assert result["unicellular_organism"].tolist() == [False, False, True]
    assert result["multifunctional_enzyme"].tolist() == [False, False, True]


@pytest.mark.postprocessing
def test_process_activity_extended__missing_dictionary(tmp_path: Path) -> None:
    search_dir = tmp_path / "exports"
    search_dir.mkdir()
    export_path = search_dir / "output.activity_20250101.csv"
    _write_activity_export(
        export_path,
        [
            {
                "activity_chembl_id": "ACT1",
                "salt_chembl_id": "SALT1",
                "molecule_chembl_id": "MOL1",
                "target_chembl_id": "TAR1",
                "assay_chembl_id": "ASSAY1",
                "document_chembl_id": "DOC1",
                "bao_endpoint": "EP1",
                "standard_type": "IC50",
                "standard_value": 1.0,
                "log_value": 5.0,
                "bao_format": "BAO_0001",
                "compound_key": "cmpd1",
                "compound_name": "Compound 1",
                "multmol_assay": pd.NA,
                "approx_cited_activity": 0,
                "shuffled_cit": 0,
                "exact_cited_activity": 0,
                "higly_correlated_cit": 0,
                "review_doc": 0,
                "rounded_data_citation": 0,
                "original_activity_approx": "approx",
                "original_activity_exact": "exact",
                "nstereo": 1,
            }
        ],
    )

    dictionary_dir = tmp_path / "dictionary"
    dictionary_dir.mkdir()

    with pytest.raises(ActivityExtendedError) as excinfo:
        process_activity_extended(search_dir=search_dir, dictionary_dir=dictionary_dir)

    assert "citation_fraction.csv" in str(excinfo.value)


@pytest.mark.postprocessing
def test_process_activity_extended__missing_column(tmp_path: Path) -> None:
    search_dir = tmp_path / "exports"
    search_dir.mkdir()
    export_path = search_dir / "output.activity_20250101.csv"
    # ``log_value`` column intentionally omitted.
    _write_activity_export(
        export_path,
        [
            {
                "activity_chembl_id": "ACT1",
                "salt_chembl_id": "SALT1",
                "molecule_chembl_id": "MOL1",
                "target_chembl_id": "TAR1",
                "assay_chembl_id": "ASSAY1",
                "document_chembl_id": "DOC1",
                "bao_endpoint": "EP1",
                "standard_type": "IC50",
                "standard_value": 1.0,
                "bao_format": "BAO_0001",
                "compound_key": "cmpd1",
                "compound_name": "Compound 1",
                "multmol_assay": pd.NA,
                "approx_cited_activity": 0,
                "shuffled_cit": 0,
                "exact_cited_activity": 0,
                "higly_correlated_cit": 0,
                "review_doc": 0,
                "rounded_data_citation": 0,
                "original_activity_approx": "approx",
                "original_activity_exact": "exact",
                "nstereo": 1,
            }
        ],
    )

    dictionary_dir = tmp_path / "dictionary"
    (dictionary_dir / "_activity").mkdir(parents=True)
    (dictionary_dir / "_target").mkdir(parents=True)
    _write_citation_fraction(dictionary_dir / "_activity" / "citation_fraction.csv")
    _write_targets_dictionary(dictionary_dir / "_target" / "targets_type.csv")

    with pytest.raises(ActivityExtendedError) as excinfo:
        process_activity_extended(search_dir=search_dir, dictionary_dir=dictionary_dir)

    assert "log_value" in str(excinfo.value)
