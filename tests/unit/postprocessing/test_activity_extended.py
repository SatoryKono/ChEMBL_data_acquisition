from pathlib import Path

import pandas as pd

from library.postprocessing.activity_extended import (
    _augment_activity_frame,
    _load_target_metadata,
)


def test_augment_activity_frame__fills_missing_log_value():
    frame = pd.DataFrame(
        {
            "log_value": pd.Series(["not-a-number", pd.NA], dtype="string"),
            "standard_value": [1.0, 1.0],
            "standard_units": ["nM", "nM"],
        }
    )

    result, _ = _augment_activity_frame(frame)

    expected = pd.Series([9.0, 9.0], dtype="Float64", name="log_value")

    assert str(result["log_value"].dtype) == "Float64"
    pd.testing.assert_series_equal(result["log_value"], expected)


def test_augment_activity_frame__creates_log_value_from_standard_measurements():
    frame = pd.DataFrame(
        {
            "standard_value": [1.0, 10.0, pd.NA],
            "standard_units": ["nM", "nM", "nM"],
        }
    )

    result, filled = _augment_activity_frame(frame)

    expected = pd.Series([9.0, 8.0, pd.NA], dtype="Float64", name="log_value")

    assert "log_value" in filled
    assert str(result["log_value"].dtype) == "Float64"
    pd.testing.assert_series_equal(result["log_value"], expected)


def test_load_target_metadata__infers_unicellular_from_type(tmp_path: Path) -> None:
    csv_path = tmp_path / "targets_type.csv"
    frame = pd.DataFrame(
        {
            "target_chembl_id": ["CHEMBL1", "CHEMBL2"],
            "target_sort_order": ["001", "002"],
            "multifunctional_enzyme": ["TRUE", "FALSE"],
            "IUPHAR_class": ["0100", "0200"],
            "IUPHAR_subclass": ["0100-1", "0200-1"],
            "genus": ["GenusA", "GenusB"],
            "superkingdom": ["Sk1", "Sk2"],
            "phylum": ["Ph1", "Ph2"],
            "taxon_id": [1, 2],
            "gene_index": ["GENE1", "GENE2"],
            "taxon_index": ["IDX1", "IDX2"],
            "type": ["Unicellular organism", "Multicellular organism"],
        }
    )
    frame.to_csv(csv_path, index=False, encoding="cp1252")

    loaded = _load_target_metadata(csv_path)

    assert "unicellular_organism" in loaded.columns
    assert loaded["unicellular_organism"].dtype == "boolean"
    assert loaded["unicellular_organism"].tolist() == [True, False]


def test_load_target_metadata__normalises_existing_boolean_column(tmp_path: Path) -> None:
    csv_path = tmp_path / "targets_type.csv"
    frame = pd.DataFrame(
        {
            "target_chembl_id": ["CHEMBL1", "CHEMBL2"],
            "target_sort_order": ["001", "002"],
            "multifunctional_enzyme": ["FALSE", "TRUE"],
            "IUPHAR_class": ["0100", "0200"],
            "IUPHAR_subclass": ["0100-1", "0200-1"],
            "genus": ["GenusA", "GenusB"],
            "superkingdom": ["Sk1", "Sk2"],
            "phylum": ["Ph1", "Ph2"],
            "taxon_id": [1, 2],
            "gene_index": ["GENE1", "GENE2"],
            "taxon_index": ["IDX1", "IDX2"],
            "Unicellular organism": ["TRUE", "FALSE"],
        }
    )
    frame.to_csv(csv_path, index=False, encoding="cp1252")

    loaded = _load_target_metadata(csv_path)

    assert loaded["unicellular_organism"].dtype == "boolean"
    assert loaded["unicellular_organism"].tolist() == [True, False]


def test_load_target_metadata__fills_missing_unicellular_column(tmp_path: Path) -> None:
    csv_path = tmp_path / "targets_type.csv"
    frame = pd.DataFrame(
        {
            "target_chembl_id": ["CHEMBL1", "CHEMBL2"],
            "target_sort_order": ["001", "002"],
            "multifunctional_enzyme": ["FALSE", "TRUE"],
            "IUPHAR_class": ["0100", "0200"],
            "IUPHAR_subclass": ["0100-1", "0200-1"],
            "genus": ["GenusA", "GenusB"],
            "superkingdom": ["Sk1", "Sk2"],
            "phylum": ["Ph1", "Ph2"],
            "taxon_id": [1, 2],
            "gene_index": ["GENE1", "GENE2"],
            "taxon_index": ["IDX1", "IDX2"],
        }
    )
    frame.to_csv(csv_path, index=False, encoding="cp1252")

    loaded = _load_target_metadata(csv_path)

    assert "unicellular_organism" in loaded.columns
    assert loaded["unicellular_organism"].dtype == "boolean"
    assert loaded["unicellular_organism"].isna().all()
