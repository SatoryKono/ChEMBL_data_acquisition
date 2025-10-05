"""Integration tests for the activity extended post-processing pipeline."""

from __future__ import annotations

import csv
import hashlib
import logging
import os
import random
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from library.postprocessing.activity_extended import process_activity_extended


def _fix_seed(seed: int = 42) -> None:
    os.environ["PYTHONHASHSEED"] = str(seed)
    random.seed(seed)
    np.random.seed(seed)


def _prepare_dictionary(dictionary_root: Path) -> None:
    activity_subdir = dictionary_root / "_activity"
    activity_subdir.mkdir(parents=True, exist_ok=True)
    with (activity_subdir / "citation_fraction.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["N", "K_min_significant"])
        writer.writeheader()
        writer.writerow({"N": "1", "K_min_significant": "1"})

    targets_subdir = dictionary_root / "_target"
    targets_subdir.mkdir(parents=True, exist_ok=True)
    with (targets_subdir / "targets_type.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "target_chembl_id",
                "target_sort_order",
                "multifunctional_enzyme",
                "IUPHAR_class",
                "IUPHAR_subclass",
                "genus",
                "superkingdom",
                "phylum",
                "taxon_id",
                "gene_index",
                "taxon_index",
            ],
        )
        writer.writeheader()
        writer.writerow(
            {
                "target_chembl_id": "CHEMBLT1",
                "target_sort_order": "1",
                "multifunctional_enzyme": "false",
                "IUPHAR_class": "ClassA",
                "IUPHAR_subclass": "SubclassA",
                "genus": "Candida",
                "superkingdom": "Eukaryota",
                "phylum": "Ascomycota",
                "taxon_id": "1234",
                "gene_index": "GENE1",
                "taxon_index": "TAX1",
            }
        )


@pytest.mark.integration
def test_process_activity_extended__fills_missing_optional_columns(tmp_path, caplog):
    """Process minimal export missing optional columns and emit backfill warning."""

    _fix_seed()

    activity_dir = tmp_path / "input"
    dictionary_root = tmp_path / "dictionary"
    activity_dir.mkdir()

    activity_path = activity_dir / "output.activities_20240101.csv"
    activity_columns = [
        "activity_id",
        "molecule_chembl_id",
        "target_chembl_id",
        "assay_chembl_id",
        "document_chembl_id",
        "bao_endpoint",
        "standard_type",
        "standard_value",
        "bao_format",
        "parent_molecule_chembl_id",
        "molecule_pref_name",
        "pchembl_value",
    ]
    with activity_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=activity_columns)
        writer.writeheader()
        writer.writerow(
            {
                "activity_id": "1",
                "molecule_chembl_id": "CHEMBL1",
                "target_chembl_id": "CHEMBLT1",
                "assay_chembl_id": "ASSAY1",
                "document_chembl_id": "DOC1",
                "bao_endpoint": "endpoint",
                "standard_type": "IC50",
                "standard_value": "1.0",
                "bao_format": "format",
                "parent_molecule_chembl_id": "CHEMBL1",
                "molecule_pref_name": "Example compound",
                "pchembl_value": "7.0",
            }
        )

    _prepare_dictionary(dictionary_root)

    with caplog.at_level(logging.WARNING):
        output_path = process_activity_extended(
            input_path=activity_path,
            dictionary_dir=dictionary_root,
        )

    assert output_path.exists()

    filled_warnings = [
        record for record in caplog.records if "activity_extended_missing_columns_filled" in record.message
    ]
    assert filled_warnings, "Expected warning about backfilled columns not emitted"

    message = filled_warnings[0].message
    for column_name in ("activity_chembl_id", "compound_key", "log_value"):
        assert column_name in message


@pytest.mark.integration
def test_process_activity_extended__deduplicates_and_deterministic(tmp_path, caplog):
    """Remove duplicate activities and produce deterministic output."""

    _fix_seed()

    activity_dir = tmp_path / "input"
    dictionary_root = tmp_path / "dictionary"
    activity_dir.mkdir()

    _prepare_dictionary(dictionary_root)

    activity_path = activity_dir / "output.activities_20240102.csv"
    activity_columns = [
        "activity_id",
        "molecule_chembl_id",
        "target_chembl_id",
        "assay_chembl_id",
        "document_chembl_id",
        "bao_endpoint",
        "standard_type",
        "standard_value",
        "bao_format",
        "parent_molecule_chembl_id",
        "molecule_pref_name",
        "pchembl_value",
    ]
    records = [
        {
            "activity_id": "101",
            "molecule_chembl_id": "CHEMBL1",
            "target_chembl_id": "CHEMBLT1",
            "assay_chembl_id": "ASSAY1",
            "document_chembl_id": "DOC1",
            "bao_endpoint": "endpoint",
            "standard_type": "IC50",
            "standard_value": "5.0",
            "bao_format": "format",
            "parent_molecule_chembl_id": "CHEMBL1",
            "molecule_pref_name": "Compound One",
            "pchembl_value": "7.0",
        },
        {
            "activity_id": "101",
            "molecule_chembl_id": "CHEMBL1",
            "target_chembl_id": "CHEMBLT1",
            "assay_chembl_id": "ASSAY1",
            "document_chembl_id": "DOC1",
            "bao_endpoint": "endpoint",
            "standard_type": "IC50",
            "standard_value": "5.0",
            "bao_format": "format",
            "parent_molecule_chembl_id": "CHEMBL1",
            "molecule_pref_name": "Compound Duplicate",
            "pchembl_value": "7.1",
        },
        {
            "activity_id": "102",
            "molecule_chembl_id": "CHEMBL2",
            "target_chembl_id": "CHEMBLT1",
            "assay_chembl_id": "ASSAY2",
            "document_chembl_id": "DOC2",
            "bao_endpoint": "endpoint",
            "standard_type": "EC50",
            "standard_value": "1.5",
            "bao_format": "format",
            "parent_molecule_chembl_id": "CHEMBL2",
            "molecule_pref_name": "Compound Two",
            "pchembl_value": "6.0",
        },
    ]
    with activity_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=activity_columns)
        writer.writeheader()
        writer.writerows(records)

    with caplog.at_level(logging.INFO):
        output_path = process_activity_extended(
            input_path=activity_path,
            dictionary_dir=dictionary_root,
        )

    assert output_path.exists()

    output_df = pd.read_csv(output_path)
    subset = [
        "activity_chembl_id",
        "assay_chembl_id",
        "document_chembl_id",
        "standard_type",
        "standard_value",
    ]
    assert set(subset).issubset(output_df.columns)
    assert len(output_df) == 2
    assert not output_df.duplicated(subset=subset).any()

    expected_sorted = output_df.sort_values(by=subset, kind="mergesort").reset_index(drop=True)
    pd.testing.assert_frame_equal(output_df.reset_index(drop=True), expected_sorted)

    dedupe_records = [
        record for record in caplog.records if "activity_extended_deduplicated" in record.message
    ]
    assert dedupe_records, "Expected deduplication log not emitted"
    assert "removed=1" in dedupe_records[0].message

    first_hash = hashlib.sha256(output_path.read_bytes()).hexdigest()

    caplog.clear()
    second_output = process_activity_extended(
        input_path=activity_path,
        dictionary_dir=dictionary_root,
    )
    assert second_output == output_path
    second_hash = hashlib.sha256(output_path.read_bytes()).hexdigest()
    assert second_hash == first_hash
