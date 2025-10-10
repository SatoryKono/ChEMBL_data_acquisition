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
    with (activity_subdir / "citation_fraction.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=["N", "K_min_significant"])
        writer.writeheader()
        writer.writerow({"N": "1", "K_min_significant": "1"})

    targets_subdir = dictionary_root / "_target"
    targets_subdir.mkdir(parents=True, exist_ok=True)
    with (targets_subdir / "targets_type.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
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
                "unicellular_organism",
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
                "unicellular_organism": "false",
            }
        )

    document_subdir = dictionary_root / "_document"
    document_subdir.mkdir(parents=True, exist_ok=True)
    with (document_subdir / "document.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.DictWriter(
            handle, fieldnames=["document_chembl_id", "completed", "review"]
        )
        writer.writeheader()
        writer.writerow(
            {
                "document_chembl_id": "DOC1",
                "completed": "1980-08-15",
                "review": "false",
            }
        )
        writer.writerow(
            {
                "document_chembl_id": "DOC2",
                "completed": "2008-01-11",
                "review": "false",
            }
        )

    assay_subdir = dictionary_root / "_assay"
    assay_subdir.mkdir(parents=True, exist_ok=True)
    with (assay_subdir / "assay.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=["assay_chembl_id", "assay_with_same_target"]
        )
        writer.writeheader()
        writer.writerow({"assay_chembl_id": "ASSAY1", "assay_with_same_target": "3"})
        writer.writerow({"assay_chembl_id": "ASSAY2", "assay_with_same_target": "2"})

    testitem_subdir = dictionary_root / "_testitem"
    testitem_subdir.mkdir(parents=True, exist_ok=True)
    with (testitem_subdir / "testitem.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.DictWriter(
            handle, fieldnames=["molecule_chembl_id", "standard_inchi_skeleton"]
        )
        writer.writeheader()
        writer.writerow(
            {"molecule_chembl_id": "CHEMBL1", "standard_inchi_skeleton": "InChI=1"}
        )
        writer.writerow(
            {"molecule_chembl_id": "CHEMBL2", "standard_inchi_skeleton": "InChI=2"}
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

    with caplog.at_level(logging.INFO):
        output_path = process_activity_extended(
            input_path=activity_path,
            dictionary_dir=dictionary_root,
        )

    assert output_path.exists()

    filled_records = [
        record
        for record in caplog.records
        if "activity_extended_missing_columns_filled" in record.message
    ]
    assert filled_records, "Expected log about backfilled columns not emitted"

    message = filled_records[0].message
    assert filled_records[0].levelno == logging.INFO
    for column_name in ("activity_chembl_id", "compound_key", "log_value"):
        assert column_name in message

    unresolved_messages = [
        record.message
        for record in caplog.records
        if "activity_extended_missing_columns_unresolved" in record.message
    ]
    assert len(unresolved_messages) <= 1
    for message in unresolved_messages:
        assert "nstereo" in message


@pytest.mark.integration
def test_process_activity_extended__warns_on_unresolved_parent(tmp_path, caplog):
    """Warn when parent identifiers remain unresolved after enrichment."""

    _fix_seed()

    activity_dir = tmp_path / "input"
    dictionary_root = tmp_path / "dictionary"
    activity_dir.mkdir()

    _prepare_dictionary(dictionary_root)

    hierarchy_path = dictionary_root / "_testitem" / "molecule_hierarchy.csv"
    hierarchy_path.parent.mkdir(parents=True, exist_ok=True)
    with hierarchy_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=["molecule_chembl_id", "parent_molecule_chembl_id"]
        )
        writer.writeheader()

    activity_path = activity_dir / "output.activities_20240103.csv"
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
        "molecule_pref_name",
        "pchembl_value",
    ]
    with activity_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=activity_columns)
        writer.writeheader()
        writer.writerow(
            {
                "activity_id": "201",
                "molecule_chembl_id": "CHEMBL1",
                "target_chembl_id": "CHEMBLT1",
                "assay_chembl_id": "ASSAY1",
                "document_chembl_id": "DOC1",
                "bao_endpoint": "endpoint",
                "standard_type": "IC50",
                "standard_value": "9.0",
                "bao_format": "format",
                "molecule_pref_name": "Example compound",
                "pchembl_value": "6.5",
            }
        )

    with caplog.at_level(logging.INFO):
        output_path = process_activity_extended(
            input_path=activity_path,
            dictionary_dir=dictionary_root,
        )

    assert output_path.exists()

    unresolved_logs = [
        record
        for record in caplog.records
        if "activity_extended_missing_columns_unresolved" in record.message
    ]
    assert unresolved_logs, "Expected warning about unresolved columns not emitted"
    assert unresolved_logs[0].levelno == logging.WARNING
    assert "parent_molecule_chembl_id" in unresolved_logs[0].message

    filled_logs = [
        record
        for record in caplog.records
        if "activity_extended_missing_columns_filled" in record.message
    ]
    assert filled_logs, "Expected info log about filled columns not emitted"
    assert filled_logs[0].levelno == logging.INFO
    assert "compound_key" in filled_logs[0].message


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

    expected_sorted = output_df.sort_values(by=subset, kind="mergesort").reset_index(
        drop=True
    )
    pd.testing.assert_frame_equal(output_df.reset_index(drop=True), expected_sorted)

    dedupe_records = [
        record
        for record in caplog.records
        if "activity_extended_deduplicated" in record.message
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
