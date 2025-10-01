"""Tests for merging IUPHAR classifications with existing target data."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
from pytest import MonkeyPatch

import library.target_postprocessing as tp
from library import protein_classification as pc
from library.config import Config
from scripts import get_target_data as gtd


class DummyRecord:
    """Lightweight classification record for deterministic predictions."""

    def __init__(self, status: str = "target_id") -> None:
        self.STATUS = status
        self.IUPHAR_class = "ClassA"
        self.IUPHAR_subclass = "SubclassA"
        self.IUPHAR_type = "TypeA"


class DummyClassifier:
    """Minimal classifier returning deterministic values for tests."""

    def __init__(self) -> None:
        self.calls: list[tuple[str, tuple[str, ...]]] = []

    def get(
        self, target_id: str, family_id: str, ec_number: str, name: str
    ) -> DummyRecord:
        self.calls.append(("get", (target_id, family_id, ec_number, name)))
        return DummyRecord()

    def by_molecular_function(self, molecular_function: str) -> DummyRecord:
        self.calls.append(("by_molecular_function", (molecular_function,)))
        return DummyRecord(status="molecular_function")


def test_iuphar_merge_preserves_ec_number(
    tmp_path: Path, cfg: Config, monkeypatch: MonkeyPatch
) -> None:
    """Ensure merging does not duplicate the ``ec_number`` column."""
    combined_df = pd.DataFrame(
        {
            "target_chembl_id": ["CHEMBL1"],
            "uniprotkb_Id": ["P12345"],
            "uniprot_id": ["P12345"],
            "gene": ["GN1"],
            "gene_name": ["GN1"],
            "synonyms": ["foo"],
            "ec_number": ["1.1.1.1"],
            "genus": ["Homo"],
            "lineage_superkingdom": ["Eukaryota"],
            "lineage_phylum": ["Chordata"],
            "lineage_class": ["Mammalia"],
        }
    )
    iuphar_df = pd.DataFrame(
        {
            "uniprot_id": ["P12345"],
            "IUPHAR_class": ["Enzyme"],
            "target_id": ["T1"],
        }
    )

    monkeypatch.setattr(tp, "postprocess_targets", lambda df: df)
    classifier = DummyClassifier()
    monkeypatch.setattr(pc, "classifier_from_config", lambda cfg: classifier)

    orig_finalise = tp.finalise_targets

    def patched_finalise(df: pd.DataFrame, **kwargs: object) -> pd.DataFrame:
        df = df.drop(
            columns=[col for col in ("type", "target_type") if col in df.columns]
        )
        return orig_finalise(df, **kwargs)

    monkeypatch.setattr(tp, "finalise_targets", patched_finalise)
    merged = gtd.merge_results(combined_df, iuphar_df, cfg)

    assert "protein_class_pred_L1" in merged.columns
    assert merged.loc[0, "protein_class_pred_L1"] == "ClassA"
    assert merged.loc[0, "target_type"] == "Multicellular organism"
    assert classifier.calls


def test_fetch_iuphar_prefers_uniprot_values(
    tmp_path: Path, cfg: Config, monkeypatch: MonkeyPatch
) -> None:
    """Ensure UniProt values override ChEMBL data for overlapping fields."""

    chembl_df = pd.DataFrame(
        {
            "target_chembl_id": ["CHEMBL1"],
            "uniprot_id": ["P11111"],
            "isoform_ids": ["CHEMBL_ISO"],
            "GuidetoPHARMACOLOGY": ["CHEMBL_GUIDE"],
            "geneName": ["CHEMBL_GENE"],
            "gene": ["CHEMBL_GENE"],
            "pref_name": ["Chembl Preferred"],
            "component_description": ["Chembl Component"],
            "names": ["Chembl Name"],
            "synonyms": ["Chembl Synonym"],
        }
    )
    uniprot_df = pd.DataFrame(
        {
            "uniprot_id": ["P11111"],
            "isoform_ids": ["UNIPROT_ISO"],
            "GuidetoPHARMACOLOGY": ["UNIPROT_GUIDE"],
            "geneName": ["UniGene"],
            "gene": ["UniGene"],
            "pref_name": ["UniProt Preferred"],
            "component_description": ["UniProt Component"],
            "names": ["UniProt Name"],
            "synonyms": ["UniProt Synonym"],
            "original_id": ["P11111"],
        }
    )

    output_csv = tmp_path / "iuphar.csv"

    def fake_run_iuphar(config: Config, args: argparse.Namespace) -> int:
        pd.DataFrame(
            {
                "uniprot_id": ["P11111"],
                "IUPHAR_class": ["ClassA"],
            }
        ).to_csv(
            args.output_csv,
            index=False,
            sep=config.io.csv_sep,
            encoding=config.io.csv_encoding,
        )
        return 0

    monkeypatch.setattr(gtd, "run_iuphar", fake_run_iuphar)

    combined_df, iuphar_df = gtd.fetch_iuphar(cfg, chembl_df, uniprot_df, output_csv)

    assert combined_df.loc[0, "isoform_ids"] == "UNIPROT_ISO"
    assert combined_df.loc[0, "GuidetoPHARMACOLOGY"] == "UNIPROT_GUIDE"
    assert combined_df.loc[0, "geneName"] == "UniGene"
    assert combined_df.loc[0, "pref_name"] == "UniProt Preferred"

    conflicting = {"isoform_ids", "GuidetoPHARMACOLOGY", "geneName", "pref_name"}
    for base in conflicting:
        assert f"{base}_chembl" not in combined_df.columns
        assert f"{base}_uniprot" not in combined_df.columns

    merged = combined_df.merge(iuphar_df, on="uniprot_id", how="left")
    processed = tp.postprocess_targets(merged)

    assert processed.loc[0, "isoform_ids"] == "UNIPROT_ISO"
    assert processed.loc[0, "GuidetoPHARMACOLOGY"] == "UNIPROT_GUIDE"
    assert processed.loc[0, "geneName"] == "UniGene"

    for base in conflicting:
        assert f"{base}_chembl" not in processed.columns
        assert f"{base}_uniprot" not in processed.columns
