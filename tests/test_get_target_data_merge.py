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


def test_fetch_iuphar_prioritises_uniprot_columns(
    monkeypatch: MonkeyPatch, tmp_path: Path, cfg: Config
) -> None:
    """Ensure overlapping columns prefer UniProt data and remain unsuffixed."""

    chembl_df = pd.DataFrame(
        {
            "target_chembl_id": ["C1"],
            "uniprot_id": ["P1"],
            "isoform_ids": ["chembl_iso"],
            "isoform_names": ["ChemblIso"],
            "GuidetoPHARMACOLOGY": ["chembl_gtop"],
            "gene": ["CHEMBLGENE"],
        }
    )
    uniprot_df = pd.DataFrame(
        {
            "uniprot_id": ["P1"],
            "original_id": ["P1"],
            "isoform_ids": ["UNI_ISO"],
            "isoform_names": [""],
            "GuidetoPHARMACOLOGY": ["uni_gtop"],
            "gene": ["UNIPROTGENE"],
        }
    )
    out = tmp_path / "iuphar.csv"

    def fake_run_iuphar(cfg: Config, args: argparse.Namespace) -> int:
        pd.DataFrame({"uniprot_id": ["P1"], "IUPHAR_class": ["Enzyme"]}).to_csv(
            args.output_csv, index=False
        )
        return 0

    monkeypatch.setattr(gtd, "run_iuphar", fake_run_iuphar)

    combined_df, iuphar_df = gtd.fetch_iuphar(cfg, chembl_df, uniprot_df, out)

    assert all(
        not column.endswith("_chembl") and not column.endswith("_uniprot")
        for column in combined_df.columns
    )
    assert combined_df.loc[0, "isoform_ids"] == "UNI_ISO"
    assert combined_df.loc[0, "isoform_names"] == "ChemblIso"
    assert combined_df.loc[0, "GuidetoPHARMACOLOGY"] == "uni_gtop"
    assert combined_df.loc[0, "gene"] == "UNIPROTGENE"

    merged = combined_df.merge(iuphar_df, on="uniprot_id", how="left")
    processed = tp.postprocess_targets(merged)

    assert processed.loc[0, "isoform_ids"].lower() == "uni_iso"
    assert processed.loc[0, "isoform_names"].lower() == "chembliso"
    assert processed.loc[0, "GuidetoPHARMACOLOGY"] == "uni_gtop"
    assert processed.loc[0, "gene"] == "UNIPROTGENE"
