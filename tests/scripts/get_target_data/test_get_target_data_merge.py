"""Tests for merging IUPHAR classifications with existing target data."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
import pytest
from pytest import MonkeyPatch

import library.pipelines.target.postprocessing as tp
from library.pipelines.target import protein_classification as pc
from library.config import Config
from library.cli_utils import PipelineError
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
            "geneName": ["ChemblGeneName"],
            "family": ["chembl_family"],
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
            "geneName": ["UniGeneName"],
            "family": ["uni_family"],
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

    allowed_suffix_columns = {"xref_chembl", "xref_uniprot"}
    disallowed_combined = [
        column
        for column in combined_df.columns
        if (
            column.endswith("_chembl") and column not in allowed_suffix_columns
        )
        or (
            column.endswith("_uniprot") and column not in allowed_suffix_columns
        )
    ]
    assert not disallowed_combined
    assert combined_df.loc[0, "isoform_ids"] == "UNI_ISO"
    assert combined_df.loc[0, "isoform_names"] == "ChemblIso"
    assert combined_df.loc[0, "GuidetoPHARMACOLOGY"] == "uni_gtop"
    assert combined_df.loc[0, "gene"] == "UNIPROTGENE"
    assert combined_df.loc[0, "geneName"] == "UniGeneName"
    assert combined_df.loc[0, "family"] == "uni_family"

    merged = combined_df.merge(iuphar_df, on="uniprot_id", how="left")
    processed = tp.postprocess_targets(merged)

    assert processed.loc[0, "isoform_ids"].lower() == "uni_iso"
    assert processed.loc[0, "isoform_names"].lower() == "chembliso"
    assert processed.loc[0, "GuidetoPHARMACOLOGY"] == "uni_gtop"
    gene_value = processed.loc[0, "gene"]
    assert gene_value.startswith("UNIPROTGENE")
    assert "UNIGENENAME" in gene_value
    assert processed.loc[0, "geneName"] == "UniGeneName"
    assert processed.loc[0, "family"] == "uni_family"
    disallowed_processed = [
        column
        for column in processed.columns
        if (
            column.endswith("_chembl") and column not in allowed_suffix_columns
        )
        or (
            column.endswith("_uniprot") and column not in allowed_suffix_columns
        )
    ]
    assert not disallowed_processed


def test_fetch_iuphar_missing_configured_uniprot_column(
    monkeypatch: MonkeyPatch, tmp_path: Path, cfg: Config
) -> None:
    """Ensure fetch step aborts when the configured column disappears."""

    cfg.target.all.uniprot_column = "unexpected_column"

    chembl_df = pd.DataFrame(
        {
            "target_chembl_id": ["CHEMBL1"],
            "uniprot_id": ["P12345"],
            "unexpected_column": ["P12345"],
        }
    )
    uniprot_df = pd.DataFrame(
        {
            "uniprot_id": ["P12345"],
            "original_id": ["P12345"],
        }
    )
    out = tmp_path / "iuphar.csv"

    def fake_run_iuphar(cfg: Config, args: argparse.Namespace) -> int:
        pd.DataFrame({"uniprot_id": ["P12345"], "IUPHAR_class": ["Enzyme"]}).to_csv(
            args.output_csv,
            index=False,
        )
        return 0

    monkeypatch.setattr(gtd, "run_iuphar", fake_run_iuphar)

    original_merge = pd.merge

    def dropping_merge(*args: object, **kwargs: object) -> pd.DataFrame:
        result = original_merge(*args, **kwargs)
        if "unexpected_column" in result.columns:
            result = result.drop(columns=["unexpected_column"])
        return result

    monkeypatch.setattr(pd, "merge", dropping_merge)

    error_events: list[str] = []
    original_error = gtd.logger.error

    def tracking_error(event: str, *args: object, **kwargs: object) -> object:
        error_events.append(event)
        return original_error(event, *args, **kwargs)

    monkeypatch.setattr(gtd.logger, "error", tracking_error)

    with pytest.raises(PipelineError, match="unexpected_column"):
        gtd.fetch_iuphar(cfg, chembl_df, uniprot_df, out)

    assert "missing_configured_uniprot_column" in error_events
