"""Tests for :mod:`library.target_postprocessing`.

The suite covers :func:`library.target_postprocessing.postprocess_targets` and
its convenience wrapper :func:`library.target_postprocessing.postprocess_file`.
"""

from __future__ import annotations

import warnings
from pathlib import Path

import pandas as pd

from library import target_postprocessing as tp
from library.config import Config, IoCfg
from schemas import TargetsSchema


def test_postprocess_targets_merges_and_normalises() -> None:
    """``postprocess_targets`` merges fields and normalises tokens."""

    df = pd.DataFrame(
        {
            "chembl_id": ["CHEMBL1"],
            "uniProtkbId": ["P12345-2_SOMETHING"],
            "uniprot_id": ["P12345"],
            "secondaryAccessions": ["S12345|S67890"],
            "pref_name": ["Protein kinase"],
            "gene_name_x": ["g1"],
            "geneName": [""],
            "gene": ["g2|g3"],
            "secondaryAccessionNames": ["name2|name3"],
            "component_description": ["desc"],
            "chembl_alternative_name": ["alt"],
            "recommendedName": ["Rec Name"],
            "names": ["name1|name4"],
            "synonyms_x": ["synX"],
            "synonyms": ["syn1|syn2"],
            "ec_number": ["1.1.1.1"],
            "ec_code": ["2.2.2.2"],
        }
    )
    out = tp.postprocess_targets(df, chembl_col="chembl_id")

    row = out.iloc[0]
    assert row["uniprotkb_Id"] == "P12345"
    assert row["secondary_uniprot_id"] == "S12345|S67890"
    assert row["recommended_name"] == "Protein kinase"
    assert row["gene_name"] == "G1"
    assert row["ec_number"] == "1.1.1.1|2.2.2.2"
    assert row["isoform_names"] == "-"
    assert row["synonyms"].startswith("g2|g3|g1|name2|name3")


def test_postprocess_targets_orders_columns() -> None:
    """Columns from the schema appear first and others alphabetically."""

    df = pd.DataFrame(
        {
            "chembl_id": ["CHEMBL1"],
            "uniprot_id": ["P1"],
            "gene": ["g1"],
            "b": ["1"],
            "a": ["2"],
        }
    )
    out = tp.postprocess_targets(df, chembl_col="chembl_id")

    schema_cols = list(TargetsSchema.columns)
    schema_cols[schema_cols.index("target_chembl_id")] = "chembl_id"
    extra_cols = sorted(c for c in out.columns if c not in schema_cols)
    assert list(out.columns) == schema_cols + extra_cols


def test_postprocess_file_roundtrip(tmp_path: Path) -> None:
    """``postprocess_file`` respects ``IoCfg`` defaults for CSV options."""

    df = pd.DataFrame(
        {
            "chembl_id": ["CHEMBL1"],
            "uniProtkbId": ["P12345-2_SOMETHING"],
            "uniprot_id": ["P12345"],
            "secondaryAccessions": ["S12345|S67890"],
            "pref_name": ["Protein kinase"],
            "gene_name_x": ["g1"],
            "geneName": [""],
            "gene": ["g2|g3"],
            "secondaryAccessionNames": ["name2|name3"],
            "component_description": ["desc"],
            "chembl_alternative_name": ["alt"],
            "recommendedName": ["Rec Name"],
            "names": ["name1|name4"],
            "synonyms_x": ["synX"],
            "synonyms": ["syn1|syn2"],
            "ec_number": ["1.1.1.1"],
            "ec_code": ["2.2.2.2"],
        }
    )

    cfg = IoCfg(csv_sep=";", csv_encoding="utf8")
    input_path = tmp_path / "in.csv"
    df.to_csv(input_path, index=False, sep=cfg.csv_sep, encoding=cfg.csv_encoding)
    output_path = tmp_path / "out.csv"

    tp.postprocess_file(input_path, output_path, cfg=cfg, chembl_col="chembl_id")

    expected = tp.postprocess_targets(df, chembl_col="chembl_id").astype(str)
    result = pd.read_csv(
        output_path,
        dtype=str,
        sep=cfg.csv_sep,
        encoding=cfg.csv_encoding,
        keep_default_na=False,
    )
    pd.testing.assert_frame_equal(result, expected)


def test_finalise_targets_filters_duplicates_and_merges() -> None:
    """``finalise_targets`` drops invalid rows and merges organism data."""

    df = pd.DataFrame(
        {
            "chembl_id": ["CHEMBL1", "CHEMBL1", "CHEMBL2"],
            "uniprot": ["P12345", "P12345", "nan"],
            "organism": ["Homo", "Homo", "Mus"],
            "isoform_names": ["IsoA", "IsoB", "IsoC"],
            "synonyms": ["SynA", "SynB", "SynC"],
            "SUPFAM": ["s1", "s2", "s3"],
            "transmembrane": ["True", "True", "False"],
        }
    )
    organism = pd.DataFrame({"organism": ["Homo"], "type": ["Mammal"]})

    out = tp.finalise_targets(
        df,
        organism,
        chembl_col="chembl_id",
        uniprot_col="uniprot",
        genus_col="organism",
    )

    assert list(out["chembl_id"]) == ["CHEMBL1"]
    assert "SUPFAM" not in out.columns
    assert out.loc[0, "isoform_names"] == "isoa"
    assert out.loc[0, "type"] == "Mammal"
    assert out.loc[0, "organism"] == "Homo"
    assert out["transmembrane"].dtype == "boolean"


def test_finalise_targets_orders_columns() -> None:
    """Columns from the schema appear first and others alphabetically."""

    df = pd.DataFrame(
        {
            "uniprotkb_Id": ["P1"],
            "a": ["1"],
            "target_chembl_id": ["CHEMBL1"],
            "b": ["2"],
            "genus": ["Homo"],
        }
    )
    organism = pd.DataFrame({"genus": ["Homo"], "type": ["Mammal"]})

    out = tp.finalise_targets(df, organism)

    schema_cols = list(TargetsSchema.columns)
    extra_cols = sorted(c for c in out.columns if c not in schema_cols)
    assert list(out.columns) == schema_cols + extra_cols


def test_finalise_file_roundtrip(tmp_path: Path, cfg: Config) -> None:
    """``finalise_file`` reads, processes and writes the expected table."""

    df = pd.DataFrame(
        {
            "chembl_id": ["CHEMBL1", "CHEMBL1", "CHEMBL2"],
            "uniprot": ["P12345", "P12345", "nan"],
            "organism": ["Homo", "Homo", "Mus"],
            "synonyms": ["SynA", "SynB", "SynC"],
            "SUPFAM": ["s1", "s2", "s3"],
        }
    )
    organism = pd.DataFrame({"organism": ["Homo"], "type": ["Mammal"]})

    input_path = tmp_path / "in.csv"
    df.to_csv(input_path, index=False)
    organism_path = tmp_path / "org.csv"
    organism.to_csv(organism_path, index=False)
    output_path = tmp_path / "out.csv"

    cfg.resources.organism_csv = organism_path
    tp.finalise_file(
        input_path,
        output_path,
        cfg=cfg,
        chembl_col="chembl_id",
        uniprot_col="uniprot",
        genus_col="organism",
    )

    expected = (
        tp.finalise_targets(
            df,
            organism,
            chembl_col="chembl_id",
            uniprot_col="uniprot",
            genus_col="organism",
        )
        .fillna("")
        .astype(str)
    )
    result = pd.read_csv(output_path, dtype=str, keep_default_na=False)
    pd.testing.assert_frame_equal(result, expected)


def test_finalise_targets_no_downcast_warning() -> None:
    """``finalise_targets`` should not emit downcast warnings during replace."""

    df = pd.DataFrame(
        {
            "chembl_id": ["CHEMBL1"],
            "uniprot": ["P12345"],
            "organism": ["Homo"],
            "transmembrane": ["True"],
        }
    )
    organism = pd.DataFrame({"organism": ["Homo"], "type": ["Mammal"]})

    with warnings.catch_warnings():
        warnings.simplefilter("error", FutureWarning)
        tp.finalise_targets(
            df,
            organism,
            chembl_col="chembl_id",
            uniprot_col="uniprot",
            genus_col="organism",
        )


def test_finalise_targets_uses_target_chembl_id_by_default() -> None:
    """Default column name ``target_chembl_id`` is preserved after finalisation."""

    df = pd.DataFrame(
        {
            "target_chembl_id": ["CHEMBL1"],
            "uniprotkb_Id": ["P12345"],
            "genus": ["Homo"],
        }
    )
    organism = pd.DataFrame({"genus": ["Homo"], "type": ["Mammal"]})

    out = tp.finalise_targets(df, organism)

    assert "target_chembl_id" in out.columns
    assert "chembl_id" not in out.columns
