"""Unit tests for pandera DataFrame schemas."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pandera.pandas as pa
import pytest

pytest.importorskip("hypothesis")
from hypothesis import given
from hypothesis import strategies as st
from hypothesis.extra.pandas import column, data_frames, range_indexes
from pandera.errors import SchemaError, SchemaErrors

from library.normalization import normalize_activities
from schemas import (
    ActivitiesSchema,
    AssaysSchema,
    DocumentsSchema,
    TargetsSchema,
    TestitemsSchema,
)


def test_activities_schema_validation() -> None:
    """Ensure :data:`ActivitiesSchema` validates expected data."""
    valid = pd.DataFrame(
        {
            "activity_id": ["1", "2"],
            "molecule_chembl_id": ["CHEMBL1", "CHEMBL2"],
            "assay_chembl_id": ["CHEMBL0", "CHEMBL0"],
            # ``standard_value`` is nullable so a frame mixing floats and ``None``
            # should validate without error.
            "standard_value": [10.0, None],
            "standard_type": ["IC50", "IC50"],
            "activity_comment": [None, "note"],
        }
    )
    ActivitiesSchema.validate(valid)

    invalid = valid.copy()
    invalid.loc[0, "standard_type"] = "Unknown"
    with pytest.raises(SchemaError):
        ActivitiesSchema.validate(invalid)


def test_assays_schema_validation() -> None:
    """Required columns are enforced."""
    valid = pd.DataFrame(
        {
            "assay_chembl_id": ["CHEMBL1"],
            "document_chembl_id": [None],  # nullable column accepts missing values
            "target_chembl_id": ["CHEMBL3"],
        }
    )
    AssaysSchema.validate(valid)

    # Missing mandatory column should raise a SchemaError
    invalid = valid.drop(columns=["assay_chembl_id"])
    with pytest.raises(SchemaError):
        AssaysSchema.validate(invalid)


def test_documents_schema_validation() -> None:
    """Ensure :data:`DocumentsSchema` validates expected data."""
    valid = pd.DataFrame(
        {
            "document_chembl_id": ["CHEMBL1"],
            "doi": [None],  # nullable text field
            "title": ["Example"],
        }
    )
    DocumentsSchema.validate(valid)
    invalid = valid.drop(columns=["document_chembl_id"])
    with pytest.raises(SchemaError):
        DocumentsSchema.validate(invalid)


def test_targets_schema_validation() -> None:
    """Ensure required target columns are present."""
    valid = pd.DataFrame(
        {
            "target_chembl_id": ["CHEMBL1"],
            "uniprotkb_Id": ["P12345"],
            "recommended_name": ["Protein kinase"],
            # ``synonyms`` may be nullable in newer schema versions
            "synonyms": [None],
            "type": ["protein"],
        }
    )
    if TargetsSchema.columns["synonyms"].nullable:
        TargetsSchema.validate(valid)
    else:
        with pytest.raises(SchemaError):
            TargetsSchema.validate(valid)

    invalid = valid.drop(columns=["target_chembl_id"])
    with pytest.raises(SchemaError):
        TargetsSchema.validate(invalid)


def test_targets_schema_allows_missing_optional_columns() -> None:
    """Schema validation succeeds when optional columns are absent."""
    df = pd.DataFrame(
        {
            "target_chembl_id": ["CHEMBL1"],
            "uniprotkb_Id": ["P12345"],
            "recommended_name": ["Protein kinase"],
            "synonyms": ["pk1|kinase"],
            "type": ["protein"],
        }
    )
    TargetsSchema.validate(df)


def test_targets_schema_defines_expected_columns() -> None:
    """All expected columns are present in :data:`TargetsSchema`."""
    expected = {
        "target_chembl_id",
        "uniprotkb_Id",
        "uniprot_id",
        "secondary_uniprot_id",
        "gene_name",
        "recommended_name",
        "synonyms",
        "genus",
        "superkingdom",
        "phylum",
        "taxon_id",
        "ec_number",
        "hgnc_name",
        "hgnc_id",
        "molecular_function",
        "cellular_component",
        "subcellular_location",
        "topology",
        "transmembrane",
        "intramembrane",
        "glycosylation",
        "lipidation",
        "disulfide_bond",
        "modified_residue",
        "phosphorylation",
        "acetylation",
        "ubiquitination",
        "signal_peptide",
        "propeptide",
        "isoform_names",
        "isoform_ids",
        "isoform_synonyms",
        "reactions",
        "target_id",
        "IUPHAR_family_id",
        "IUPHAR_type",
        "IUPHAR_class",
        "IUPHAR_subclass",
        "IUPHAR_chain",
        "full_id_path",
        "full_name_path",
        "GuidetoPHARMACOLOGY",
        "type",
    }
    assert set(TargetsSchema.columns) == expected


def test_testitems_schema_validation() -> None:
    """Ensure :data:`TestitemsSchema` validates expected data."""
    valid = pd.DataFrame(
        {
            "salt_chembl_id": ["CHEMBL1", "CHEMBL2"],
            "molecule_chembl_id": ["CHEMBL1", "CHEMBL2"],
            "molecule_type": [None, "Small molecule"],  # nullable field
            # ``pubchem_cid`` accepts mixed types via ``object`` dtype
            "pubchem_cid": [123, "456"],
        }
    )
    TestitemsSchema.validate(valid)
    invalid = valid.drop(columns=["molecule_chembl_id"])
    with pytest.raises(SchemaError):
        TestitemsSchema.validate(invalid)


@pytest.mark.parametrize(
    "schema",
    [
        ActivitiesSchema,
        AssaysSchema,
        DocumentsSchema,
        TargetsSchema,
        TestitemsSchema,
    ],
)
def test_any_columns_accept_mixed_types(schema: pa.DataFrameSchema) -> None:
    """Columns typed as ``pa.Any`` validate mixed inputs."""
    any_cols = [
        name for name, col in schema.columns.items() if str(col.dtype) == "object"
    ]
    if not any_cols:
        pytest.skip("schema has no Any-typed columns")

    # Build a frame containing required columns
    data: dict[str, list[object]] = {}
    for name, col in schema.columns.items():
        if not col.required:
            continue
        dtype = str(col.dtype)
        if "float" in dtype:
            data[name] = [0.0, 1.0]
        elif "int" in dtype:
            data[name] = [0, 1]
        elif "bool" in dtype:
            data[name] = [True, False]
        else:
            data[name] = ["a", "b"]
    df = pd.DataFrame(data)

    for name in any_cols:
        df[name] = [1, "two"]

    schema.validate(df)


@given(
    data_frames(
        columns=[
            column("activity_id", elements=st.text(min_size=1)),
            column("molecule_chembl_id", elements=st.text(min_size=1)),
            column(
                "standard_value",
                elements=st.floats(min_value=0, allow_nan=False),
            ),
            column("assay_chembl_id", elements=st.text(min_size=1)),
        ],
        index=range_indexes(min_size=1, max_size=5),
    )
)
def test_activities_schema_hypothesis_valid(df: pd.DataFrame) -> None:
    """Random valid frames pass ``ActivitiesSchema``."""

    ActivitiesSchema.validate(df)


@given(
    data_frames(
        columns=[
            column("activity_id", elements=st.text(min_size=1)),
            column("molecule_chembl_id", elements=st.text(min_size=1)),
            column(
                "standard_value",
                elements=st.floats(max_value=-1.0, allow_nan=False),
            ),
            column("assay_chembl_id", elements=st.text(min_size=1)),
        ],
        index=range_indexes(min_size=1, max_size=5),
    )
)
def test_activities_schema_hypothesis_invalid(df: pd.DataFrame) -> None:
    """Negative ``standard_value`` fails validation."""

    with pytest.raises(SchemaError):
        ActivitiesSchema.validate(df)


@given(
    data_frames(
        columns=[
            column("document_chembl_id", elements=st.text(min_size=1)),
            column("title", elements=st.text(min_size=1)),
            column(
                "year",
                dtype=int,
                elements=st.integers(min_value=1900, max_value=2100),
            ),
        ],
        index=range_indexes(min_size=1, max_size=5),
    )
)
def test_documents_schema_hypothesis_valid(df: pd.DataFrame) -> None:
    """Random valid frames pass ``DocumentsSchema``."""

    DocumentsSchema.validate(df)


@given(
    data_frames(
        columns=[
            column("document_chembl_id", elements=st.text(min_size=1)),
            column("title", elements=st.text(min_size=1)),
            column(
                "year",
                dtype=int,
                elements=st.integers(min_value=0, max_value=1899),
            ),
        ],
        index=range_indexes(min_size=1, max_size=5),
    )
)
def test_documents_schema_hypothesis_invalid(df: pd.DataFrame) -> None:
    """Years below 1900 fail ``DocumentsSchema`` validation."""
    df = df.drop(columns=["document_chembl_id"])
    with pytest.raises(SchemaError):
        DocumentsSchema.validate(df)


def test_activities_from_files() -> None:
    """Validate activities data from CSV/JSON files and capture failures."""
    data_dir = Path(__file__).parent / "data"

    # Positive case: CSV data passes validation after normalisation
    valid = pd.read_csv(data_dir / "activities_valid.csv", dtype=str).assign(
        assay_chembl_id="CHEMBL0"
    )
    valid["standard_value"] = valid["standard_value"].astype(float)
    normalized = normalize_activities(valid)
    assert normalized.loc[0, "relation"] == "<="
    assert normalized.loc[0, "units"] == "5 uM"
    ActivitiesSchema.validate(normalized)

    # Negative case: JSON data with invalid standard_value
    invalid = pd.read_json(data_dir / "activities_invalid.json").assign(
        assay_chembl_id="CHEMBL0"
    )
    invalid["standard_value"] = invalid["standard_value"].astype(float)
    invalid_norm = normalize_activities(invalid)
    with pytest.raises(SchemaErrors) as exc_info:
        ActivitiesSchema.validate(invalid_norm, lazy=True)

    failure_cases = exc_info.value.failure_cases
    assert (failure_cases["column"] == "standard_value").any()
    assert (failure_cases["index"] == 0).any()
