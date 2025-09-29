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
from library.constants import (
    ActivitiesSchema,
    AssaysSchema,
    DocumentsSchema,
    TargetsSchema,
    TestitemsSchema,
)
from library.constants import TARGETS_COLUMN_ORDER


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


def test_activities_schema_accepts_configured_standard_types() -> None:
    """Configured standard types are permitted by the schema."""

    df = pd.DataFrame(
        {
            "activity_id": ["1", "2", "3"],
            "molecule_chembl_id": ["CHEMBL1", "CHEMBL2", "CHEMBL3"],
            "assay_chembl_id": ["CHEMBL0", "CHEMBL0", "CHEMBL1"],
            "standard_value": [1.0, 2.0, 3.0],
            "standard_type": ["IC50", "EC50", "KD"],
        }
    )

    ActivitiesSchema.validate(df)


def test_activities_schema_accepts_object_dtypes() -> None:
    """``standard_value`` validates with object dtype."""

    df = pd.DataFrame(
        {
            # numeric IDs are coerced to ``object`` to allow flexible typing
            "activity_id": pd.Series([1], dtype=object),
            "molecule_chembl_id": ["CHEMBL1"],
            "assay_chembl_id": ["CHEMBL0"],
            "standard_value": pd.Series([1.0], dtype=object),
            "pchembl_value": pd.Series([5.0], dtype=object),
            "src_assay_id": pd.Series([2], dtype=object),
            "src_id": pd.Series([3], dtype=object),
            "value": pd.Series([7.0], dtype=object),
        }
    )
    ActivitiesSchema.validate(df)


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


def test_assays_schema_allows_varied_dtypes_and_nulls() -> None:
    """Int and bool fields accept flexible representations."""
    df = pd.DataFrame(
        {
            "assay_chembl_id": [None],
            "acts_per_assay_step5": ["10"],
            "cited_assay_corr": [1],
            "month": ["5"],
            "shuffled_cit": ["false"],
            "version": [None],
            "year": ["2024"],
        }
    )
    AssaysSchema.validate(df)


def test_documents_schema_validation() -> None:
    """Ensure :data:`DocumentsSchema` validates expected data."""
    valid = pd.DataFrame(
        {
            "document_chembl_id": ["CHEMBL1"],
            "title": ["Example"],
            "PubMed.PMID": [12345],
            "OpenAlex.Error": [None],
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
            "uniprot_id_primary": ["P12345"],
        }
    )
    TargetsSchema.validate(valid)

    invalid = valid.drop(columns=["target_chembl_id"])
    with pytest.raises(SchemaError):
        TargetsSchema.validate(invalid)


def test_targets_schema_allows_missing_optional_columns() -> None:
    """Schema validation succeeds when optional columns are absent."""
    df = pd.DataFrame({"target_chembl_id": ["CHEMBL1"]})
    TargetsSchema.validate(df)


def test_targets_schema_defines_expected_columns() -> None:
    """All expected columns are present in :data:`TargetsSchema`."""
    assert set(TargetsSchema.columns) == set(TARGETS_COLUMN_ORDER)


def test_targets_schema_nullable_and_any_columns() -> None:
    """All columns are nullable and bool/int types use ``object`` dtype."""

    # Every column should be nullable
    for col in TargetsSchema.columns.values():
        assert col.nullable is True

    # Columns previously typed as ``bool`` or ``int`` now accept any object
    any_columns = [
        "taxon_id",
        "features_signal_peptide",
        "features_transmembrane",
        "ptm_glycosylation",
        "ptm_lipidation",
        "ptm_disulfide_bond",
        "ptm_modified_residue",
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
    ]
    for name in any_columns:
        assert str(TargetsSchema.columns[name].dtype) == "object"

    # Schema validation with null values in these columns should succeed
    df = pd.DataFrame(
        {"target_chembl_id": ["CHEMBL1"], **{c: [None] for c in any_columns}}
    )
    TargetsSchema.validate(df)


def test_testitems_schema_validation() -> None:
    """Ensure :data:`TestitemsSchema` validates expected data."""
    valid = pd.DataFrame(
        {
            "molecule_chembl_id": ["CHEMBL1"],
            "parent_molecule_id": ["CHEMBL0"],
            "first_approval": [1950],
            "black_box_warning": [0],
            "oral": [True],
            "parenteral": [False],
            "topical": [False],
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
                "PubMed.PMID",
                dtype=int,
                elements=st.integers(min_value=0, max_value=2**63 - 1),
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
            column(
                "PubMed.PMID",
                dtype=int,
                elements=st.integers(min_value=0, max_value=2**63 - 1),
            ),
        ],
        index=range_indexes(min_size=1, max_size=5),
    )
)
def test_documents_schema_hypothesis_invalid(df: pd.DataFrame) -> None:
    """Missing required columns fail ``DocumentsSchema`` validation."""
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
