import pytest
from scripts import get_document_data

from library.pipelines.document import pipeline as document_pipeline
from library.schemas import DocumentsSchema
from library.schemas.document_spec import (
    DOCUMENT_COLUMN_GROUPS,
    DOCUMENT_COLUMN_SPECS,
    DOCUMENT_EXPORT_COLUMNS,
    DOCUMENT_SCHEMA_COLUMNS,
)


@pytest.mark.unit
def test_document_schema_groups__pipeline_alignment() -> None:
    assert document_pipeline.CH_EMBL_COLUMNS == list(DOCUMENT_COLUMN_GROUPS["chembl"])
    assert document_pipeline.PUBMED_COLUMNS == list(DOCUMENT_COLUMN_GROUPS["pubmed"])
    assert document_pipeline.SEMANTIC_SCHOLAR_COLUMNS == list(
        DOCUMENT_COLUMN_GROUPS["scholar"]
    )
    assert document_pipeline.OPENALEX_COLUMNS == list(
        DOCUMENT_COLUMN_GROUPS["openalex"]
    )
    assert document_pipeline.CROSSREF_COLUMNS == list(
        DOCUMENT_COLUMN_GROUPS["crossref"]
    )
    assert document_pipeline.DOCUMENT_SCHEMA_COLUMNS == list(DOCUMENT_SCHEMA_COLUMNS)


@pytest.mark.unit
def test_document_schema_export__shared_constant() -> None:
    assert get_document_data._EXPORT_COLUMNS == list(DOCUMENT_EXPORT_COLUMNS)


@pytest.mark.unit
@pytest.mark.parametrize(
    "column_name, expected_dtype, required, nullable, coerce",
    [
        ("document_chembl_id", str, True, True, False),
        ("PubMed.PMID", object, False, True, True),
        ("publication_type_score_review", int, False, True, True),
    ],
)
def test_document_schema_columns__attributes(
    column_name: str,
    expected_dtype: type,
    required: bool,
    nullable: bool,
    coerce: bool,
) -> None:
    spec = DOCUMENT_COLUMN_SPECS[column_name]
    assert spec.name == column_name
    assert spec.required is required
    assert spec.nullable is nullable
    assert spec.coerce is coerce
    column = DocumentsSchema.columns[column_name]
    assert str(column.dtype).lower().startswith(expected_dtype.__name__)
    assert column.required is required
    assert column.nullable is nullable
    assert getattr(column, "coerce", False) is coerce
