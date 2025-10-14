"""Declarative Pandera schema for document metadata."""

from __future__ import annotations

from collections import OrderedDict
from typing import Any, ClassVar, Mapping

import pandas as pd

from library._compat.pandera import DataFrameSchema, pa

from .common import int_column, object_column, string_column

__all__ = [
    "DocumentSchema",
    "DOCUMENT_SCHEMA",
    "DOCUMENT_SCHEMA_COLUMNS",
    "DOCUMENT_COLUMN_GROUPS",
    "DOCUMENT_EXPORT_COLUMNS",
    "validate_document_frame",
]


_DOCUMENT_COLUMN_DEFINITIONS = (
    # chembl
    ("document_chembl_id", string_column(required=True, nullable=True, coerce=False)),
    ("title", string_column(required=False, nullable=True, coerce=False)),
    ("abstract", string_column(required=False, nullable=True, coerce=False)),
    ("doi", string_column(required=False, nullable=True, coerce=False)),
    ("year", object_column(dtype=object, required=False, nullable=True, coerce=True)),
    ("journal", string_column(required=False, nullable=True, coerce=False)),
    ("journal_abbrev", string_column(required=False, nullable=True, coerce=False)),
    ("volume", object_column(dtype=object, required=False, nullable=True, coerce=True)),
    ("issue", object_column(dtype=object, required=False, nullable=True, coerce=True)),
    ("first_page", object_column(dtype=object, required=False, nullable=True, coerce=True)),
    ("last_page", object_column(dtype=object, required=False, nullable=True, coerce=True)),
    ("pubmed_id", object_column(dtype=object, required=False, nullable=True, coerce=True)),
    ("authors", string_column(required=False, nullable=True, coerce=False)),
    ("source", string_column(required=False, nullable=True, coerce=False)),
    ("doc_type", string_column(required=False, nullable=True, coerce=False)),
    # derived
    ("doi_normalised", string_column(required=False, nullable=True, coerce=False)),
    (
        "publication_types_normalised",
        string_column(required=False, nullable=True, coerce=False),
    ),
    (
        "publication_type_score_review",
        int_column(required=False, nullable=True, coerce=True),
    ),
    (
        "publication_type_score_experimental",
        int_column(required=False, nullable=True, coerce=True),
    ),
    (
        "publication_type_score_unknown",
        int_column(required=False, nullable=True, coerce=True),
    ),
    ("publication_class", string_column(required=False, nullable=True, coerce=False)),
    # pipeline_status
    ("fetch_status", string_column(required=False, nullable=True, coerce=False)),
    ("error_source", string_column(required=False, nullable=True, coerce=False)),
    # pubmed
    ("PubMed.PMID", object_column(dtype=object, required=False, nullable=True, coerce=True)),
    ("PubMed.DOI", object_column(dtype=object, required=False, nullable=True, coerce=True)),
    ("PubMed.ArticleTitle", string_column(required=False, nullable=True, coerce=False)),
    ("PubMed.Abstract", string_column(required=False, nullable=True, coerce=False)),
    ("PubMed.JournalTitle", string_column(required=False, nullable=True, coerce=False)),
    ("PubMed.Volume", object_column(dtype=object, required=False, nullable=True, coerce=True)),
    ("PubMed.Issue", object_column(dtype=object, required=False, nullable=True, coerce=True)),
    (
        "PubMed.StartPage",
        object_column(dtype=object, required=False, nullable=True, coerce=True),
    ),
    (
        "PubMed.EndPage",
        object_column(dtype=object, required=False, nullable=True, coerce=True),
    ),
    ("PubMed.PublicationType", string_column(required=False, nullable=True, coerce=False)),
    (
        "PubMed.MeSH_Descriptors",
        string_column(required=False, nullable=True, coerce=False),
    ),
    (
        "PubMed.MeSH_Qualifiers",
        string_column(required=False, nullable=True, coerce=False),
    ),
    ("PubMed.ChemicalList", string_column(required=False, nullable=True, coerce=False)),
    (
        "PubMed.DayRevised",
        object_column(dtype=object, required=False, nullable=True, coerce=True),
    ),
    (
        "PubMed.MonthRevised",
        object_column(dtype=object, required=False, nullable=True, coerce=True),
    ),
    (
        "PubMed.YearRevised",
        object_column(dtype=object, required=False, nullable=True, coerce=True),
    ),
    (
        "PubMed.YearCompleted",
        object_column(dtype=object, required=False, nullable=True, coerce=True),
    ),
    (
        "PubMed.MonthCompleted",
        object_column(dtype=object, required=False, nullable=True, coerce=True),
    ),
    (
        "PubMed.DayCompleted",
        object_column(dtype=object, required=False, nullable=True, coerce=True),
    ),
    ("PubMed.Error", object_column(dtype=object, required=False, nullable=True, coerce=True)),
    ("PubMed.ISSN", string_column(required=False, nullable=True, coerce=False)),
    # scholar
    ("scholar.PMID", object_column(dtype=object, required=False, nullable=True, coerce=True)),
    ("scholar.Venue", object_column(dtype=object, required=False, nullable=True, coerce=True)),
    (
        "scholar.PublicationTypes",
        string_column(required=False, nullable=True, coerce=False),
    ),
    (
        "scholar.SemanticScholarId",
        object_column(dtype=object, required=False, nullable=True, coerce=True),
    ),
    (
        "scholar.ExternalIds",
        object_column(dtype=object, required=False, nullable=True, coerce=True),
    ),
    ("scholar.DOI", object_column(dtype=object, required=False, nullable=True, coerce=True)),
    ("scholar.Error", string_column(required=False, nullable=True, coerce=False)),
    # openalex
    (
        "OpenAlex.PublicationTypes",
        object_column(dtype=object, required=False, nullable=True, coerce=True),
    ),
    ("OpenAlex.TypeCrossref", string_column(required=False, nullable=True, coerce=False)),
    ("OpenAlex.Genre", object_column(dtype=object, required=False, nullable=True, coerce=True)),
    ("OpenAlex.Id", string_column(required=False, nullable=True, coerce=False)),
    ("OpenAlex.Venue", object_column(dtype=object, required=False, nullable=True, coerce=True)),
    (
        "OpenAlex.MeshDescriptors",
        string_column(required=False, nullable=True, coerce=False),
    ),
    (
        "OpenAlex.MeshQualifiers",
        object_column(dtype=object, required=False, nullable=True, coerce=True),
    ),
    ("OpenAlex.Error", object_column(dtype=object, required=False, nullable=True, coerce=True)),
    # crossref
    ("crossref.Type", object_column(dtype=object, required=False, nullable=True, coerce=True)),
    ("crossref.Subtype", object_column(dtype=object, required=False, nullable=True, coerce=True)),
    ("crossref.Title", object_column(dtype=object, required=False, nullable=True, coerce=True)),
    ("crossref.Subtitle", object_column(dtype=object, required=False, nullable=True, coerce=True)),
    ("crossref.Subject", object_column(dtype=object, required=False, nullable=True, coerce=True)),
    ("crossref.Error", string_column(required=False, nullable=True, coerce=False)),
    # pipeline_runtime
    ("date_code", string_column(required=False, nullable=True, coerce=False)),
    ("Index", object_column(dtype=object, required=False, nullable=True, coerce=True)),
    ("PubMed.is_review", object_column(dtype=object, required=False, nullable=True, coerce=True)),
    ("scholar.is_review", object_column(dtype=object, required=False, nullable=True, coerce=True)),
    ("OpenAlex.is_review", object_column(dtype=object, required=False, nullable=True, coerce=True)),
    ("pipeline_version", string_column(required=False, nullable=True, coerce=False)),
    ("timestamp_utc", string_column(required=False, nullable=True, coerce=False)),
)


_DOCUMENT_COLUMNS: "OrderedDict[str, pa.Column]" = OrderedDict(_DOCUMENT_COLUMN_DEFINITIONS)

_DOCUMENT_COLUMN_GROUPS: dict[str, tuple[str, ...]] = {
    "chembl": (
        "document_chembl_id",
        "title",
        "abstract",
        "doi",
        "year",
        "journal",
        "journal_abbrev",
        "volume",
        "issue",
        "first_page",
        "last_page",
        "pubmed_id",
        "authors",
        "source",
        "doc_type",
    ),
    "derived": (
        "doi_normalised",
        "publication_types_normalised",
        "publication_type_score_review",
        "publication_type_score_experimental",
        "publication_type_score_unknown",
        "publication_class",
    ),
    "pipeline_status": (
        "fetch_status",
        "error_source",
    ),
    "pubmed": (
        "PubMed.PMID",
        "PubMed.DOI",
        "PubMed.ArticleTitle",
        "PubMed.Abstract",
        "PubMed.JournalTitle",
        "PubMed.Volume",
        "PubMed.Issue",
        "PubMed.StartPage",
        "PubMed.EndPage",
        "PubMed.PublicationType",
        "PubMed.MeSH_Descriptors",
        "PubMed.MeSH_Qualifiers",
        "PubMed.ChemicalList",
        "PubMed.DayRevised",
        "PubMed.MonthRevised",
        "PubMed.YearRevised",
        "PubMed.YearCompleted",
        "PubMed.MonthCompleted",
        "PubMed.DayCompleted",
        "PubMed.Error",
        "PubMed.ISSN",
    ),
    "scholar": (
        "scholar.PMID",
        "scholar.Venue",
        "scholar.PublicationTypes",
        "scholar.SemanticScholarId",
        "scholar.ExternalIds",
        "scholar.DOI",
        "scholar.Error",
    ),
    "openalex": (
        "OpenAlex.PublicationTypes",
        "OpenAlex.TypeCrossref",
        "OpenAlex.Genre",
        "OpenAlex.Id",
        "OpenAlex.Venue",
        "OpenAlex.MeshDescriptors",
        "OpenAlex.MeshQualifiers",
        "OpenAlex.Error",
    ),
    "crossref": (
        "crossref.Type",
        "crossref.Subtype",
        "crossref.Title",
        "crossref.Subtitle",
        "crossref.Subject",
        "crossref.Error",
    ),
    "pipeline_runtime": (
        "date_code",
        "Index",
        "PubMed.is_review",
        "scholar.is_review",
        "OpenAlex.is_review",
        "pipeline_version",
        "timestamp_utc",
    ),
}

_DOCUMENT_EXPORT_COLUMNS: tuple[str, ...] = (
    "PubMed.PMID",
    "PubMed.DOI",
    "PubMed.ArticleTitle",
    "PubMed.Abstract",
    "PubMed.JournalTitle",
    "PubMed.JournalISOAbbrev",
    "PubMed.Volume",
    "PubMed.Issue",
    "PubMed.StartPage",
    "PubMed.EndPage",
    "PubMed.ISSN",
    "PubMed.PublicationType",
    "PubMed.MeSH_Descriptors",
    "PubMed.MeSH_Qualifiers",
    "PubMed.ChemicalList",
    "PubMed.YearCompleted",
    "PubMed.MonthCompleted",
    "PubMed.DayCompleted",
    "PubMed.YearRevised",
    "PubMed.MonthRevised",
    "PubMed.DayRevised",
    "PubMed.Error",
    "scholar.PMID",
    "scholar.DOI",
    "scholar.PublicationTypes",
    "scholar.Venue",
    "scholar.SemanticScholarId",
    "scholar.ExternalIds",
    "scholar.Error",
    "OpenAlex.PMID",
    "OpenAlex.DOI",
    "OpenAlex.PublicationTypes",
    "OpenAlex.TypeCrossref",
    "OpenAlex.Genre",
    "OpenAlex.Venue",
    "OpenAlex.MeshDescriptors",
    "OpenAlex.MeshQualifiers",
    "OpenAlex.Id",
    "OpenAlex.Error",
    "crossref.DOI",
    "crossref.Type",
    "crossref.Subtype",
    "crossref.Title",
    "crossref.Subtitle",
    "crossref.Subject",
    "crossref.Error",
    "publication_types_normalised",
    "publication_review_score",
    "publication_experimental_score",
    "publication_class",
    "ChEMBL.document_chembl_id",
    "ChEMBL.title",
    "ChEMBL.abstract",
    "ChEMBL.doi",
    "ChEMBL.year",
    "ChEMBL.journal",
    "ChEMBL.journal_abbrev",
    "ChEMBL.volume",
    "ChEMBL.issue",
    "ChEMBL.first_page",
    "ChEMBL.last_page",
    "ChEMBL.pubmed_id",
    "ChEMBL.authors",
    "ChEMBL.source",
)


class DocumentSchema(DataFrameSchema):
    """Concrete :class:`pandera.DataFrameSchema` for document metadata."""

    BASE_COLUMNS: ClassVar[Mapping[str, pa.Column]] = _DOCUMENT_COLUMNS
    COLUMN_GROUPS: ClassVar[Mapping[str, tuple[str, ...]]] = _DOCUMENT_COLUMN_GROUPS
    ORDERED_COLUMNS: ClassVar[tuple[str, ...]] = tuple(_DOCUMENT_COLUMNS.keys())
    EXPORT_COLUMNS: ClassVar[tuple[str, ...]] = _DOCUMENT_EXPORT_COLUMNS
    BACKEND_REGISTRY: ClassVar = DataFrameSchema.BACKEND_REGISTRY

    def __init__(self) -> None:
        super().__init__(
            columns=OrderedDict((name, column) for name, column in self.BASE_COLUMNS.items()),
            coerce=True,
            ordered=True,
            strict=False,
        )

    @classmethod
    def get_backend(
        cls,
        check_obj: object | None = None,
        check_type: type[object] | None = None,
    ) -> Any:
        return DataFrameSchema.get_backend(check_obj, check_type)


DOCUMENT_SCHEMA = DocumentSchema()
"""Singleton instance of :class:`DocumentSchema` for reuse across modules."""

DOCUMENT_SCHEMA_COLUMNS: tuple[str, ...] = DocumentSchema.ORDERED_COLUMNS
"""Ordered column names defined by :class:`DocumentSchema`."""

DOCUMENT_COLUMN_GROUPS: Mapping[str, tuple[str, ...]] = DocumentSchema.COLUMN_GROUPS
"""Mapping of logical column groups to their ordered member columns."""

DOCUMENT_EXPORT_COLUMNS: tuple[str, ...] = DocumentSchema.EXPORT_COLUMNS
"""Default export projection for downstream CLI tools."""


def validate_document_frame(frame: pd.DataFrame) -> pd.DataFrame:
    """Validate a dataframe with the canonical document schema.

    Args:
        frame: DataFrame containing document metadata.

    Returns:
        The validated DataFrame with coerced dtypes where applicable.
    """

    return DOCUMENT_SCHEMA.validate(frame)


if __name__ == "__main__":  # pragma: no cover - developer utility
    for column_name in DOCUMENT_SCHEMA_COLUMNS:
        print(column_name)
