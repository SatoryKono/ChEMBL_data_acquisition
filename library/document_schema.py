"""Shared column definitions for the document pipeline.

The original Power Query implementation relied on a number of hard coded
column lists to keep the output schema deterministic.  The data acquisition
scripts translate the logic to ``pandas`` and re-use the same collections to
guarantee stable ordering across the different execution modes.
"""

from __future__ import annotations

from typing import Final

# Column order expected by the downstream schema validation.  The names are
# derived from the dictionary files used by the legacy ETL and extended with
# the publication type metrics calculated during post-processing.
DOCUMENT_SCHEMA_COLUMNS: Final[tuple[str, ...]] = (
    "Index",
    "OpenAlex.Error",
    "OpenAlex.Genre",
    "OpenAlex.Id",
    "OpenAlex.MeshDescriptors",
    "OpenAlex.MeshQualifiers",
    "OpenAlex.PublicationTypes",
    "OpenAlex.TypeCrossref",
    "OpenAlex.Venue",
    "OpenAlex.is_review",
    "PubMed.Abstract",
    "PubMed.ArticleTitle",
    "PubMed.ChemicalList",
    "PubMed.DOI",
    "PubMed.DayCompleted",
    "PubMed.DayRevised",
    "PubMed.EndPage",
    "PubMed.Error",
    "PubMed.ISSN",
    "PubMed.Issue",
    "PubMed.JournalTitle",
    "PubMed.MeSH_Descriptors",
    "PubMed.MeSH_Qualifiers",
    "PubMed.MonthCompleted",
    "PubMed.MonthRevised",
    "PubMed.PMID",
    "PubMed.PublicationType",
    "PubMed.StartPage",
    "PubMed.Volume",
    "PubMed.YearCompleted",
    "PubMed.YearRevised",
    "PubMed.is_review",
    "abstract",
    "authors",
    "crossref.Error",
    "crossref.Subject",
    "crossref.Subtitle",
    "crossref.Subtype",
    "crossref.Title",
    "crossref.Type",
    "date_code",
    "document_chembl_id",
    "doi",
    "first_page",
    "issue",
    "journal",
    "journal_abbrev",
    "last_page",
    "pubmed_id",
    "scholar.DOI",
    "scholar.Error",
    "scholar.ExternalIds",
    "scholar.PMID",
    "scholar.PublicationTypes",
    "scholar.SemanticScholarId",
    "scholar.Venue",
    "scholar.is_review",
    "source",
    "title",
    "volume",
    "year",
    "publication_types_normalised",
    "publication_type_score_review",
    "publication_type_score_experimental",
    "publication_type_score_unknown",
    "publication_class",
)

# Columns provided by the ChEMBL API.  ``get_document_data`` relies on the list
# to reindex partial responses and guarantee a deterministic output order when
# ChEMBL is the only source.
CH_EMBL_COLUMNS: Final[tuple[str, ...]] = (
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
)

# Subset of columns originating from the Semantic Scholar API.  The combined
# dataset aligns the field names with the legacy export to ease downstream
# comparisons.
SEMANTIC_SCHOLAR_COLUMNS: Final[tuple[str, ...]] = (
    "scholar.PMID",
    "scholar.Venue",
    "scholar.SemanticScholarId",
    "scholar.ExternalIds",
    "scholar.DOI",
    "scholar.PublicationTypes",
    "scholar.is_review",
    "scholar.Error",
)

__all__ = [
    "DOCUMENT_SCHEMA_COLUMNS",
    "CH_EMBL_COLUMNS",
    "SEMANTIC_SCHOLAR_COLUMNS",
]
