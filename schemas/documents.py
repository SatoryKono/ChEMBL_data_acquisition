"""Schema definitions for document data."""

from __future__ import annotations

import pandera.pandas as pa

from library.document_schema import DOCUMENT_SCHEMA_COLUMNS

# ``pa.Any`` is not available in the current pandera version; ``object`` with
# ``coerce=True`` permits arbitrary dtypes.

_COLUMN_DEFINITIONS: dict[str, pa.Column] = {
    "document_chembl_id": pa.Column(str, required=True, nullable=True),
    "doi": pa.Column(str, required=False, nullable=True),
    "title": pa.Column(str, required=False, nullable=True),
    "Index": pa.Column(object, required=False, nullable=True, coerce=True),
    "OpenAlex.Error": pa.Column(object, required=False, nullable=True, coerce=True),
    "OpenAlex.Genre": pa.Column(object, required=False, nullable=True, coerce=True),
    "OpenAlex.Id": pa.Column(str, required=False, nullable=True),
    "OpenAlex.MeshDescriptors": pa.Column(str, required=False, nullable=True),
    "OpenAlex.MeshQualifiers": pa.Column(
        object, required=False, nullable=True, coerce=True
    ),
    "OpenAlex.PublicationTypes": pa.Column(
        str, required=False, nullable=True
    ),
    "OpenAlex.TypeCrossref": pa.Column(str, required=False, nullable=True),
    "OpenAlex.Venue": pa.Column(object, required=False, nullable=True, coerce=True),
    "OpenAlex.is_review": pa.Column(
        object, required=False, nullable=True, coerce=True
    ),
    "PubMed.Abstract": pa.Column(str, required=False, nullable=True),
    "PubMed.ArticleTitle": pa.Column(str, required=False, nullable=True),
    "PubMed.ChemicalList": pa.Column(str, required=False, nullable=True),
    "PubMed.DOI": pa.Column(object, required=False, nullable=True, coerce=True),
    "PubMed.DayCompleted": pa.Column(
        object, required=False, nullable=True, coerce=True
    ),
    "PubMed.DayRevised": pa.Column(
        object, required=False, nullable=True, coerce=True
    ),
    "PubMed.EndPage": pa.Column(object, required=False, nullable=True, coerce=True),
    "PubMed.Error": pa.Column(object, required=False, nullable=True, coerce=True),
    "PubMed.ISSN": pa.Column(str, required=False, nullable=True),
    "PubMed.Issue": pa.Column(object, required=False, nullable=True, coerce=True),
    "PubMed.JournalTitle": pa.Column(str, required=False, nullable=True),
    "PubMed.MeSH_Descriptors": pa.Column(str, required=False, nullable=True),
    "PubMed.MeSH_Qualifiers": pa.Column(str, required=False, nullable=True),
    "PubMed.MonthCompleted": pa.Column(
        object, required=False, nullable=True, coerce=True
    ),
    "PubMed.MonthRevised": pa.Column(
        object, required=False, nullable=True, coerce=True
    ),
    "PubMed.PMID": pa.Column(object, required=False, nullable=True, coerce=True),
    "PubMed.PublicationType": pa.Column(
        str, required=False, nullable=True
    ),
    "PubMed.StartPage": pa.Column(
        object, required=False, nullable=True, coerce=True
    ),
    "PubMed.Volume": pa.Column(object, required=False, nullable=True, coerce=True),
    "PubMed.YearCompleted": pa.Column(
        object, required=False, nullable=True, coerce=True
    ),
    "PubMed.YearRevised": pa.Column(
        object, required=False, nullable=True, coerce=True
    ),
    "PubMed.is_review": pa.Column(
        object, required=False, nullable=True, coerce=True
    ),
    "abstract": pa.Column(str, required=False, nullable=True),
    "authors": pa.Column(str, required=False, nullable=True),
    "crossref.Error": pa.Column(str, required=False, nullable=True),
    "crossref.Subject": pa.Column(
        object, required=False, nullable=True, coerce=True
    ),
    "crossref.Subtitle": pa.Column(
        object, required=False, nullable=True, coerce=True
    ),
    "crossref.Subtype": pa.Column(
        object, required=False, nullable=True, coerce=True
    ),
    "crossref.Title": pa.Column(object, required=False, nullable=True, coerce=True),
    "crossref.Type": pa.Column(object, required=False, nullable=True, coerce=True),
    "date_code": pa.Column(str, required=False, nullable=True),
    "first_page": pa.Column(object, required=False, nullable=True, coerce=True),
    "issue": pa.Column(object, required=False, nullable=True, coerce=True),
    "journal": pa.Column(str, required=False, nullable=True),
    "journal_abbrev": pa.Column(str, required=False, nullable=True),
    "last_page": pa.Column(object, required=False, nullable=True, coerce=True),
    "pubmed_id": pa.Column(object, required=False, nullable=True, coerce=True),
    "scholar.DOI": pa.Column(object, required=False, nullable=True, coerce=True),
    "scholar.Error": pa.Column(str, required=False, nullable=True),
    "scholar.ExternalIds": pa.Column(
        object, required=False, nullable=True, coerce=True
    ),
    "scholar.PMID": pa.Column(object, required=False, nullable=True, coerce=True),
    "scholar.PublicationTypes": pa.Column(
        str, required=False, nullable=True
    ),
    "scholar.SemanticScholarId": pa.Column(
        object, required=False, nullable=True, coerce=True
    ),
    "scholar.Venue": pa.Column(object, required=False, nullable=True, coerce=True),
    "scholar.is_review": pa.Column(
        object, required=False, nullable=True, coerce=True
    ),
    "source": pa.Column(str, required=False, nullable=True),
    "volume": pa.Column(object, required=False, nullable=True, coerce=True),
    "year": pa.Column(object, required=False, nullable=True, coerce=True),
    "publication_types_normalised": pa.Column(
        str, required=False, nullable=True
    ),
    "publication_type_score_review": pa.Column(
        object, required=False, nullable=True, coerce=True
    ),
    "publication_type_score_experimental": pa.Column(
        object, required=False, nullable=True, coerce=True
    ),
    "publication_type_score_unknown": pa.Column(
        object, required=False, nullable=True, coerce=True
    ),
    "publication_class": pa.Column(str, required=False, nullable=True),
}

DocumentsSchema: pa.DataFrameSchema = pa.DataFrameSchema(
    {name: _COLUMN_DEFINITIONS[name] for name in DOCUMENT_SCHEMA_COLUMNS}
)

"""pa.DataFrameSchema: Validation schema for documents."""
