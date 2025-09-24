import pandas as pd

from library.document_pipeline import (
    DOCUMENT_SCHEMA_COLUMNS,
    build_dataframe,
    build_quality_report,
    merge_metadata,
    merge_with_chembl,
)


def test_merge_metadata_normalises_fields() -> None:
    """``merge_metadata`` normalises DOI and computes publication scores."""

    pubmed = {
        "PubMed.DOI": "DOI:10.1000/ABCDEF",
        "PubMed.PublicationType": "Review",
        "PubMed.PMID": "123",
    }
    scholar = {
        "scholar.DOI": "https://doi.org/10.1000/XYZ",
        "scholar.PublicationTypes": "Comparative Study",
    }
    openalex = {"OpenAlex.PublicationTypes": "journal-article"}
    crossref = {"crossref.Type": "journal-article"}

    merged = merge_metadata(pubmed, scholar, openalex, crossref)

    assert merged["doi_normalised"] == "10.1000/abcdef"
    assert merged["doi"] == "10.1000/abcdef"
    assert merged["scholar.DOI"] == "10.1000/xyz"
    assert (
        merged["publication_types_normalised"]
        == "comparative study; journal-article; review"
    )
    assert merged["publication_type_score_review"] > 0
    assert merged["publication_type_score_experimental"] > 0
    assert merged["publication_type_score_unknown"] == 0
    assert merged["publication_class"] == "review"


def test_build_dataframe_orders_columns() -> None:
    """``build_dataframe`` respects the schema column ordering."""

    data = [{"extra": "x", "document_chembl_id": "CHEMBL1", "doi": "10/1"}]
    df = build_dataframe(data)
    expected_head = [c for c in DOCUMENT_SCHEMA_COLUMNS if c in df.columns]
    assert list(df.columns)[: len(expected_head)] == expected_head


def test_merge_with_chembl_aligns_pubmed_ids() -> None:
    """``merge_with_chembl`` joins metadata on normalised PubMed identifiers."""

    chembl_df = pd.DataFrame(
        {
            "document_chembl_id": ["CHEMBL1"],
            "pubmed_id": [123],
            "title": ["T"],
        }
    )
    meta_df = pd.DataFrame(

        {
            "document_chembl_id": ["IGNORED"],
            "title": ["Other"],
            "PubMed.PMID": ["123"],
            "PubMed.DOI": ["10.1/abc"],
        }

    )

    merged = merge_with_chembl(chembl_df, meta_df)
    assert merged["PubMed.DOI"].iloc[0] == "10.1/abc"

    assert "document_chembl_id_x" not in merged.columns
    assert "document_chembl_id_y" not in merged.columns
    assert merged["document_chembl_id"].iloc[0] == "CHEMBL1"



def test_build_quality_report_counts() -> None:
    """Quality report exposes coverage and error statistics."""

    df = pd.DataFrame(
        {
            "doi": ["10.1/abc", ""],
            "publication_class": ["review", ""],
            "PubMed.Error": ["", "err"],
            "scholar.Error": ["", ""],
            "OpenAlex.Error": ["", ""],
            "crossref.Error": ["oops", ""],
        }
    )

    report = build_quality_report(df)
    assert report["rows_total"] == 2
    assert report["doi_coverage"] == 0.5
    assert report["publication_class_counts"]["review"] == 1
    assert report["error_counts"]["pubmed"] == 1
    assert report["error_counts"]["crossref"] == 1
