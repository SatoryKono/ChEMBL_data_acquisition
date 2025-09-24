"""Parsing helpers for PubMed records."""

from __future__ import annotations

from collections.abc import Iterable
from typing import Any
from xml.etree import ElementTree as ET

__all__ = [
    "text_or_none",
    "combine",
    "find_one",
    "find_all",
    "parse_pubmed_article",
    "EMPTY_PUBMED",
]


def text_or_none(node: ET.Element | None) -> str | None:
    """Return stripped text of an XML node if present."""
    if node is not None and node.text is not None:
        return node.text.strip()
    return None


def combine(items: Iterable[str]) -> str:
    """Combine non-empty items into a pipe-separated string without spaces."""
    return "|".join(x for x in items if x)


def find_one(node: ET.Element | None, xpath: str) -> ET.Element | None:
    """Safe wrapper around :func:`Element.find`."""
    return node.find(xpath) if node is not None else None


def find_all(node: ET.Element | None, xpath: str) -> list[ET.Element]:
    """Safe wrapper around :func:`Element.findall` returning a list."""
    return node.findall(xpath) if node is not None else []


def _clean_doi(raw: str | None) -> str | None:
    """Normalise DOI strings by stripping known prefixes."""

    if not raw:
        return None

    doi = raw.strip()
    lower = doi.lower()
    if lower.startswith("doi:"):
        doi = doi[4:].strip()
        lower = doi.lower()

    for prefix in ("https://doi.org/", "http://doi.org/", "doi.org/"):
        if lower.startswith(prefix):
            doi = doi[len(prefix) :].strip()
            break
        lower = doi.lower()

    return doi or None


def _extract_doi(
    article: ET.Element | None, pubmed_data: ET.Element | None
) -> str | None:
    """Extract DOI from ``Article`` or ``PubmedData`` blocks."""

    candidates: list[str] = []

    if article is not None:
        for node in find_all(article, "./ArticleIdList/ArticleId"):
            if node.get("IdType") == "doi" and node.text:
                candidates.append(node.text)
        for node in find_all(article, "./ELocationID"):
            if node.get("EIdType") == "doi" and node.text:
                candidates.append(node.text)

    if pubmed_data is not None:
        for node in find_all(pubmed_data, "./ArticleIdList/ArticleId"):
            if node.get("IdType") == "doi" and node.text:
                candidates.append(node.text)

    for candidate in candidates:
        doi = _clean_doi(candidate)
        if doi:
            return doi

    return None


def parse_pubmed_article(art: ET.Element) -> dict[str, Any]:
    """Parse ``PubMedArticle`` into a dictionary of selected fields."""
    mc = find_one(art, "./MedlineCitation")
    article = find_one(mc, "./Article") if mc is not None else None
    journal = find_one(article, "./Journal") if article is not None else None
    journal_issue = find_one(journal, "./JournalIssue") if journal is not None else None
    pagination = find_one(article, "./Pagination") if article is not None else None
    pubmed_data = find_one(art, "./PubmedData")

    pmid = text_or_none(find_one(mc, "./PMID")) if mc is not None else None

    # DOI: prefer Article block, fallback to PubmedData block
    doi = _extract_doi(article, pubmed_data)

    article_title = (
        text_or_none(find_one(article, "./ArticleTitle"))
        if article is not None
        else None
    )

    article_abstract = None
    if article is not None:
        segments = find_all(article, "./Abstract/AbstractText")
        if segments:
            parts: list[str] = []
            for seg in segments:
                seg_text = text_or_none(seg)
                if seg_text:
                    label = seg.get("Label")
                    parts.append(f"{label}: {seg_text}" if label else seg_text)
            article_abstract = " ".join(parts) if parts else None
        if article_abstract is None:
            article_abstract = text_or_none(
                find_one(article, "./Abstract/AbstractText")
            )

    journal_title = text_or_none(find_one(journal, "./Title"))
    issn = text_or_none(find_one(journal, "./ISSN"))
    volume = text_or_none(find_one(journal_issue, "./Volume"))
    issue = text_or_none(find_one(journal_issue, "./Issue"))

    start_page = text_or_none(find_one(pagination, "./StartPage"))
    end_page = text_or_none(find_one(pagination, "./EndPage"))

    pubtypes = [
        p
        for p in (
            text_or_none(x)
            for x in find_all(article, "./PublicationTypeList/PublicationType")
        )
        if p
    ]

    mh_list = find_one(mc, "./MeshHeadingList")
    mesh_descriptors: list[str] = []
    mesh_qualifiers: list[str] = []
    if mh_list is not None:
        for mh in mh_list.findall("./MeshHeading"):
            d = text_or_none(find_one(mh, "./DescriptorName"))
            if d:
                mesh_descriptors.append(d)
            for q in mh.findall("./QualifierName"):
                qt = text_or_none(q)
                if qt:
                    mesh_qualifiers.append(qt)

    chemical_list: list[str] = []
    chem_list_node = find_one(mc, "./ChemicalList")
    if chem_list_node is not None:
        for chem in chem_list_node.findall("./Chemical"):
            nos = text_or_none(find_one(chem, "./NameOfSubstance"))
            if nos:
                chemical_list.append(nos)

    dr = find_one(mc, "./DateRevised")
    year_revised = text_or_none(find_one(dr, "./Year")) if dr is not None else ""
    month_revised = text_or_none(find_one(dr, "./Month")) if dr is not None else ""
    day_revised = text_or_none(find_one(dr, "./Day")) if dr is not None else ""

    dc = find_one(mc, "./DateCompleted")
    year_completed = text_or_none(find_one(dc, "./Year")) if dc is not None else ""
    month_completed = text_or_none(find_one(dc, "./Month")) if dc is not None else ""
    day_completed = text_or_none(find_one(dc, "./Day")) if dc is not None else ""

    return {
        "PubMed.PMID": pmid or "",
        "PubMed.DOI": doi or "",
        "PubMed.ArticleTitle": article_title or "-",
        "PubMed.Abstract": article_abstract or "-",
        "PubMed.JournalTitle": journal_title or "-",
        "PubMed.ISSN": issn or "",
        "PubMed.Volume": volume or "0",
        "PubMed.Issue": issue or "0",
        "PubMed.StartPage": start_page or "0",
        "PubMed.EndPage": end_page or "0",
        "PubMed.PublicationType": combine(pubtypes) or "unknown",
        "PubMed.MeSH_Descriptors": combine(mesh_descriptors) or "unknown",
        "PubMed.MeSH_Qualifiers": combine(mesh_qualifiers) or "unknown",
        "PubMed.ChemicalList": combine(chemical_list) or "unknown",
        "PubMed.DayRevised": day_revised or "0",
        "PubMed.MonthRevised": month_revised or "0",
        "PubMed.YearRevised": year_revised or "0",
        "PubMed.YearCompleted": year_completed or "0",
        "PubMed.MonthCompleted": month_completed or "0",
        "PubMed.DayCompleted": day_completed or "0",
    }


EMPTY_PUBMED: dict[str, str] = {
    "PubMed.PMID": "",
    "PubMed.DOI": "",
    "PubMed.ArticleTitle": "",
    "PubMed.Abstract": "",
    "PubMed.JournalTitle": "",
    "PubMed.ISSN": "",
    "PubMed.Volume": "",
    "PubMed.Issue": "",
    "PubMed.StartPage": "",
    "PubMed.EndPage": "",
    "PubMed.PublicationType": "",
    "PubMed.MeSH_Descriptors": "",
    "PubMed.MeSH_Qualifiers": "",
    "PubMed.ChemicalList": "",
    "PubMed.DayRevised": "",
    "PubMed.MonthRevised": "",
    "PubMed.YearRevised": "",
    "PubMed.YearCompleted": "",
    "PubMed.MonthCompleted": "",
    "PubMed.DayCompleted": "",
    "PubMed.Error": "",
}
