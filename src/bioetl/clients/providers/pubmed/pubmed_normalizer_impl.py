from __future__ import annotations

from datetime import datetime
from typing import Any, Mapping, Sequence

from bioetl.clients.base.normalizers import INormalizer


def _concat_abstract(abstract: Any) -> str | None:
    if abstract is None:
        return None
    if isinstance(abstract, str):
        return abstract
    if isinstance(abstract, Mapping):
        text = abstract.get("AbstractText")
        return _concat_abstract(text)
    if isinstance(abstract, Sequence) and not isinstance(abstract, (str, bytes)):
        return "\n".join(filter(None, (_concat_abstract(item) or "" for item in abstract))).strip() or None
    return str(abstract)


def _format_author(author: Mapping[str, Any]) -> str:
    if "CollectiveName" in author:
        return str(author.get("CollectiveName"))
    last = author.get("LastName")
    fore = author.get("ForeName")
    if last and fore:
        return f"{fore} {last}"
    if last:
        return str(last)
    if fore:
        return str(fore)
    return ""


def _extract_pub_date(article: Mapping[str, Any]) -> str | None:
    article_dates = article.get("ArticleDate")
    if isinstance(article_dates, Sequence) and article_dates:
        return _compose_date(article_dates[0])

    journal_issue = article.get("Journal", {}).get("JournalIssue", {}) if isinstance(article, Mapping) else {}
    pub_date = journal_issue.get("PubDate") if isinstance(journal_issue, Mapping) else None
    if isinstance(pub_date, Mapping):
        date = _compose_date(pub_date)
        if date:
            return date
        medline_date = pub_date.get("MedlineDate")
        if isinstance(medline_date, str):
            return medline_date
    return None


def _compose_date(parts: Mapping[str, Any]) -> str | None:
    year = parts.get("Year")
    month = parts.get("Month")
    day = parts.get("Day")
    if year and month and day:
        try:
            return datetime(int(year), int(month), int(day)).date().isoformat()
        except ValueError:
            return None
    if year and month:
        return f"{year}-{str(month).zfill(2)}"
    if year:
        return str(year)
    return None


class PubMedNormalizerImpl(INormalizer):
    """Нормализует записи PubMed в унифицированный формат."""

    def normalize(self, record: Mapping[str, Any] | Any) -> Mapping[str, Any]:
        if not isinstance(record, Mapping):
            record = {"id": record}

        pmid = self._extract_pmid(record)
        article = self._extract_article(record)

        authors_list = self._extract_authors(article)
        journal = self._extract_journal(article)
        abstract = _concat_abstract(self._extract_abstract(article))
        pub_date = _extract_pub_date(article) if article else None

        title = None
        if article:
            raw_title = article.get("ArticleTitle")
            title = raw_title if not isinstance(raw_title, Mapping) else raw_title.get("_text")

        return {
            "id": pmid,
            "title": title,
            "authors": authors_list,
            "journal": journal,
            "abstract": abstract,
            "pub_date": pub_date,
        }

    def _extract_pmid(self, record: Mapping[str, Any]) -> str | None:
        if "PMID" in record:
            return str(record.get("PMID")) if record.get("PMID") is not None else None
        if "id" in record:
            return str(record["id"])
        if "value" in record:
            return str(record["value"])

        medline = self._extract_medline(record)
        if medline and "PMID" in medline:
            return str(medline.get("PMID")) if medline.get("PMID") is not None else None
        return None

    def _extract_medline(self, record: Mapping[str, Any]) -> Mapping[str, Any] | None:
        if "MedlineCitation" in record and isinstance(record.get("MedlineCitation"), Mapping):
            return record.get("MedlineCitation")  # type: ignore[return-value]
        if "PubmedArticle" in record and isinstance(record.get("PubmedArticle"), Mapping):
            return record.get("PubmedArticle", {}).get("MedlineCitation")
        if "PubmedArticleSet" in record:
            article_set = record.get("PubmedArticleSet")
            if isinstance(article_set, Mapping):
                articles = article_set.get("PubmedArticle")
                if isinstance(articles, Mapping):
                    return articles.get("MedlineCitation")
                if isinstance(articles, Sequence) and articles:
                    first = articles[0]
                    if isinstance(first, Mapping):
                        return first.get("MedlineCitation")
        return None

    def _extract_article(self, record: Mapping[str, Any]) -> Mapping[str, Any] | None:
        medline = self._extract_medline(record)
        if medline and isinstance(medline.get("Article"), Mapping):
            return medline.get("Article")  # type: ignore[return-value]
        return None

    def _extract_authors(self, article: Mapping[str, Any] | None) -> list[str]:
        if not article:
            return []
        authors = article.get("AuthorList")
        if not isinstance(authors, Sequence):
            return []
        formatted: list[str] = []
        for raw_author in authors:
            if isinstance(raw_author, Mapping):
                name = _format_author(raw_author)
                if name:
                    formatted.append(name)
        return formatted

    def _extract_journal(self, article: Mapping[str, Any] | None) -> str | None:
        if not article:
            return None
        journal = article.get("Journal")
        if isinstance(journal, Mapping):
            title = journal.get("Title")
            return str(title) if title is not None else None
        return None

    def _extract_abstract(self, article: Mapping[str, Any] | None) -> Any:
        if not article:
            return None
        abstract = article.get("Abstract")
        if isinstance(abstract, Mapping):
            return abstract.get("AbstractText")
        return abstract
