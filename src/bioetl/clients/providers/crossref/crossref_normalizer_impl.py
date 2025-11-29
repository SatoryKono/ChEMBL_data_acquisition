from __future__ import annotations

from typing import Any, Mapping, Sequence

from bioetl.clients.base.normalizers import INormalizer


def _normalize_doi(raw_doi: str | None) -> str | None:
    if not raw_doi:
        return None
    return f"doi:{raw_doi}" if not raw_doi.startswith("doi:") else raw_doi


def _first_title(value: Any) -> str | None:
    if isinstance(value, str):
        return value
    if isinstance(value, Sequence) and value:
        first = value[0]
        return str(first) if first is not None else None
    return None


def _collect_authors(authors: Any) -> list[Mapping[str, Any]]:
    if not isinstance(authors, Sequence):
        return []

    normalized: list[Mapping[str, Any]] = []
    for author in authors:
        if not isinstance(author, Mapping):
            continue
        given = str(author.get("given")) if author.get("given") is not None else ""
        family = str(author.get("family")) if author.get("family") is not None else ""
        full_name = " ".join(part for part in (given, family) if part).strip()
        affiliations = [
            aff.get("name")
            for aff in author.get("affiliation", [])
            if isinstance(aff, Mapping) and aff.get("name")
        ]
        entry: dict[str, Any] = {"name": full_name or given or family}
        if affiliations:
            entry["affiliations"] = affiliations
        normalized.append(entry)
    return normalized


def _format_date(parts: Sequence[int] | None) -> str | None:
    if not parts:
        return None
    padded = []
    for idx, value in enumerate(parts):
        if idx == 0:
            padded.append(f"{int(value):04d}")
        else:
            padded.append(f"{int(value):02d}")
    return "-".join(padded)


def _extract_date(message: Mapping[str, Any]) -> str | None:
    for key in ("published-print", "published-online", "issued"):
        container = message.get(key)
        if isinstance(container, Mapping):
            date_parts = container.get("date-parts")
            if isinstance(date_parts, Sequence) and date_parts:
                first = date_parts[0]
                if isinstance(first, Sequence):
                    return _format_date(list(first))
    return None


class CrossrefNormalizerImpl(INormalizer):
    """Нормализатор записей Crossref в упрощенную доменную модель."""

    def normalize(self, record: Mapping[str, Any]) -> Mapping[str, Any]:
        message = record.get("message") if isinstance(record, Mapping) else None
        payload = message if isinstance(message, Mapping) else record

        doi = _normalize_doi(payload.get("DOI")) if isinstance(payload, Mapping) else None
        title = _first_title(payload.get("title") if isinstance(payload, Mapping) else None)
        authors = _collect_authors(payload.get("author")) if isinstance(payload, Mapping) else []
        journal = None
        if isinstance(payload, Mapping):
            container_titles = payload.get("container-title")
            journal = _first_title(container_titles)
        published_date = _extract_date(payload) if isinstance(payload, Mapping) else None

        return {
            "id": doi,
            "title": title,
            "authors": authors,
            "journal": journal,
            "published_date": published_date,
            "references_count": payload.get("references-count") if isinstance(payload, Mapping) else None,
            "abstract": payload.get("abstract") if isinstance(payload, Mapping) else None,
        }
