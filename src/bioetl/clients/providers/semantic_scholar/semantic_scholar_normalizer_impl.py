from __future__ import annotations

from typing import Any, Iterable, Iterator, Mapping

from bioetl.clients.base.normalizers import INormalizer


class SemanticScholarNormalizerImpl(INormalizer):
    """Нормализатор записей Semantic Scholar в доменную модель."""

    def normalize(self, record: Mapping[str, Any]) -> Mapping[str, Any]:
        authors = self._normalize_authors(record.get("authors"))
        return {
            "id": record.get("paperId") or "",
            "title": record.get("title") or "",
            "authors": authors,
            "year": record.get("year"),
            "abstract": record.get("abstract") or "",
            "citation_count": record.get("citationCount"),
        }

    def normalize_batch(
        self, records: Iterable[Mapping[str, Any]]
    ) -> Iterator[Mapping[str, Any]]:
        for record in records:
            yield self.normalize(record)

    @staticmethod
    def _normalize_authors(value: Any) -> list[Mapping[str, Any]]:
        if not isinstance(value, Iterable) or isinstance(value, (str, bytes)):
            return []

        normalized: list[Mapping[str, Any]] = []
        for author in value:
            if isinstance(author, Mapping):
                normalized.append(
                    {
                        "name": author.get("name") or "",
                        "authorId": author.get("authorId") or None,
                    }
                )
        return normalized
