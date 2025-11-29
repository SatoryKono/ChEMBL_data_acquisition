from __future__ import annotations

from typing import Any, Iterable, Mapping, MutableMapping

from bioetl.clients.base import INormalizer


class OpenAlexNormalizerImpl(INormalizer):
    """Нормализует записи OpenAlex в упрощённый доменный формат."""

    def __init__(self) -> None:
        self._fields_map = {
            "id": "id",
            "display_name": "title",
            "title": "title",
            "authorships": "authors",
        }

    def _extract_authors(self, authorship_entries: Iterable[Any]) -> list[str]:
        authors: list[str] = []
        for authorship in authorship_entries:
            if not isinstance(authorship, Mapping):
                continue
            author = authorship.get("author")
            display_name: Any | None = None
            if isinstance(author, Mapping):
                display_name = author.get("display_name") or author.get("name")
            if isinstance(display_name, str) and display_name:
                authors.append(display_name)
        return authors

    def normalize(self, record: Mapping[str, Any]) -> MutableMapping[str, Any]:
        normalized: MutableMapping[str, Any] = {"id": "", "title": "", "authors": []}

        if not isinstance(record, Mapping):
            return normalized

        for key, normalized_key in self._fields_map.items():
            if key not in record:
                continue
            value = record.get(key)
            if normalized_key == "authors":
                normalized["authors"] = (
                    self._extract_authors(value)
                    if isinstance(value, Iterable) and not isinstance(value, (str, bytes))
                    else []
                )
                continue

            if value is None:
                continue
            normalized[normalized_key] = value if not isinstance(value, Mapping) else dict(value)

        if normalized["id"]:
            normalized["id"] = str(normalized["id"])
        if normalized["title"]:
            normalized["title"] = str(normalized["title"])

        return normalized
