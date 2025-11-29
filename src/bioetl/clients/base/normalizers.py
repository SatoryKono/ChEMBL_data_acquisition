from __future__ import annotations

from typing import Any, Iterable, Iterator, Mapping, Protocol


class INormalizer(Protocol):
    """Интерфейс нормализаторов записей."""

    def normalize(self, record: Mapping[str, Any]) -> Mapping[str, Any]:
        """Нормализовать одну запись."""

    def normalize_batch(
        self, records: Iterable[Mapping[str, Any]]
    ) -> Iterator[Mapping[str, Any]]:
        """Нормализовать последовательность записей."""

        for record in records:
            yield self.normalize(record)


class IdentityNormalizerImpl(INormalizer):
    """Нормализатор, возвращающий копию исходной записи без изменений."""

    def normalize(self, record: Mapping[str, Any]) -> Mapping[str, Any]:
        return dict(record)

    def normalize_batch(
        self, records: Iterable[Mapping[str, Any]]
    ) -> Iterator[Mapping[str, Any]]:
        for record in records:
            yield self.normalize(record)
