from __future__ import annotations

from collections.abc import Iterator, Mapping
from dataclasses import dataclass
from typing import Any, Protocol, TypedDict, runtime_checkable


class Record(TypedDict, total=False):
    """Минимально типизированная запись, определяется конфигурацией источника."""


@dataclass(slots=True)
class RequestContext:
    """Контекст запроса, прокидываемый через транспорт."""

    trace_id: str | None = None
    options: Mapping[str, Any] | None = None


@dataclass(slots=True)
class Pagination:
    """Параметры пагинации, передаваемые клиенту."""

    page_size: int | None = None
    cursor: str | None = None


@dataclass(slots=True)
class ClientRequest:
    """Единое описание запроса к внешнему источнику."""

    route: str
    params: Mapping[str, Any] | None = None
    context: RequestContext | None = None
    pagination: Pagination | None = None


@dataclass(slots=True)
class Page:
    """Пакет элементов с информацией о продолжении."""

    items: list[Record]
    next_cursor: str | None = None
    raw: Any | None = None


@runtime_checkable
class ExternalDataClient(Protocol):
    """Единый контракт для клиентов внешних источников данных."""

    def fetch_one(self, request: ClientRequest) -> Record | None:
        ...

    def fetch_many(self, request: ClientRequest) -> Iterator[Record]:
        ...

    def iter_pages(self, request: ClientRequest) -> Iterator[Page]:
        ...

    def metadata(self) -> Mapping[str, Any]:
        ...

    def close(self) -> None:
        ...
