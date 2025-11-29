from __future__ import annotations

import abc
from typing import Any, Iterator, Mapping

from bioetl.clients.base.http import BaseHttpClientABC


class PaginationStrategyABC(abc.ABC):
    """Абстракция стратегии пагинации для HTTP-транспортов."""

    @abc.abstractmethod
    def iter_pages(
        self,
        initial_response_or_args: Any,
        transport: BaseHttpClientABC,
        *,
        endpoint: str,
        params: Mapping[str, Any] | None = None,
        page_key: str | None = None,
        next_key: str | None = None,
        page_param: str | None = None,
    ) -> Iterator[Any]:
        """Итерирует страницы, возвращаемые транспортом."""


class TransportPaginationStrategyImpl(PaginationStrategyABC):
    """Стратегия пагинации, делегирующая логику транспорту."""

    def iter_pages(
        self,
        initial_response_or_args: Any,
        transport: BaseHttpClientABC,
        *,
        endpoint: str,
        params: Mapping[str, Any] | None = None,
        page_key: str | None = None,
        next_key: str | None = None,
        page_param: str | None = None,
    ) -> Iterator[Any]:
        yield from transport.paginate_json(
            endpoint,
            params=params,
            page_key=page_key,
            next_key=next_key,
            page_param=page_param,
        )
