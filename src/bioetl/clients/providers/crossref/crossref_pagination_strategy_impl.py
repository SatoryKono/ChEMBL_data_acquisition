from __future__ import annotations

from typing import Any, Iterator, Mapping

from bioetl.clients.base.http import BaseHttpClientABC
from bioetl.clients.base.interfaces import Page
from bioetl.clients.base.pagination import PaginationStrategyABC


class CrossrefPaginationStrategyImpl(PaginationStrategyABC):
    """Пагинация для ответов Crossref с вложенным ``message.items``."""

    def __init__(
        self,
        *,
        page_key: str = "message",
        items_key: str = "items",
        next_key: str | None = "next-cursor",
        page_param: str | None = "cursor",
    ) -> None:
        self.page_key = page_key
        self.items_key = items_key
        self.next_key = next_key
        self.page_param = page_param

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
        effective_page_key = page_key or self.page_key
        effective_next_key = next_key or self.next_key
        effective_page_param = page_param or self.page_param
        current_params = dict(params or {})

        while True:
            response = transport.get_json(endpoint, params=current_params or None)
            message = response.get(effective_page_key) if effective_page_key else response
            if not isinstance(message, Mapping):
                yield Page(items=[], next_cursor=None, raw=response)
                break

            items = list(message.get(self.items_key, []) or [])
            next_cursor = message.get(effective_next_key) if effective_next_key else None
            yield Page(items=items, next_cursor=next_cursor, raw=response)

            if not next_cursor or effective_page_param is None:
                break

            current_params = dict(current_params)
            current_params[effective_page_param] = next_cursor
