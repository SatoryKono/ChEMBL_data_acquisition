from __future__ import annotations

import abc
from typing import Any, Callable, Iterator, Mapping, TypeVar

from bioetl.clients.base.http import BaseHttpClientABC
from bioetl.clients.base.interfaces import (
    DataProviderError,
    Page,
    PaginationParams,
    RequestContext,
)

T = TypeVar("T")


class BaseDataProviderABC(abc.ABC):
    def __init__(
        self,
        http: BaseHttpClientABC,
        *,
        logger: Any | None = None,
        options: Mapping[str, Any] | None = None,
    ) -> None:
        self._http = http
        self._logger = logger
        self._options = options or {}

    def _wrap_callable(self, fn: Callable[..., T], context: RequestContext | None = None) -> T:
        try:
            return fn()
        except DataProviderError:
            raise
        except Exception as exc:  # noqa: BLE001
            if self._logger is not None:
                self._logger.exception("Data provider call failed", extra={"context": context})
            raise DataProviderError(str(exc)) from exc

    def _wrap_iterator(self, it: Iterator[T], context: RequestContext | None = None) -> Iterator[T]:
        try:
            for item in it:
                yield item
        except DataProviderError:
            raise
        except Exception as exc:  # noqa: BLE001
            if self._logger is not None:
                self._logger.exception("Data provider iteration failed", extra={"context": context})
            raise DataProviderError(str(exc)) from exc

    @abc.abstractmethod
    def fetch_one(
        self,
        identifier: str | Mapping[str, Any],
        *,
        params: Mapping[str, Any] | None = None,
        context: RequestContext | None = None,
    ) -> Iterator[Mapping[str, Any]]:
        ...

    @abc.abstractmethod
    def iter_pages(
        self,
        *,
        params: Mapping[str, Any] | None = None,
        page_size: int | None = None,
        context: RequestContext | None = None,
    ) -> Iterator[Page]:
        ...

    @abc.abstractmethod
    def fetch_many(
        self,
        *,
        params: Mapping[str, Any] | None = None,
        page_size: int | None = None,
        context: RequestContext | None = None,
    ) -> Iterator[Mapping[str, Any]]:
        ...

    def metadata(self) -> Mapping[str, Any]:
        return getattr(self._http, "metadata", {}) or {}

    def close(self) -> None:
        self._http.close()


class PagedDataProviderABC(BaseDataProviderABC):
    @abc.abstractmethod
    def _iterate_pages_impl(
        self,
        params: Mapping[str, Any] | None,
        pagination: PaginationParams,
        context: RequestContext | None = None,
    ) -> Iterator[Any]:
        ...

    def iter_pages(
        self,
        *,
        params: Mapping[str, Any] | None = None,
        page_size: int | None = None,
        context: RequestContext | None = None,
    ) -> Iterator[Page]:
        pagination = PaginationParams(
            page_size=page_size if page_size is not None else self._options.get("page_size"),
            page_key=self._options.get("page_key", "items"),
            next_key=self._options.get("next_key", "next"),
            page_param=self._options.get("page_param", "page"),
        )

        def generator() -> Iterator[Page]:
            for raw in self._iterate_pages_impl(params, pagination=pagination, context=context):
                yield self._normalize_page_payload(raw, pagination.page_key, pagination.next_key)

        return self._wrap_iterator(generator(), context=context)

    def fetch_many(
        self,
        *,
        params: Mapping[str, Any] | None = None,
        page_size: int | None = None,
        context: RequestContext | None = None,
    ) -> Iterator[Mapping[str, Any]]:
        for page in self.iter_pages(params=params, page_size=page_size, context=context):
            for item in page.items:
                yield item

    def _normalize_page_payload(self, payload: Any, page_key: str | None, next_key: str | None) -> Page:
        if isinstance(payload, Page):
            return payload

        items: list[Any] = []
        next_cursor: str | None = None

        if isinstance(payload, Mapping):
            if page_key is None:
                items = list(payload.values()) if payload else []
            else:
                raw_items = payload.get(page_key, [])
                items = list(raw_items) if raw_items is not None else []
            if next_key is not None:
                next_cursor_val = payload.get(next_key)
                if isinstance(next_cursor_val, (str, type(None))):
                    next_cursor = next_cursor_val
                else:
                    next_cursor = str(next_cursor_val) if next_cursor_val is not None else None
            raw_payload: Any | None = payload
        else:
            items = list(payload) if payload is not None else []
            raw_payload = payload

        return Page(items=items, next_cursor=next_cursor, raw=raw_payload)
