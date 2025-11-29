from __future__ import annotations

from typing import Any, Iterator, Mapping
from unittest.mock import Mock

import pytest

from bioetl.clients.base.interfaces import DataProviderError, Page, PaginationParams, RequestContext
from bioetl.clients.providers.base_provider import BaseDataProviderABC, PagedDataProviderABC


class DummyHttpClient:
    def __init__(self, metadata: Mapping[str, Any] | None = None) -> None:
        self.metadata = metadata or {}
        self.closed = False

    def close(self) -> None:
        self.closed = True


class StubProvider(BaseDataProviderABC):
    def fetch_one(
        self,
        identifier: str | Mapping[str, Any],
        *,
        params: Mapping[str, Any] | None = None,
        context: RequestContext | None = None,
    ) -> Iterator[Mapping[str, Any]]:
        return iter([])

    def iter_pages(
        self,
        *,
        params: Mapping[str, Any] | None = None,
        page_size: int | None = None,
        context: RequestContext | None = None,
    ) -> Iterator[Page]:
        return iter([])

    def fetch_many(
        self,
        *,
        params: Mapping[str, Any] | None = None,
        page_size: int | None = None,
        context: RequestContext | None = None,
    ) -> Iterator[Mapping[str, Any]]:
        return iter([])


class DummyPagedProvider(PagedDataProviderABC):
    def __init__(self, http: DummyHttpClient, *, options: Mapping[str, Any] | None = None) -> None:
        super().__init__(http, options=options)

    def fetch_one(
        self,
        identifier: str | Mapping[str, Any],
        *,
        params: Mapping[str, Any] | None = None,
        context: RequestContext | None = None,
    ) -> Iterator[Mapping[str, Any]]:
        return iter([])

    def _iterate_pages_impl(
        self,
        params: Mapping[str, Any] | None,
        pagination: PaginationParams,
        context: RequestContext | None = None,
    ) -> Iterator[Any]:
        yield {pagination.page_key: [{"id": 1}], pagination.next_key: "next-cursor"}
        if params and params.get("fail"):
            raise ValueError("boom")
        yield Page(items=[{"id": 2}], next_cursor=None, raw=None)


class GeneratorError(Exception):
    pass


def test_wrap_callable_success_and_error_logging() -> None:
    http = DummyHttpClient()
    logger = Mock()
    provider = StubProvider(http=http, logger=logger)

    assert provider._wrap_callable(lambda: "ok") == "ok"

    with pytest.raises(DataProviderError) as excinfo:
        provider._wrap_callable(lambda: (_ for _ in ()).throw(ValueError("fail")))
    assert "fail" in str(excinfo.value)
    logger.exception.assert_called_once()


def test_wrap_iterator_wraps_errors() -> None:
    http = DummyHttpClient()
    logger = Mock()
    provider = StubProvider(http=http, logger=logger)

    def good_iter() -> Iterator[int]:
        yield 1
        yield 2

    assert list(provider._wrap_iterator(good_iter())) == [1, 2]

    def bad_iter() -> Iterator[int]:
        yield 1
        raise GeneratorError("iterator failed")

    with pytest.raises(DataProviderError) as excinfo:
        list(provider._wrap_iterator(bad_iter()))
    assert "iterator failed" in str(excinfo.value)
    logger.exception.assert_called_once()


def test_metadata_and_close_proxied() -> None:
    http = DummyHttpClient(metadata={"foo": "bar"})
    provider = StubProvider(http=http)

    assert provider.metadata() == {"foo": "bar"}
    provider.close()
    assert http.closed is True


def test_paged_provider_iter_pages_normalizes_and_wraps_errors() -> None:
    provider = DummyPagedProvider(DummyHttpClient())

    pages = list(provider.iter_pages())
    assert len(pages) == 2

    first, second = pages
    assert isinstance(first, Page)
    assert first.items == [{"id": 1}]
    assert first.next_cursor == "next-cursor"
    assert isinstance(second, Page)
    assert second.items == [{"id": 2}]

    with pytest.raises(DataProviderError):
        list(provider.iter_pages(params={"fail": True}))


def test_fetch_many_flattens_items() -> None:
    provider = DummyPagedProvider(DummyHttpClient())
    items = list(provider.fetch_many())
    assert items == [{"id": 1}, {"id": 2}]
