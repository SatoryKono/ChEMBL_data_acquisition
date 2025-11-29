from __future__ import annotations

from typing import Any, Iterator, Mapping

from bioetl.clients.base import PaginationParams
from bioetl.clients.base.pagination import PaginationStrategyABC
from bioetl.clients.providers._provider_template import (
    ProviderDataClientImpl,
    RouteConfig,
)


class DummyHttpClient:
    def __init__(self) -> None:
        self.requests: list[tuple[str, Mapping[str, Any] | None]] = []

    def get_json(
        self,
        url: str,
        *,
        params: Mapping[str, Any] | None = None,
        context: Mapping[str, Any] | None = None,
    ) -> Any:
        self.requests.append((url, params))
        return {"id": "123", "context": context}

    def paginate_json(
        self,
        url: str,
        *,
        params: Mapping[str, Any] | None = None,
        page_key: str | None = None,
        next_key: str | None = None,
        page_param: str | None = None,
        context: Mapping[str, Any] | None = None,
    ) -> Iterator[Any]:
        yield {page_key or "items": [{"id": 1}, {"id": 2}], next_key or "next": None}

    def close(self) -> None:  # pragma: no cover - not used in these tests
        ...


class DummyNormalizer:
    def normalize(self, record: Mapping[str, Any]) -> Mapping[str, Any]:
        return {**record, "normalized": True}


class StubPaginationStrategy(PaginationStrategyABC):
    def __init__(self) -> None:
        self.calls: list[dict[str, Any]] = []

    def iter_pages(
        self,
        initial_response_or_args: Any,
        transport: DummyHttpClient,
        *,
        endpoint: str,
        params: Mapping[str, Any] | None = None,
        page_key: str | None = None,
        next_key: str | None = None,
        page_param: str | None = None,
    ) -> Iterator[Any]:
        self.calls.append(
            {
                "endpoint": endpoint,
                "params": params,
                "page_key": page_key,
                "next_key": next_key,
                "page_param": page_param,
            }
        )
        yield {page_key or "items": [{"id": 10}], next_key or "next": "cursor-1"}
        yield {page_key or "items": [{"id": 20}], next_key or "next": None}


def test_fetch_one_resolves_route_and_normalizes() -> None:
    http = DummyHttpClient()
    provider = ProviderDataClientImpl(
        http,
        routes=[RouteConfig(name="fetch", path="/entities/{value}", query_param="ref")],
        normalizer=DummyNormalizer(),
    )

    result = list(provider.fetch_one("CHEMBL1", params={"extra": "ok"}, context={"trace": "x"}))

    assert result == [{"id": "123", "context": {"trace": "x"}, "normalized": True}]
    assert http.requests == [("/entities/CHEMBL1", {"ref": "CHEMBL1", "extra": "ok"})]


def test_iter_pages_uses_strategy_and_pagination_defaults() -> None:
    strategy = StubPaginationStrategy()
    provider = ProviderDataClientImpl(
        DummyHttpClient(),
        routes=[RouteConfig(name="fetch", path="/entities/{value}")],
        default_pagination=PaginationParams(page_key="results", next_key="next", page_param="page"),
        pagination_strategy=strategy,
    )

    pages = list(provider.iter_pages(params={"_path": "/items", "foo": "bar"}))

    assert [p.items for p in pages] == [[{"id": 10}], [{"id": 20}]]
    assert pages[0].next_cursor == "cursor-1"
    assert strategy.calls == [
        {
            "endpoint": "/items",
            "params": {"foo": "bar"},
            "page_key": "results",
            "next_key": "next",
            "page_param": "page",
        }
    ]


def test_fetch_many_flattens_pages_with_normalizer() -> None:
    provider = ProviderDataClientImpl(
        DummyHttpClient(),
        routes=[RouteConfig(name="fetch", path="/entities/{value}")],
    )

    items = list(provider.fetch_many(params={"_path": "/items"}))

    assert items == [{"id": 1}, {"id": 2}]


def test_configure_updates_transport_and_pagination() -> None:
    http1 = DummyHttpClient()
    http2 = DummyHttpClient()
    provider = ProviderDataClientImpl(
        http1,
        routes=[RouteConfig(name="fetch", path="/entities/{value}")],
    )

    provider.configure(transport=http2, pagination=PaginationParams(page_key="data", next_key="cursor"))

    pages = list(provider.iter_pages(params={"_path": "/items"}))
    assert [p.items for p in pages] == [[{"id": 1}, {"id": 2}]]
    assert provider._http is http2
