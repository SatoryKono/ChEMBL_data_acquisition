from __future__ import annotations

from typing import Any, Iterator, Mapping

from bioetl.clients.base.pagination import PaginationStrategyABC
from bioetl.clients.providers.semantic_scholar import (
    SemanticScholarDataClientImpl,
    SemanticScholarNormalizerImpl,
)


class DummyHttpClient:
    def __init__(self) -> None:
        self.calls: list[dict[str, Any]] = []

    def get_json(
        self,
        url: str,
        *,
        params: Mapping[str, Any] | None = None,
        context: Mapping[str, Any] | None = None,
    ) -> Mapping[str, Any]:
        self.calls.append({"url": url, "params": params, "context": context})
        return {
            "paperId": "S2:123",
            "title": "Example Paper",
            "authors": [{"name": "Author One", "authorId": "A1"}],
            "year": 2024,
            "abstract": "Abstract text",
            "citationCount": 10,
        }

    def paginate_json(
        self,
        url: str,
        *,
        params: Mapping[str, Any] | None = None,
        page_key: str | None = None,
        next_key: str | None = None,
        page_param: str | None = None,
        context: Mapping[str, Any] | None = None,
    ) -> Iterator[Mapping[str, Any]]:
        self.calls.append(
            {
                "url": url,
                "params": params,
                "page_key": page_key,
                "next_key": next_key,
                "page_param": page_param,
            }
        )
        yield {page_key or "data": [{"paperId": "S2:1"}], next_key or "next": "cursor"}
        yield {page_key or "data": [{"paperId": "S2:2"}], next_key or "next": None}

    def close(self) -> None:  # pragma: no cover - not used in these tests
        ...


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
        yield {page_key or "data": [{"paperId": "S2:10"}], next_key or "next": "n1"}
        yield {page_key or "data": [{"paperId": "S2:20"}], next_key or "next": None}


def test_fetch_one() -> None:
    http = DummyHttpClient()
    client = SemanticScholarDataClientImpl(http)

    result = list(client.fetch_one("S2:42", params={"fields": "title"}))

    assert result[0]["id"] == "S2:123"
    assert http.calls[0]["url"] == "/paper/S2:42"
    assert http.calls[0]["params"] == {"fields": "title"}


def test_fetch_many() -> None:
    http = DummyHttpClient()
    client = SemanticScholarDataClientImpl(http)

    items = list(client.fetch_many(params={"_path": "/paper/search", "query": "covid"}))

    assert [item["id"] for item in items] == ["S2:1", "S2:2"]
    assert http.calls[0]["url"] == "/paper/search"
    assert http.calls[0]["params"] == {"query": "covid"}


def test_iter_pages() -> None:
    strategy = StubPaginationStrategy()
    http = DummyHttpClient()
    client = SemanticScholarDataClientImpl(http, pagination_strategy=strategy)

    pages = list(client.iter_pages(params={"_path": "/paper/search", "query": "ml"}))

    assert [page.items for page in pages] == [[{"paperId": "S2:10"}], [{"paperId": "S2:20"}]]
    assert pages[0].next_cursor == "n1"
    assert strategy.calls == [
        {
            "endpoint": "/paper/search",
            "params": {"query": "ml"},
            "page_key": "data",
            "next_key": "next",
            "page_param": "offset",
        }
    ]


def test_title_search_alias() -> None:
    http = DummyHttpClient()
    client = SemanticScholarDataClientImpl(http)

    items = list(client.title_search("graph"))

    assert [item["id"] for item in items] == ["S2:1", "S2:2"]
    assert http.calls[0]["url"] == "/paper/search"
    assert http.calls[0]["params"] == {"query": "graph"}


def test_normalizer_fields() -> None:
    normalizer = SemanticScholarNormalizerImpl()
    record = {
        "paperId": "S2:99",
        "title": "Sample",
        "authors": [{"name": "Alice"}, {"name": "Bob", "authorId": "A2"}],
        "year": 2023,
        "abstract": "Text",
        "citationCount": 7,
    }

    normalized = normalizer.normalize(record)

    assert normalized == {
        "id": "S2:99",
        "title": "Sample",
        "authors": [
            {"name": "Alice", "authorId": None},
            {"name": "Bob", "authorId": "A2"},
        ],
        "year": 2023,
        "abstract": "Text",
        "citation_count": 7,
    }
