from __future__ import annotations

from typing import Any, Iterator, Mapping

import pytest

from bioetl.clients.base.http import BaseHttpClientABC
from bioetl.clients.base.interfaces import RequestContext
from bioetl.clients.base.pagination import TransportPaginationStrategyImpl
from bioetl.clients.factories.crossref_client_factory_impl import CrossrefClientFactoryImpl
from bioetl.clients.providers.crossref import CrossrefPaginationStrategyImpl

pytestmark = pytest.mark.unit


class _FakeTransport(BaseHttpClientABC):
    def __init__(self, responses: Mapping[tuple[str, tuple[tuple[str, Any], ...] | None], Any]) -> None:
        self.responses = responses
        self.paginated_payloads: list[Any] = []

    def request(
        self,
        method: str,
        url: str,
        *,
        params: Mapping[str, Any] | None = None,
        json: Any | None = None,
        headers: Mapping[str, str] | None = None,
        context: RequestContext | None = None,
    ) -> Any:
        raise NotImplementedError

    def get_json(
        self,
        url: str,
        *,
        params: Mapping[str, Any] | None = None,
        context: RequestContext | None = None,
    ) -> Any:
        key = (url, tuple(sorted(params.items())) if params else None)
        return self.responses[key]

    def paginate_json(
        self,
        url: str,
        *,
        params: Mapping[str, Any] | None = None,
        page_key: str | None = None,
        next_key: str | None = None,
        page_param: str | None = None,
        context: RequestContext | None = None,
    ) -> Iterator[Any]:
        yield from self.paginated_payloads

    def close(self) -> None:  # pragma: no cover - nothing to cleanup
        return None


def _factory(transport: _FakeTransport, pagination: Mapping[str, Any] | None = None) -> CrossrefDataClientImpl:
    factory = CrossrefClientFactoryImpl()
    return factory.create({"transport": transport, "pagination": pagination or {}})


def test_fetch_one() -> None:
    doi = "10.1000/xyz"
    transport = _FakeTransport(
        {("/works/10.1000/xyz", None): {"message": {"DOI": doi, "title": ["Sample"]}}}
    )
    client = _factory(transport)

    result = list(client.fetch_one(doi))

    assert result[0]["id"] == f"doi:{doi}"
    assert result[0]["title"] == "Sample"


def test_fetch_many_uses_search_route() -> None:
    transport = _FakeTransport(
        {
            ("/works", (("query", "chemistry"),)): {
                "message": {"items": [{"DOI": "10.1", "title": ["A"]}], "next-cursor": None}
            }
        }
    )
    client = _factory(
        transport,
        pagination={"page_key": "message", "items_key": "items", "next_key": "next-cursor", "page_param": "cursor"},
    )

    result = list(client.search("chemistry"))

    assert len(result) == 1
    assert result[0]["id"] == "doi:10.1"


def test_iter_pages_with_message_items() -> None:
    transport = _FakeTransport(
        {
            ("/works", (("cursor", "*"),)): {
                "message": {
                    "items": [{"DOI": "10.2", "title": ["First"]}],
                    "next-cursor": "abc",
                }
            },
            ("/works", (("cursor", "abc"),)): {
                "message": {
                    "items": [{"DOI": "10.3", "title": ["Second"]}],
                    "next-cursor": None,
                }
            },
        }
    )
    client = _factory(
        transport,
        pagination={"page_key": "message", "items_key": "items", "next_key": "next-cursor", "page_param": "cursor"},
    )

    pages = list(client.iter_pages(params={"_path": "/works", "cursor": "*"}))

    assert len(pages) == 2
    assert pages[0].items[0]["DOI"] == "10.2"
    assert pages[1].items[0]["DOI"] == "10.3"


def test_pagination_strategy_selection() -> None:
    transport = _FakeTransport({("/works", None): {"message": {"items": [], "next-cursor": None}}})
    client = _factory(
        transport,
        pagination={"page_key": "message", "items_key": "items", "next_key": "next-cursor", "page_param": "cursor"},
    )

    assert isinstance(client.pagination_strategy, CrossrefPaginationStrategyImpl)

    client_default = _factory(transport, pagination={"page_key": "items"})

    assert isinstance(client_default.pagination_strategy, TransportPaginationStrategyImpl)


def test_normalizer_fields() -> None:
    record = {
        "message": {
            "DOI": "10.4000/demo",
            "title": ["Normalized Title"],
            "author": [
                {"given": "John", "family": "Doe", "affiliation": [{"name": "Uni"}]},
                {"family": "Solo"},
            ],
            "container-title": ["Journal of Tests"],
            "published-print": {"date-parts": [[2020, 5, 2]]},
            "references-count": 4,
        }
    }
    transport = _FakeTransport({("/works/anything", None): record})
    client = _factory(transport)

    normalized = list(client.fetch_one("anything"))[0]

    assert normalized["id"] == "doi:10.4000/demo"
    assert normalized["title"] == "Normalized Title"
    assert normalized["journal"] == "Journal of Tests"
    assert normalized["published_date"] == "2020-05-02"
    assert normalized["references_count"] == 4
    assert normalized["authors"][0]["name"] == "John Doe"
    assert normalized["authors"][0]["affiliations"] == ["Uni"]
