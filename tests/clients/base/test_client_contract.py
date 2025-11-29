from __future__ import annotations

from collections.abc import Iterator, Mapping
from typing import Any

from bioetl.clients.base.client import (
    ClientRequest,
    ExternalDataClient,
    Page,
    Pagination,
    Record,
    RequestContext,
)


def test_client_request_defaults() -> None:
    request = ClientRequest(route="test")

    assert request.route == "test"
    assert request.params is None
    assert request.context is None
    assert request.pagination is None


def test_request_context_structure() -> None:
    context = RequestContext(trace_id="trace-1", options={"timeout": 5})

    assert hasattr(RequestContext, "__slots__")
    assert context.trace_id == "trace-1"
    assert context.options == {"timeout": 5}


class _DummyClient:
    closed = False

    def fetch_one(self, request: ClientRequest) -> Record | None:  # noqa: ARG002
        return {"id": 1}

    def fetch_many(self, request: ClientRequest) -> Iterator[Record]:  # noqa: ARG002
        yield {"id": 1}

    def iter_pages(self, request: ClientRequest) -> Iterator[Page]:  # noqa: ARG002
        yield Page(items=[{"id": 1}], next_cursor=None, raw={"meta": {"cursor": None}})

    def metadata(self) -> Mapping[str, Any]:
        return {"source": "dummy"}

    def close(self) -> None:
        self.closed = True


def test_external_data_client_protocol() -> None:
    client = _DummyClient()

    assert isinstance(client, ExternalDataClient)
    assert list(client.fetch_many(ClientRequest(route="r")))[0]["id"] == 1
    assert isinstance(next(client.iter_pages(ClientRequest(route="r"))), Page)


def test_pagination_structure() -> None:
    pagination = Pagination(page_size=10, cursor="abc")

    assert hasattr(Pagination, "__slots__")
    assert pagination.page_size == 10
    assert pagination.cursor == "abc"
