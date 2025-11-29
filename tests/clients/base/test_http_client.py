from __future__ import annotations

from unittest.mock import Mock

import pytest

from src.bioetl.clients.enricher_base import OptionsAwareApiClientImpl
from src.bioetl.clients.base.interfaces import (
    BaseApiClient,
    LoggingTransportAdapter,
    RequestContext,
    TransportError,
)


@pytest.fixture()
def base_client() -> Mock:
    return Mock(spec=BaseApiClient)


@pytest.fixture()
def adapter() -> LoggingTransportAdapter:
    return LoggingTransportAdapter()


@pytest.fixture()
def client(base_client: Mock, adapter: LoggingTransportAdapter) -> OptionsAwareApiClientImpl:
    return OptionsAwareApiClientImpl(base_client=base_client, adapter=adapter)


def test_request_proxies_arguments(client: OptionsAwareApiClientImpl, base_client: Mock) -> None:
    expected = {"ok": True}
    base_client.request.return_value = expected
    context = RequestContext(options={"trace_id": "123"})

    result = client.request(
        "GET",
        "http://example.com",
        params={"limit": 1},
        json=None,
        headers={"X-Test": "1"},
        context=context,
    )

    assert result is expected
    base_client.request.assert_called_once_with(
        "GET",
        "http://example.com",
        params={"limit": 1},
        json=None,
        headers={"X-Test": "1"},
        context=context,
    )
    assert client._adapter.current_context is context


def test_get_json_proxies_arguments(client: OptionsAwareApiClientImpl, base_client: Mock) -> None:
    payload = {"data": [1, 2, 3]}
    base_client.get_json.return_value = payload
    context = RequestContext(options={})

    result = client.get_json("http://example.com", params={"q": "x"}, context=context)

    assert result is payload
    base_client.get_json.assert_called_once_with(
        "http://example.com", params={"q": "x"}, context=context
    )


def test_paginate_json_proxies_arguments(client: OptionsAwareApiClientImpl, base_client: Mock) -> None:
    base_client.paginate_json.return_value = iter([[1], [2]])
    context = RequestContext(options=None)

    result = client.paginate_json(
        "http://example.com",
        params={"page": 1},
        page_key="items",
        next_key="next",
        page_param="page",
        context=context,
    )

    assert list(result) == [[1], [2]]
    base_client.paginate_json.assert_called_once_with(
        "http://example.com",
        params={"page": 1},
        page_key="items",
        next_key="next",
        page_param="page",
        context=context,
    )


def test_close_proxies(client: OptionsAwareApiClientImpl, base_client: Mock) -> None:
    client.close()
    base_client.close.assert_called_once_with()


def test_exceptions_are_wrapped(client: OptionsAwareApiClientImpl, base_client: Mock) -> None:
    base_client.request.side_effect = ValueError("boom")

    with pytest.raises(TransportError) as excinfo:
        client.request("GET", "http://example.com")

    assert "boom" in str(excinfo.value)
