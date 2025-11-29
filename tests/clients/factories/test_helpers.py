from __future__ import annotations

from typing import Any, Iterator, Mapping

import pytest

from bioetl.clients.base.http import BaseHttpClientABC
from bioetl.clients.base.interfaces import RequestContext
from bioetl.clients.base.normalizers import INormalizer, IdentityNormalizerImpl
from bioetl.clients.factories.helpers import build_normalizer, build_transport


class DummyTransport(BaseHttpClientABC):
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
        return {
            "method": method,
            "url": url,
            "params": params,
            "json": json,
            "headers": headers,
            "context": context,
        }

    def get_json(
        self,
        url: str,
        *,
        params: Mapping[str, Any] | None = None,
        context: RequestContext | None = None,
    ) -> Any:
        return {"url": url, "params": params, "context": context}

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
        yield {"url": url, "params": params, "context": context}

    def close(self) -> None:
        return None


class DummyNormalizer(INormalizer):
    def normalize(self, record: Mapping[str, Any]) -> Mapping[str, Any]:
        return {**record, "normalized": True}


def test_build_transport_from_instance() -> None:
    transport = DummyTransport()

    assert build_transport({"transport": transport}) is transport


def test_build_transport_from_callable() -> None:
    transport = build_transport({"transport": DummyTransport})

    assert isinstance(transport, DummyTransport)


def test_build_transport_without_value_raises_value_error() -> None:
    with pytest.raises(ValueError):
        build_transport({})


def test_build_normalizer_defaults_to_identity() -> None:
    normalizer = build_normalizer({})

    assert isinstance(normalizer, IdentityNormalizerImpl)


def test_build_normalizer_from_instance() -> None:
    normalizer = DummyNormalizer()

    assert build_normalizer({"normalizer": normalizer}) is normalizer


def test_build_normalizer_from_callable() -> None:
    normalizer = build_normalizer({"normalizer": DummyNormalizer})

    assert isinstance(normalizer, DummyNormalizer)


def test_build_normalizer_with_invalid_type_raises_value_error() -> None:
    with pytest.raises(ValueError):
        build_normalizer({"normalizer": "unexpected"})
