from __future__ import annotations

from typing import Any, Mapping

from bioetl.clients.base import LoggingTransportAdapter, PaginationParams
from bioetl.clients.enricher_base import OptionsAwareApiClientImpl
from bioetl.clients.factories import ClientFactory
from bioetl.clients.factories.helpers import build_normalizer
from bioetl.clients.providers.openalex import OpenAlexDataClientImpl, OpenAlexNormalizerImpl


def _build_adapter(config: Mapping[str, Any]) -> LoggingTransportAdapter:
    adapter = config.get("adapter")
    if adapter is None:
        return LoggingTransportAdapter()
    if not isinstance(adapter, LoggingTransportAdapter):
        if callable(adapter):
            built = adapter()
            if isinstance(built, LoggingTransportAdapter):
                return built
        raise ValueError("Adapter has invalid type")
    return adapter


def _build_base_client(config: Mapping[str, Any]) -> Any:
    base_client = config.get("base_client")
    if callable(base_client):
        base_client = base_client()
    if base_client is None:
        raise ValueError("Base client is not configured")
    for attr in ("request", "get_json", "paginate_json", "close"):
        if not hasattr(base_client, attr):
            raise ValueError("Base client has invalid type")
    return base_client


class OpenAlexClientFactoryImpl(ClientFactory):
    """Фабрика клиентов OpenAlex."""

    def create(self, config: Mapping[str, Any]) -> OpenAlexDataClientImpl:
        base_client = _build_base_client(config)
        adapter = _build_adapter(config)
        transport = OptionsAwareApiClientImpl(base_client=base_client, adapter=adapter)

        normalizer_cfg = config.get("normalizer")
        normalizer = (
            build_normalizer({"normalizer": normalizer_cfg})
            if normalizer_cfg is not None
            else OpenAlexNormalizerImpl()
        )

        pagination = config.get("pagination")
        pagination_params = pagination if isinstance(pagination, PaginationParams) else None

        return OpenAlexDataClientImpl(
            http=transport,
            routes=config.get("routes") or OpenAlexDataClientImpl.ROUTES,
            default_pagination=pagination_params or OpenAlexDataClientImpl.DEFAULT_PAGINATION,
            pagination_strategy=config.get("pagination_strategy"),
            normalizer=normalizer,
            logger=config.get("logger"),
            options=config.get("options"),
        )
