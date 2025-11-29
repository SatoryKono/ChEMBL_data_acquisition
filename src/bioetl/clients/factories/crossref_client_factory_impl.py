from __future__ import annotations

from typing import Any, Mapping

from bioetl.clients.base.interfaces import PaginationParams
from bioetl.clients.base.pagination import TransportPaginationStrategyImpl
from bioetl.clients.enricher_base import OptionsAwareApiClientImpl
from bioetl.clients.factories.client_factory import ClientFactory
from bioetl.clients.factories.helpers import build_normalizer, build_transport
from bioetl.clients.providers.crossref.crossref_data_client_impl import CrossrefDataClientImpl
from bioetl.clients.providers.crossref.crossref_normalizer_impl import CrossrefNormalizerImpl
from bioetl.clients.providers.crossref.crossref_pagination_strategy_impl import (
    CrossrefPaginationStrategyImpl,
)


class CrossrefClientFactoryImpl(ClientFactory):
    """Фабрика клиентов для провайдера Crossref."""

    def create(self, config: Mapping[str, Any]) -> CrossrefDataClientImpl:
        pagination_cfg = config.get("pagination") or {}
        pagination_params = PaginationParams(
            page_size=pagination_cfg.get("page_size"),
            page_key=pagination_cfg.get("page_key"),
            next_key=pagination_cfg.get("next_key"),
            page_param=pagination_cfg.get("page_param"),
        )

        if pagination_params.page_key == "message":
            pagination_strategy = CrossrefPaginationStrategyImpl(
                page_key=pagination_params.page_key or "message",
                items_key=pagination_cfg.get("items_key", "items"),
                next_key=pagination_params.next_key or "next-cursor",
                page_param=pagination_params.page_param or "cursor",
            )
        else:
            pagination_strategy = config.get("pagination_strategy") or TransportPaginationStrategyImpl()

        transport = config.get("http") or config.get("transport")
        if transport is None:
            base_client = config.get("base_client")
            adapter = config.get("adapter")
            if base_client and adapter:
                transport = OptionsAwareApiClientImpl(base_client=base_client, adapter=adapter)
            else:
                transport = build_transport(config)

        normalizer = config.get("normalizer") or build_normalizer(config)
        if not isinstance(normalizer, CrossrefNormalizerImpl):
            if callable(normalizer):  # type: ignore[unreachable]
                normalizer = normalizer()  # pragma: no cover - defensive
            else:
                normalizer = CrossrefNormalizerImpl()

        options = config.get("options")

        return CrossrefDataClientImpl(
            http=transport,
            routes=config.get("routes") or CrossrefDataClientImpl.ROUTES,
            default_pagination=pagination_params,
            pagination_strategy=pagination_strategy,
            normalizer=normalizer,
            logger=config.get("logger"),
            options=options,
        )
