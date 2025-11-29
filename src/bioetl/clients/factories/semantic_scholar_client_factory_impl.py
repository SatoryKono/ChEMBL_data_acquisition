from __future__ import annotations

from typing import Any, Mapping

from bioetl.clients.base.interfaces import PaginationParams
from bioetl.clients.base.pagination import PaginationStrategyABC
from bioetl.clients.factories import ClientFactory, build_normalizer, build_transport
from bioetl.clients.providers.semantic_scholar import (
    SemanticScholarDataClientImpl,
    SemanticScholarNormalizerImpl,
)


class SemanticScholarClientFactoryImpl(ClientFactory):
    """Factory creating configured SemanticScholar data client instances."""

    def create(self, config: Mapping[str, Any]) -> SemanticScholarDataClientImpl:
        transport = build_transport(config)
        normalizer = config.get("normalizer")
        if normalizer is None:
            normalizer_instance = SemanticScholarNormalizerImpl()
        else:
            normalizer_instance = build_normalizer(config)

        pagination = config.get("pagination")
        if isinstance(pagination, PaginationParams):
            pagination_params = pagination
        elif isinstance(pagination, Mapping):
            pagination_params = PaginationParams(**pagination)
        else:
            pagination_params = None

        pagination_strategy = config.get("pagination_strategy")
        pagination_strategy_instance: PaginationStrategyABC | None
        if isinstance(pagination_strategy, PaginationStrategyABC):
            pagination_strategy_instance = pagination_strategy
        elif callable(pagination_strategy):
            candidate = pagination_strategy()
            pagination_strategy_instance = (
                candidate if isinstance(candidate, PaginationStrategyABC) else None
            )
        else:
            pagination_strategy_instance = None

        return SemanticScholarDataClientImpl(
            transport,
            normalizer=normalizer_instance,
            default_pagination=pagination_params,
            pagination_strategy=pagination_strategy_instance,
            logger=config.get("logger"),
            options=config.get("options"),
        )
