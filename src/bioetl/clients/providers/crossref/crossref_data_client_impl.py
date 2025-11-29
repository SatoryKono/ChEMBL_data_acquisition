from __future__ import annotations

from typing import Any, Mapping, Sequence

from bioetl.clients.base.interfaces import PaginationParams, RequestContext
from bioetl.clients.base.pagination import PaginationStrategyABC, TransportPaginationStrategyImpl
from bioetl.clients.providers._provider_template import ProviderDataClientImpl, RouteConfig
from bioetl.clients.providers.crossref.crossref_normalizer_impl import CrossrefNormalizerImpl
from bioetl.clients.providers.crossref.crossref_pagination_strategy_impl import (
    CrossrefPaginationStrategyImpl,
)


class CrossrefDataClientImpl(ProviderDataClientImpl):
    """Клиент для провайдера Crossref."""

    ROUTES: Sequence[RouteConfig] = (
        RouteConfig("fetch", "/works/{value}"),
        RouteConfig("search", "/works", query_param="query"),
    )

    DEFAULT_HEADERS: Mapping[str, str] = {}
    OPTIONAL_CONFIG: Mapping[str, Any] = {}

    def __init__(
        self,
        http: Any,
        *,
        routes: Sequence[RouteConfig] | None = None,
        default_pagination: PaginationParams | None = None,
        pagination_strategy: PaginationStrategyABC | None = None,
        normalizer: CrossrefNormalizerImpl | None = None,
        logger: Any | None = None,
        options: Mapping[str, Any] | None = None,
    ) -> None:
        super().__init__(
            http,
            routes=routes or self.ROUTES,
            default_pagination=default_pagination,
            pagination_strategy=pagination_strategy or TransportPaginationStrategyImpl(),
            normalizer=normalizer or CrossrefNormalizerImpl(),
            logger=logger,
            options=options,
        )

    def search(
        self,
        query: str,
        *,
        params: Mapping[str, Any] | None = None,
        page_size: int | None = None,
        pagination: PaginationParams | None = None,
        context: RequestContext | None = None,
    ) -> Any:
        path, params_with_value = self._resolve_route("search", value=query, params=params)
        combined_params = {"_path": path, **(params_with_value or {})}
        return self.fetch_many(
            params=combined_params,
            page_size=page_size,
            pagination=pagination,
            context=context,
        )

    @staticmethod
    def build_pagination_strategy(config: Mapping[str, Any] | None) -> PaginationStrategyABC:
        if not config:
            return TransportPaginationStrategyImpl()

        if config.get("page_key") == "message":
            return CrossrefPaginationStrategyImpl(
                page_key=config.get("page_key", "message"),
                items_key=config.get("items_key", "items"),
                next_key=config.get("next_key", "next-cursor"),
                page_param=config.get("page_param", "cursor"),
            )

        return TransportPaginationStrategyImpl()
