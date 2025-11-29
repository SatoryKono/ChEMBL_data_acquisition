from __future__ import annotations

from typing import Any, Mapping

from bioetl.clients.base.pagination import PaginationStrategyABC, TransportPaginationStrategyImpl
from bioetl.clients.base.normalizers import INormalizer
from bioetl.clients.base.interfaces import PaginationParams, RequestContext
from bioetl.clients.providers._provider_template import (
    ProviderDataClientImpl,
    RouteConfig,
)
from bioetl.clients.providers.semantic_scholar.semantic_scholar_normalizer_impl import (
    SemanticScholarNormalizerImpl,
)
from bioetl.clients.base.http import BaseHttpClientABC


class SemanticScholarDataClientImpl(ProviderDataClientImpl):
    """Клиент Semantic Scholar, совместимый с интерфейсом провайдера данных."""

    ROUTES = [
        RouteConfig("fetch", "/paper/{value}"),
        RouteConfig("search", "/paper/search", query_param="query"),
    ]

    def __init__(
        self,
        http: BaseHttpClientABC,
        *,
        normalizer: INormalizer | None = None,
        default_pagination: PaginationParams | None = None,
        pagination_strategy: PaginationStrategyABC | None = None,
        logger: Any | None = None,
        options: Mapping[str, Any] | None = None,
    ) -> None:
        pagination = default_pagination or PaginationParams(
            page_key="data", next_key="next", page_param="offset"
        )
        super().__init__(
            http,
            routes=self.ROUTES,
            default_pagination=pagination,
            pagination_strategy=pagination_strategy or TransportPaginationStrategyImpl(),
            normalizer=normalizer or SemanticScholarNormalizerImpl(),
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
        route = self._routes["search"]
        params_with_query = dict(params or {})
        if route.query_param:
            params_with_query = {route.query_param: query, **params_with_query}
        merged_params: Mapping[str, Any] = {"_path": route.path, **params_with_query}
        return self.fetch_many(
            params=merged_params,
            page_size=page_size,
            pagination=pagination,
            context=context,
        )

    def title_search(self, title: str, **kwargs: Any) -> Any:
        """Alias для поиска по заголовку."""

        return self.search(title, **kwargs)
