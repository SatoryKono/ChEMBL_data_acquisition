from __future__ import annotations

from typing import Any, Mapping, Sequence

from bioetl.clients.base import (
    BaseHttpClientABC,
    PaginationParams,
    RequestContext,
    TransportPaginationStrategyImpl,
)
from bioetl.clients.base.normalizers import INormalizer
from bioetl.clients.providers._provider_template import (
    ProviderDataClientImpl,
    RouteConfig,
)

from .uniprot_normalizer_impl import UniprotNormalizerImpl


class UniprotDataClientImpl(ProviderDataClientImpl):
    """Провайдер UniProt с поддержкой поиска и получения записей."""

    ROUTES: Sequence[RouteConfig] = (
        RouteConfig("fetch", "/uniprot/{value}.json"),
        RouteConfig("search", "/uniprot/search", query_param="query"),
    )

    def __init__(
        self,
        http: BaseHttpClientABC,
        *,
        pagination_strategy: TransportPaginationStrategyImpl | None = None,
        normalizer: INormalizer | None = None,
        options: Mapping[str, Any] | None = None,
    ) -> None:
        default_pagination = PaginationParams(
            page_key="results",
            next_key="next",
            page_param="cursor",
        )
        super().__init__(
            http,
            routes=self.ROUTES,
            default_pagination=default_pagination,
            pagination_strategy=pagination_strategy or TransportPaginationStrategyImpl(),
            normalizer=normalizer or UniprotNormalizerImpl(),
            options=options,
        )

    def fetch_one(
        self,
        ref: str | Mapping[str, Any],
        *,
        params: Mapping[str, Any] | None = None,
        context: RequestContext | None = None,
    ) -> Any:
        return super().fetch_one(ref, params=params, context=context)

    def fetch_many(
        self,
        *,
        params: Mapping[str, Any] | None = None,
        page_size: int | None = None,
        context: RequestContext | None = None,
    ) -> Any:
        return super().fetch_many(params=params, page_size=page_size, context=context)

    def iter_pages(
        self,
        *,
        params: Mapping[str, Any] | None = None,
        page_size: int | None = None,
        pagination: PaginationParams | None = None,
        context: RequestContext | None = None,
    ) -> Any:
        effective_pagination = pagination or PaginationParams(
            page_key="results",
            next_key="next",
            page_param="cursor",
        )
        return super().iter_pages(
            params=params,
            page_size=page_size,
            pagination=effective_pagination,
            context=context,
        )
