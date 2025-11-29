from __future__ import annotations

from typing import Any, Iterable, Iterator, Mapping

from bioetl.clients.base import PaginationParams, RequestContext
from bioetl.clients.base.pagination import PaginationStrategyABC
from bioetl.clients.providers._provider_template import (
    ProviderDataClientImpl,
    RouteConfig,
)

from .openalex_normalizer_impl import OpenAlexNormalizerImpl


class OpenAlexDataClientImpl(ProviderDataClientImpl):
    """Конкретная реализация провайдера для OpenAlex."""

    ROUTES: tuple[RouteConfig, ...] = (
        RouteConfig("fetch", "/works/{value}"),
        RouteConfig("search", "/works", query_param="search"),
    )

    DEFAULT_PAGINATION = PaginationParams(page_key="results", next_key="meta.next", page_param="page")

    def __init__(
        self,
        http: Any,
        *,
        routes: Iterable[RouteConfig] | None = None,
        default_pagination: PaginationParams | None = None,
        pagination_strategy: PaginationStrategyABC | None = None,
        normalizer: OpenAlexNormalizerImpl | None = None,
        logger: Any | None = None,
        options: Mapping[str, Any] | None = None,
    ) -> None:
        super().__init__(
            http,
            routes=tuple(routes) if routes is not None else self.ROUTES,
            default_pagination=default_pagination or self.DEFAULT_PAGINATION,
            pagination_strategy=pagination_strategy,
            normalizer=normalizer or OpenAlexNormalizerImpl(),
            logger=logger,
            options=options,
        )

    def _iterate_pages_impl(
        self,
        params: Mapping[str, Any] | None,
        pagination: PaginationParams,
        context: RequestContext | None = None,
    ) -> Iterator[Any]:
        params = params or {}
        default_path = self._routes.get("search", RouteConfig("search", "/")).path
        path = params.get("_path") or default_path
        query_params = {k: v for k, v in params.items() if k != "_path"}
        yield from super()._iterate_pages_impl({"_path": path, **query_params}, pagination, context)

    def _normalize_page_payload(
        self, payload: Any, page_key: str | None, next_key: str | None
    ) -> Any:
        if isinstance(payload, Mapping) and next_key and "." in next_key:
            nested = payload
            for part in next_key.split("."):
                if isinstance(nested, Mapping):
                    nested = nested.get(part)
                else:
                    nested = None
                    break
            flattened_payload = dict(payload)
            flattened_payload[next_key] = nested
            return super()._normalize_page_payload(flattened_payload, page_key, next_key)

        return super()._normalize_page_payload(payload, page_key, next_key)

    def search(
        self,
        query: str,
        *,
        params: Mapping[str, Any] | None = None,
        page_size: int | None = None,
        pagination: PaginationParams | None = None,
        context: RequestContext | None = None,
    ) -> Iterator[Mapping[str, Any]]:
        path, search_params = self._resolve_route("search", value=query, params=params)
        merged_params: dict[str, Any] = {"_path": path}
        if search_params:
            merged_params.update(search_params)
        return self.fetch_many(
            params=merged_params,
            page_size=page_size,
            pagination=pagination,
            context=context,
        )

    # Deprecated aliases to keep backward compatibility.
    def fetch(
        self,
        ref: str | Mapping[str, Any],
        *,
        params: Mapping[str, Any] | None = None,
        context: RequestContext | None = None,
    ) -> Iterator[Mapping[str, Any]]:
        return self.fetch_one(ref, params=params, context=context)

    def fetch_search(
        self,
        query: str,
        *,
        params: Mapping[str, Any] | None = None,
        page_size: int | None = None,
        pagination: PaginationParams | None = None,
        context: RequestContext | None = None,
    ) -> Iterator[Mapping[str, Any]]:
        return self.search(
            query,
            params=params,
            page_size=page_size,
            pagination=pagination,
            context=context,
        )
