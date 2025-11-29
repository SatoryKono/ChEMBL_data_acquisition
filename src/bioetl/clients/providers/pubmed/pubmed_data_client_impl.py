from __future__ import annotations

from typing import Any, Iterator, Mapping

from bioetl.clients.base import (
    INormalizer,
    PaginationParams,
    RequestContext,
    TransportPaginationStrategyImpl,
)
from bioetl.clients.base.interfaces import Page
from bioetl.clients.providers._provider_template import (
    ProviderDataClientImpl,
    RouteConfig,
)


class PubMedDataClientImpl(ProviderDataClientImpl):
    """Провайдер PubMed с поддержкой поиска и выборки статей."""

    ROUTES = (
        RouteConfig(name="fetch", path="/efetch.fcgi", query_param="id"),
        RouteConfig(name="search", path="/esearch.fcgi", query_param="term"),
    )

    def __init__(
        self,
        http: Any,
        *,
        normalizer: INormalizer | None = None,
        logger: Any | None = None,
        options: Mapping[str, Any] | None = None,
    ) -> None:
        pagination = PaginationParams(
            page_size=None,
            page_key="esearchresult",
            next_key=None,
            page_param="retstart",
        )
        merged_options = {"items_key": "idlist"}
        if options:
            merged_options.update(options)
        super().__init__(
            http,
            routes=self.ROUTES,
            default_pagination=pagination,
            pagination_strategy=TransportPaginationStrategyImpl(),
            normalizer=normalizer,
            logger=logger,
            options=merged_options,
        )
        self._items_key = merged_options.get("items_key", "idlist")

    def _normalize_page_payload(self, payload: Any, page_key: str | None, next_key: str | None):
        if isinstance(payload, Mapping):
            items_source = payload.get(page_key) if page_key else payload
            next_cursor = payload.get(next_key) if next_key else None
            if isinstance(items_source, Mapping):
                idlist = items_source.get(self._items_key)
                if idlist is not None:
                    items: list[Any] = list(idlist)
                else:
                    items = list(items_source.values()) if items_source else []
            elif items_source is None:
                items = []
            else:
                items = list(items_source)
            return Page(items=items, next_cursor=next_cursor, raw=payload)
        return super()._normalize_page_payload(payload, page_key, next_key)

    def _iterate_pages_impl(
        self,
        params: Mapping[str, Any] | None,
        pagination: PaginationParams,
        context: RequestContext | None = None,
    ) -> Iterator[Any]:
        params = params or {}
        path = params.get("_path") or "/"
        query_params = {k: v for k, v in params.items() if k != "_path"}
        if pagination.page_size is not None and "retmax" not in query_params:
            query_params["retmax"] = pagination.page_size

        yield from self.pagination_strategy.iter_pages(
            None,
            self._http,
            endpoint=path,
            params=query_params or None,
            page_key=pagination.page_key,
            next_key=pagination.next_key,
            page_param=pagination.page_param,
        )

    def search(
        self,
        term: str,
        *,
        query_param: str | None = "term",
        params: Mapping[str, Any] | None = None,
        page_size: int | None = None,
        context: RequestContext | None = None,
    ) -> Iterator[Any]:
        route = self._routes["search"]
        query_params: dict[str, Any] = {**(params or {})}
        if query_param:
            query_params = {query_param: term, **query_params}
        return self.iter_pages(
            params={"_path": route.path, **query_params},
            page_size=page_size,
            context=context,
        )

    def search_by_title(self, title: str, **kwargs: Any) -> Iterator[Any]:
        return self.search(title, query_param="title", **kwargs)

    def fetch(
        self,
        identifier: str,
        *,
        params: Mapping[str, Any] | None = None,
        context: RequestContext | None = None,
    ) -> Iterator[Mapping[str, Any]]:
        return self.fetch_one(identifier, params=params, context=context)

    def fetch_by_pmid(
        self,
        pmid: str,
        **kwargs: Any,
    ) -> Iterator[Mapping[str, Any]]:
        return self.fetch(pmid, **kwargs)

    def fetch_many(
        self,
        *,
        params: Mapping[str, Any] | None = None,
        page_size: int | None = None,
        pagination: PaginationParams | None = None,
        context: RequestContext | None = None,
    ) -> Iterator[Mapping[str, Any]]:
        search_params = params or {}
        if pagination is None:
            pagination = self._compose_pagination(page_size=page_size, pagination=None)
        if pagination.page_size is not None and "retmax" not in search_params:
            search_params = {**search_params, "retmax": pagination.page_size}
        return super().fetch_many(
            params=search_params,
            page_size=page_size,
            pagination=pagination,
            context=context,
        )
