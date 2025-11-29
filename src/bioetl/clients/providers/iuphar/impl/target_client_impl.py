from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Iterator, Mapping

from bioetl.clients.base import (
    BaseHttpClientABC,
    DataProviderProtocol,
    PaginationParams,
    PaginationStrategyABC,
    RequestContext,
    TransportPaginationStrategyImpl,
)
from bioetl.clients.providers.base_provider import PagedDataProviderABC
from bioetl.clients.providers.iuphar.normalization import ProviderNormalizer


@dataclass(slots=True)
class RouteConfig:
    """Маршрут провайдера IUPHAR/targets (контракт RouteConfig).

    Ищите impl в ``providers/iuphar/impl`` и фабрику по умолчанию
    в ``clients/target/factories.py``.
    """

    name: str
    path: str
    query_param: str | None = None


class TargetIupharDataProviderImpl(PagedDataProviderABC, DataProviderProtocol):
    """Провайдер IUPHAR для сущности target (контракт DataProviderProtocol).

    Реализация располагается в ``providers/iuphar/impl``; фабричный метод
    ``default_target_iuphar`` доступен в ``clients/target/factories.py``.
    """

    ROUTES: tuple[RouteConfig, RouteConfig] = (
        RouteConfig(name="fetch", path="/targets/{identifier}"),
        RouteConfig(name="search", path="/targets"),
    )

    DEFAULT_PAGINATION = PaginationParams(page_key=None, next_key="next", page_param="page")

    def __init__(
        self,
        http: BaseHttpClientABC,
        *,
        normalizer: ProviderNormalizer | None = None,
        pagination_strategy: PaginationStrategyABC | None = None,
        logger: Any | None = None,
        options: Mapping[str, Any] | None = None,
    ) -> None:
        self._routes = {route.name: route for route in self.ROUTES}
        self._pagination_strategy = pagination_strategy or TransportPaginationStrategyImpl()
        self._normalizer = normalizer or ProviderNormalizer()
        merged_options: dict[str, Any] = {
            "page_size": self.DEFAULT_PAGINATION.page_size,
            "page_key": self.DEFAULT_PAGINATION.page_key,
            "next_key": self.DEFAULT_PAGINATION.next_key,
            "page_param": self.DEFAULT_PAGINATION.page_param,
        }
        if options:
            merged_options.update(options)
        super().__init__(http, logger=logger, options=merged_options)

    def _iterate_pages_impl(
        self,
        params: Mapping[str, Any] | None,
        pagination: PaginationParams,
        context: RequestContext | None = None,
    ) -> Iterator[Any]:
        params = params or {}
        path = params.get("_path") or self._routes["search"].path
        query_params = {key: value for key, value in params.items() if key != "_path"}
        if pagination.page_size is not None and "pageSize" not in query_params:
            query_params["pageSize"] = pagination.page_size

        yield from self._pagination_strategy.iter_pages(
            None,
            self._http,
            endpoint=path,
            params=query_params or None,
            page_key=pagination.page_key,
            next_key=pagination.next_key,
            page_param=pagination.page_param,
        )

    def iter_pages(
        self,
        *,
        params: Mapping[str, Any] | None = None,
        page_size: int | None = None,
        context: RequestContext | None = None,
    ) -> Iterator[Any]:
        return super().iter_pages(params=params, page_size=page_size, context=context)

    def fetch_many(
        self,
        *,
        params: Mapping[str, Any] | None = None,
        page_size: int | None = None,
        context: RequestContext | None = None,
    ) -> Iterator[Mapping[str, Any]]:
        def generator() -> Iterator[Mapping[str, Any]]:
            for page in self.iter_pages(params=params, page_size=page_size, context=context):
                for item in page.items:
                    yield self._normalizer.normalize(item)

        return self._wrap_iterator(generator(), context=context)

    def fetch_one(
        self,
        identifier: str | Mapping[str, Any],
        *,
        params: Mapping[str, Any] | None = None,
        context: RequestContext | None = None,
    ) -> Iterator[Mapping[str, Any]]:
        route = self._routes["fetch"]
        merged_params = {**(params or {})}
        if route.query_param:
            merged_params = {route.query_param: identifier, **merged_params}
        path = route.path.format(identifier=identifier)

        def generator() -> Iterator[Mapping[str, Any]]:
            payload = self._http.get_json(path, params=merged_params or None, context=context)
            if payload is None:
                return
            if isinstance(payload, Mapping):
                yield self._normalizer.normalize(payload)
            else:
                yield self._normalizer.normalize({"value": payload})

        return self._wrap_iterator(generator(), context=context)

    def search(
        self,
        query: str,
        *,
        params: Mapping[str, Any] | None = None,
        page_size: int | None = None,
        context: RequestContext | None = None,
    ) -> Iterator[Any]:
        search_route = self._routes.get("search")
        effective_params = {"_path": search_route.path if search_route else "/targets", **(params or {})}
        if "name" not in effective_params:
            effective_params["name"] = query
        return self.iter_pages(params=effective_params, page_size=page_size, context=context)

    def metadata(self) -> Mapping[str, Any]:
        return {"provider": "iuphar", "entity": "target", **super().metadata()}

    def close(self) -> None:
        super().close()
