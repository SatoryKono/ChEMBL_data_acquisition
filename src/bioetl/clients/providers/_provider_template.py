from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Iterable, Iterator, Mapping, Protocol, Sequence

from bioetl.clients.base import (
    BaseHttpClientABC,
    INormalizer,
    IdentityNormalizerImpl,
    PaginationParams,
    PaginationStrategyABC,
    RequestContext,
    TransportPaginationStrategyImpl,
)
from bioetl.clients.providers import PagedDataProviderABC


class DataProviderProtocol(Protocol):
    """Протокол, описывающий интерфейс провайдеров данных."""

    def fetch_one(
        self,
        ref: str | Mapping[str, Any],
        *,
        params: Mapping[str, Any] | None = None,
        context: RequestContext | None = None,
    ) -> Iterator[Mapping[str, Any]]:
        ...

    def iter_pages(
        self,
        *,
        params: Mapping[str, Any] | None = None,
        page_size: int | None = None,
        context: RequestContext | None = None,
    ) -> Iterator[Any]:
        ...

    def fetch_many(
        self,
        *,
        params: Mapping[str, Any] | None = None,
        page_size: int | None = None,
        context: RequestContext | None = None,
    ) -> Iterator[Mapping[str, Any]]:
        ...


class ConfigurableProvider(Protocol):
    """Протокол конфигурируемых провайдеров."""

    def configure(
        self,
        *,
        transport: BaseHttpClientABC | None = None,
        pagination: PaginationParams | None = None,
        retries: Any | None = None,
    ) -> ConfigurableProvider:
        ...


@dataclass(slots=True)
class RouteConfig:
    """Описание маршрута для провайдера."""

    name: str
    path: str
    query_param: str | None = None


def normalize_payload(payload: Any, *, page_key: str | None) -> Iterator[Mapping[str, Any]]:
    """Итерирует записи из произвольного payload."""

    if payload is None:
        return iter(())

    if isinstance(payload, Mapping):
        if page_key is None:
            return iter([dict(payload)])
        items = payload.get(page_key) or []
        return _normalize_items(items)

    if isinstance(payload, Iterable) and not isinstance(payload, (str, bytes)):
        return _normalize_items(payload)

    return iter([payload])


def _normalize_items(items: Iterable[Any]) -> Iterator[Mapping[str, Any]]:
    for item in items:
        if isinstance(item, Mapping):
            yield dict(item)
        else:
            yield {"value": item}


class ProviderDataClientImpl(PagedDataProviderABC, DataProviderProtocol, ConfigurableProvider):
    """Шаблонный провайдер с поддержкой пагинации и нормализации."""

    def __init__(
        self,
        http: BaseHttpClientABC,
        *,
        routes: Sequence[RouteConfig],
        default_pagination: PaginationParams | None = None,
        pagination_strategy: PaginationStrategyABC | None = None,
        normalizer: INormalizer | None = None,
        logger: Any | None = None,
        options: Mapping[str, Any] | None = None,
    ) -> None:
        self._routes = {route.name: route for route in routes}
        self._default_pagination = default_pagination or PaginationParams()
        self.pagination_strategy = pagination_strategy or TransportPaginationStrategyImpl()
        self.normalizer = normalizer or IdentityNormalizerImpl()

        merged_options = {
            "page_size": self._default_pagination.page_size,
            "page_key": self._default_pagination.page_key,
            "next_key": self._default_pagination.next_key,
            "page_param": self._default_pagination.page_param,
        }
        if options:
            merged_options.update(options)

        super().__init__(http, logger=logger, options=merged_options)

    def _resolve_route(
        self,
        route_name: str,
        *,
        value: str | Mapping[str, Any],
        params: Mapping[str, Any] | None,
    ) -> tuple[str, Mapping[str, Any] | None]:
        route = self._routes[route_name]
        params_with_value = dict(params or {})
        if route.query_param:
            params_with_value = {route.query_param: value, **params_with_value}
        path = route.path.format(value=value)
        return path, params_with_value or None

    def _iterate_pages_impl(
        self,
        params: Mapping[str, Any] | None,
        pagination: PaginationParams,
        context: RequestContext | None = None,
    ) -> Iterator[Any]:
        params = params or {}
        path = params.get("_path") or "/"
        query_params = {k: v for k, v in params.items() if k != "_path"}

        if self.pagination_strategy:
            yield from self.pagination_strategy.iter_pages(
                None,
                self._http,
                endpoint=path,
                params=query_params or None,
                page_key=pagination.page_key,
                next_key=pagination.next_key,
                page_param=pagination.page_param,
            )
        else:
            yield from self._http.paginate_json(
                path,
                params=query_params or None,
                page_key=pagination.page_key,
                next_key=pagination.next_key,
                page_param=pagination.page_param,
            )

    def _compose_pagination(
        self,
        *,
        page_size: int | None = None,
        pagination: PaginationParams | None = None,
    ) -> PaginationParams:
        base_pagination = pagination or self._default_pagination
        return PaginationParams(
            page_size=page_size if page_size is not None else base_pagination.page_size,
            page_key=base_pagination.page_key,
            next_key=base_pagination.next_key,
            page_param=base_pagination.page_param,
        )

    def iter_pages(
        self,
        *,
        params: Mapping[str, Any] | None = None,
        page_size: int | None = None,
        pagination: PaginationParams | None = None,
        context: RequestContext | None = None,
    ) -> Iterator[Any]:
        effective_pagination = self._compose_pagination(page_size=page_size, pagination=pagination)

        def generator() -> Iterator[Any]:
            for raw_page in self._iterate_pages_impl(params, pagination=effective_pagination, context=context):
                yield self._normalize_page_payload(
                    raw_page,
                    page_key=effective_pagination.page_key,
                    next_key=effective_pagination.next_key,
                )

        return self._wrap_iterator(generator(), context=context)

    def fetch_many(
        self,
        *,
        params: Mapping[str, Any] | None = None,
        page_size: int | None = None,
        pagination: PaginationParams | None = None,
        context: RequestContext | None = None,
    ) -> Iterator[Mapping[str, Any]]:
        for page in self.iter_pages(
            params=params,
            page_size=page_size,
            pagination=pagination,
            context=context,
        ):
            for item in page.items:
                yield self.normalizer.normalize(item)

    def fetch_one(
        self,
        ref: str | Mapping[str, Any],
        *,
        params: Mapping[str, Any] | None = None,
        context: RequestContext | None = None,
    ) -> Iterator[Mapping[str, Any]]:
        path, params_with_value = self._resolve_route("fetch", value=ref, params=params)
        payload = self._http.get_json(path, params=params_with_value, context=context)
        for record in normalize_payload(payload, page_key=None):
            yield self.normalizer.normalize(record)

    def configure(
        self,
        *,
        transport: BaseHttpClientABC | None = None,
        pagination: PaginationParams | None = None,
        retries: Any | None = None,
    ) -> ProviderDataClientImpl:
        if transport:
            self._http = transport
        if pagination:
            self._default_pagination = pagination
            self._options.update(
                {
                    "page_size": pagination.page_size,
                    "page_key": pagination.page_key,
                    "next_key": pagination.next_key,
                    "page_param": pagination.page_param,
                }
            )
        return self
