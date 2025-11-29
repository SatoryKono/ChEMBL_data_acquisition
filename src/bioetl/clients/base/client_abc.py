from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Iterator, Mapping, MutableMapping, Protocol, TypedDict, runtime_checkable

from bioetl.clients.base.types import PaginationConfig
from bioetl.clients.config.models import ClientConfig


class Record(TypedDict, total=False):
    """Минимально типизированная запись."""


@dataclass(slots=True)
class RequestContext:
    """Контекст запроса, распространяемый по цепочке вызовов."""

    trace_id: str | None = None
    options: Mapping[str, Any] | None = None


@dataclass(slots=True)
class ClientRequest:
    """Унифицированное описание запроса к клиенту."""

    route: str
    params: Mapping[str, Any] | None = None
    context: RequestContext | None = None
    pagination: PaginationConfig | None = None


@dataclass(slots=True)
class Page:
    """Коллекция элементов и вспомогательной информации о продолжении."""

    items: list[Record]
    next_cursor: str | None = None
    raw: Any | None = None


@runtime_checkable
class ExternalDataClient(Protocol):
    """Контракт для клиентов внешних источников данных."""

    def fetch_one(self, request: ClientRequest) -> Record | None:
        ...

    def fetch_many(self, request: ClientRequest) -> Iterator[Record]:
        ...

    def iter_pages(self, request: ClientRequest) -> Iterator[Page]:
        ...

    def metadata(self) -> Mapping[str, Any]:
        ...

    def close(self) -> None:
        ...


class BaseExternalDataClient(ExternalDataClient):
    """Базовая реализация клиента с общими сценариями пагинации."""

    def __init__(self, *, base_url: str, routes: Mapping[str, str], default_params: Mapping[str, Any] | None = None,
                 pagination: PaginationConfig | None = None, transport: Any | None = None) -> None:
        self._base_url = base_url.rstrip("/")
        self._routes = dict(routes)
        self._default_params = dict(default_params or {})
        self._pagination = pagination
        self._transport = transport

    def _build_url(self, route_name: str, params: Mapping[str, Any] | None) -> tuple[str, MutableMapping[str, Any]]:
        try:
            path = self._routes[route_name]
        except KeyError as exc:
            raise KeyError(f"Route '{route_name}' is not configured") from exc
        path_params = dict(params or {})
        url = f"{self._base_url}{path.format(**path_params)}"
        return url, path_params

    def iter_pages(self, request: ClientRequest) -> Iterator[Page]:
        url, params = self._build_url(request.route, request.params)
        merged_params: MutableMapping[str, Any] = {**self._default_params, **params}
        pagination = request.pagination or self._pagination
        next_cursor: str | None = None

        while True:
            effective_params = dict(merged_params)
            if pagination and pagination.page_param and next_cursor:
                effective_params[pagination.page_param] = next_cursor
            raw = self._transport.get_json(url, params=effective_params, context=request.context)
            items = self._extract_items(raw, pagination)
            next_cursor = self._extract_next_cursor(raw, pagination)
            yield Page(items=items, next_cursor=next_cursor, raw=raw)
            if not next_cursor:
                break

    def fetch_many(self, request: ClientRequest) -> Iterator[Record]:
        for page in self.iter_pages(request):
            for item in page.items:
                yield item

    def fetch_one(self, request: ClientRequest) -> Record | None:
        for page in self.iter_pages(request):
            if page.items:
                return page.items[0]
        return None

    def metadata(self) -> Mapping[str, Any]:
        return {"base_url": self._base_url, "routes": self._routes}

    def close(self) -> None:  # pragma: no cover - транспорт может не требовать закрытия
        if hasattr(self._transport, "close"):
            self._transport.close()

    @staticmethod
    def _extract_items(raw: Any, pagination: PaginationConfig | None) -> list[Record]:
        if raw is None:
            return []
        if pagination and pagination.page_key:
            items_ref = BaseExternalDataClient._resolve_path(raw, pagination.page_key)
            if items_ref is None:
                return []
            items = items_ref if isinstance(items_ref, list) else []
            return [item if isinstance(item, Mapping) else {"value": item} for item in items]
        if isinstance(raw, Mapping):
            return [raw]  # type: ignore[list-item]
        if isinstance(raw, list):
            return [item if isinstance(item, Mapping) else {"value": item} for item in raw]
        return [{"value": raw}]

    @staticmethod
    def _extract_next_cursor(raw: Any, pagination: PaginationConfig | None) -> str | None:
        if not pagination:
            return None
        if pagination.next_key:
            next_ref = BaseExternalDataClient._resolve_path(raw, pagination.next_key)
            if next_ref:
                return str(next_ref)
        return None

    @staticmethod
    def _resolve_path(payload: Any, dotted_path: str) -> Any:
        if not isinstance(payload, Mapping):
            return None
        if "." not in dotted_path:
            return payload.get(dotted_path)
        current: Any = payload
        for part in dotted_path.split("."):
            if not isinstance(current, Mapping):
                return None
            current = current.get(part)
            if current is None:
                return None
        return current


class ConfiguredExternalDataClient(BaseExternalDataClient):
    """Базовый клиент, инициализируемый pydantic-конфигурацией."""

    def __init__(self, *, config: ClientConfig, transport: Any) -> None:
        super().__init__(
            base_url=config.normalized_base_url,
            routes=config.routes,
            default_params=config.params,
            pagination=config.pagination,
            transport=transport,
        )
        self.config = config

    def metadata(self) -> Mapping[str, Any]:
        return {
            "name": self.config.name,
            "base_url": self.config.normalized_base_url,
            "routes": dict(self.config.routes),
        }
