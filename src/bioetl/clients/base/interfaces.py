from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Iterator, Mapping, Protocol

from bioetl.clients.base.client import (
    ClientRequest,
    ExternalDataClient,
    Page,
    Record,
    RequestContext,
)


class TransportError(RuntimeError):
    """Error raised when an HTTP transport operation fails."""


class DataProviderError(RuntimeError):
    """Error raised when data provider operations fail."""


@dataclass(slots=True)
class PaginationParams:
    """Pagination configuration passed to provider implementations."""

    page_size: int | None = None
    page_key: str | None = "items"
    next_key: str | None = "next"
    page_param: str | None = "page"


class LoggingTransportAdapter:
    """Simple adapter to bind a request context for logging or tracing purposes."""

    def __init__(self) -> None:
        self.current_context: RequestContext | None = None

    def set_context(self, context: RequestContext | None) -> None:
        self.current_context = context


class DataProviderProtocol(Protocol):
    """Protocol describing the expected data provider surface."""

    def fetch_one(
        self,
        identifier: str | Mapping[str, Any],
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
    ) -> Iterator[Page]:
        ...

    def fetch_many(
        self,
        *,
        params: Mapping[str, Any] | None = None,
        page_size: int | None = None,
        context: RequestContext | None = None,
    ) -> Iterator[Mapping[str, Any]]:
        ...

    def metadata(self) -> Mapping[str, Any]:
        ...

    def close(self) -> None:
        ...


class ClientProtocol(ExternalDataClient, Protocol):
    """Backward-compatible alias for the unified ExternalDataClient contract."""

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


class BaseApiClient(Protocol):
    """Protocol describing the expected Base API client surface."""

    def request(
        self,
        method: str,
        url: str,
        *,
        params: Mapping[str, Any] | None = None,
        json: Any | None = None,
        headers: Mapping[str, str] | None = None,
        context: RequestContext | None = None,
    ) -> Any:
        ...

    def get_json(
        self,
        url: str,
        *,
        params: Mapping[str, Any] | None = None,
        context: RequestContext | None = None,
    ) -> Any:
        ...

    def paginate_json(
        self,
        url: str,
        *,
        params: Mapping[str, Any] | None = None,
        page_key: str | None = None,
        next_key: str | None = None,
        page_param: str | None = None,
        context: RequestContext | None = None,
    ) -> Iterator[Any]:
        ...

    def close(self) -> None:
        ...
