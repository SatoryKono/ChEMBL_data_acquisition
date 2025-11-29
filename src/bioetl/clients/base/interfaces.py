from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Iterator, Mapping, Protocol


class TransportError(RuntimeError):
    """Error raised when an HTTP transport operation fails."""


@dataclass(slots=True)
class RequestContext:
    """Contextual data propagated through transport calls."""

    options: Mapping[str, Any] | None = None


class LoggingTransportAdapter:
    """Simple adapter to bind a request context for logging or tracing purposes."""

    def __init__(self) -> None:
        self.current_context: RequestContext | None = None

    def set_context(self, context: RequestContext | None) -> None:
        self.current_context = context


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
