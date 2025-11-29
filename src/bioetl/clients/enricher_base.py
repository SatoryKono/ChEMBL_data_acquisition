from __future__ import annotations

from typing import Any, Iterator, Mapping

from .base.http import BaseHttpClientABC
from .base.interfaces import BaseApiClient, LoggingTransportAdapter, RequestContext, TransportError


class OptionsAwareApiClientImpl(BaseHttpClientABC):
    """Adapter exposing ``BaseApiClient`` via the transport abstraction."""

    def __init__(
        self,
        base_client: BaseApiClient,
        adapter: LoggingTransportAdapter,
    ) -> None:
        self._base_client = base_client
        self._adapter = adapter

    def _apply_context(self, context: RequestContext | None) -> None:
        self._adapter.set_context(context)

    def _proxy_call(self, func: Any, *args: Any, **kwargs: Any) -> Any:
        try:
            return func(*args, **kwargs)
        except TransportError:
            raise
        except Exception as exc:  # noqa: BLE001
            raise TransportError(str(exc)) from exc

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
        self._apply_context(context)
        return self._proxy_call(
            self._base_client.request,
            method,
            url,
            params=params,
            json=json,
            headers=headers,
            context=context,
        )

    def get_json(
        self,
        url: str,
        *,
        params: Mapping[str, Any] | None = None,
        context: RequestContext | None = None,
    ) -> Any:
        self._apply_context(context)
        return self._proxy_call(
            self._base_client.get_json,
            url,
            params=params,
            context=context,
        )

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
        self._apply_context(context)
        return self._proxy_call(
            self._base_client.paginate_json,
            url,
            params=params,
            page_key=page_key,
            next_key=next_key,
            page_param=page_param,
            context=context,
        )

    def close(self) -> None:
        self._proxy_call(self._base_client.close)
