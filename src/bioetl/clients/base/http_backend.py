from __future__ import annotations

import abc
from typing import Any, Iterator, Mapping

import backoff
import requests

from bioetl.clients.base.types import PaginationConfig, TransportError


class BaseHttpBackend(abc.ABC):
    """Абстракция HTTP-бэкенда."""

    @abc.abstractmethod
    def request(
        self,
        method: str,
        url: str,
        *,
        params: Mapping[str, Any] | None = None,
        json: Any | None = None,
        headers: Mapping[str, str] | None = None,
        context: Any | None = None,
    ) -> Any:
        ...

    def get_json(
        self,
        url: str,
        *,
        params: Mapping[str, Any] | None = None,
        headers: Mapping[str, str] | None = None,
        context: Any | None = None,
    ) -> Any:
        return self.request("GET", url, params=params, headers=headers, context=context)

    def paginate_json(
        self,
        url: str,
        *,
        params: Mapping[str, Any] | None = None,
        pagination: PaginationConfig | None = None,
        headers: Mapping[str, str] | None = None,
        context: Any | None = None,
    ) -> Iterator[Any]:
        next_cursor: str | None = None
        effective_params = dict(params or {})
        while True:
            if pagination and pagination.page_param and next_cursor:
                effective_params[pagination.page_param] = next_cursor
            payload = self.get_json(url, params=effective_params, headers=headers, context=context)
            yield payload
            if not pagination or not pagination.next_key:
                break
            if isinstance(payload, Mapping):
                next_cursor = payload.get(pagination.next_key)  # type: ignore[assignment]
            else:
                next_cursor = None
            if not next_cursor:
                break

    def close(self) -> None:  # pragma: no cover - по умолчанию делать нечего
        return


class RequestsHttpBackend(BaseHttpBackend):
    """HTTP-бэкенд на базе requests с backoff."""

    def __init__(
        self,
        *,
        base_headers: Mapping[str, str] | None = None,
        timeout: float = 30.0,
        session: requests.Session | None = None,
        max_tries: int = 5,
    ) -> None:
        self._session = session or requests.Session()
        self._timeout = timeout
        self._base_headers = dict(base_headers or {})
        self._max_tries = max_tries

    @backoff.on_exception(backoff.expo, requests.RequestException, max_tries=5)
    def request(
        self,
        method: str,
        url: str,
        *,
        params: Mapping[str, Any] | None = None,
        json: Any | None = None,
        headers: Mapping[str, str] | None = None,
        context: Any | None = None,
    ) -> Any:
        merged_headers = {**self._base_headers, **(headers or {})}
        try:
            response = self._session.request(
                method,
                url,
                params=params,
                json=json,
                headers=merged_headers,
                timeout=self._timeout,
            )
            response.raise_for_status()
            if "application/json" in response.headers.get("Content-Type", ""):
                return response.json()
            return response.text
        except requests.RequestException as exc:  # pragma: no cover - зависит от сети
            raise TransportError(str(exc)) from exc

    def close(self) -> None:  # pragma: no cover - закрывает сессию
        self._session.close()
