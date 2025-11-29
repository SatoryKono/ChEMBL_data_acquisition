from __future__ import annotations

import abc
from typing import Any, Iterator, Mapping

from .interfaces import RequestContext


class BaseHttpClientABC(abc.ABC):
    """Абстрактный HTTP-транспорт для провайдеров."""

    @abc.abstractmethod
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

    @abc.abstractmethod
    def get_json(
        self,
        url: str,
        *,
        params: Mapping[str, Any] | None = None,
        context: RequestContext | None = None,
    ) -> Any:
        ...

    @abc.abstractmethod
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
        """Optional close hook."""

