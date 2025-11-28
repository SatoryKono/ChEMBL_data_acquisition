from __future__ import annotations

from abc import ABC, abstractmethod
from collections.abc import Iterator, Mapping
from dataclasses import dataclass
from typing import Any


@dataclass(slots=True)
class BasePaginatedClient(ABC):
    """Base class for clients supporting limit/offset pagination."""

    limit_param: str = "limit"
    offset_param: str = "offset"
    page_size: int = 100

    def __post_init__(self) -> None:
        if self.page_size <= 0:
            raise ValueError("page_size must be positive")

    def _build_pagination_params(
        self,
        *,
        params: Mapping[str, Any] | None,
        limit: int,
        offset: int,
    ) -> dict[str, Any]:
        if limit <= 0:
            raise ValueError("limit must be positive")
        if offset < 0:
            raise ValueError("offset must be non-negative")

        effective_params: dict[str, Any] = {
            self.limit_param: limit,
            self.offset_param: offset,
        }
        if params:
            effective_params.update(params)
        return effective_params

    def iter_pages(
        self,
        endpoint: str,
        *,
        params: Mapping[str, Any] | None = None,
        start_offset: int = 0,
    ) -> Iterator[Mapping[str, Any]]:
        """Yield subsequent pages until the API stops providing a next link."""

        offset = start_offset
        if offset < 0:
            raise ValueError("start_offset must be non-negative")

        while True:
            page_params = self._build_pagination_params(
                params=params, limit=self.page_size, offset=offset
            )
            payload = self.request_json(endpoint, params=page_params)
            yield payload

            page_meta = payload.get("page_meta") if isinstance(payload, Mapping) else None
            next_url = page_meta.get("next") if isinstance(page_meta, Mapping) else None
            if not next_url:
                break
            offset += self.page_size

    @abstractmethod
    def request_json(
        self,
        endpoint: str,
        *,
        params: Mapping[str, Any] | None = None,
        **kwargs: Any,
    ) -> Mapping[str, Any]:
        """Perform a JSON request and return a response payload."""

