from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Iterable, Iterator, Mapping

from bioetl.clients.base.client_abc import Page
from bioetl.clients.base.types import PaginationConfig


@dataclass(slots=True)
class CursorState:
    """Состояние пагинации во время обхода."""

    cursor: str | None = None


def iter_with_pagination(
    transport: Any,
    *,
    url: str,
    params: Mapping[str, Any] | None,
    pagination: PaginationConfig | None,
    context: Any | None = None,
) -> Iterator[Page]:
    """Обёртка для итерации по страницам при помощи транспорта."""

    state = CursorState()
    base_params = dict(params or {})
    while True:
        effective_params = dict(base_params)
        if pagination and pagination.page_param and state.cursor:
            effective_params[pagination.page_param] = state.cursor
        payload = transport.get_json(url, params=effective_params, context=context)
        items = _extract_items(payload, pagination)
        state.cursor = _extract_next(payload, pagination)
        yield Page(items=items, next_cursor=state.cursor, raw=payload)
        if not state.cursor:
            break


def _extract_items(payload: Any, pagination: PaginationConfig | None) -> list[Mapping[str, Any]]:
    if payload is None:
        return []
    if pagination and pagination.page_key and isinstance(payload, Mapping):
        raw_items = payload.get(pagination.page_key, [])
        return _normalize_items(raw_items)
    if isinstance(payload, Iterable) and not isinstance(payload, (str, bytes)):
        return _normalize_items(payload)
    if isinstance(payload, Mapping):
        return [payload]
    return [{"value": payload}]


def _normalize_items(items: Iterable[Any]) -> list[Mapping[str, Any]]:
    normalized: list[Mapping[str, Any]] = []
    for item in items:
        if isinstance(item, Mapping):
            normalized.append(item)
        else:
            normalized.append({"value": item})
    return normalized


def _extract_next(payload: Any, pagination: PaginationConfig | None) -> str | None:
    if not pagination or not pagination.next_key:
        return None
    if isinstance(payload, Mapping):
        next_value = payload.get(pagination.next_key)
        if next_value:
            return str(next_value)
    return None
