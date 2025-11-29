from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Mapping, Protocol, runtime_checkable


class TransportError(RuntimeError):
    """Ошибка транспорта при обращении к внешнему API."""


class DataProviderError(RuntimeError):
    """Ошибка, возникающая при обработке данных провайдера."""


@dataclass(slots=True)
class PaginationConfig:
    """Настройки пагинации для клиентов внешних источников."""

    page_key: str | None = None
    next_key: str | None = None
    page_param: str | None = None
    page_size_param: str | None = None


@dataclass(slots=True)
class LoggingTransportAdapter:
    """Адаптер для проброса контекста в транспорт и логирование."""

    current_context: Mapping[str, Any] | None = None

    def set_context(self, context: Mapping[str, Any] | None) -> None:
        self.current_context = context


@runtime_checkable
class SupportsClose(Protocol):
    """Протокол, описывающий наличие метода close."""

    def close(self) -> None:
        ...
