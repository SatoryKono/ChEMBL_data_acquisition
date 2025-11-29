from __future__ import annotations

import abc
from typing import Any


class BaseDbBackend(abc.ABC):
    """Заготовка под будущие реализации DB/FTP бэкендов."""

    @abc.abstractmethod
    def fetch(self, query: str, **kwargs: Any) -> Any:
        ...

    def close(self) -> None:  # pragma: no cover - реализация зависит от конкретного драйвера
        return
