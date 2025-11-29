from __future__ import annotations

from typing import Any, Mapping, Protocol

from bioetl.clients.base.interfaces import DataProviderProtocol


class ClientFactory(Protocol):
    def create(self, config: Mapping[str, Any]) -> DataProviderProtocol:
        ...
