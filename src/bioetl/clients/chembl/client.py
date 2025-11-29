from __future__ import annotations

from typing import Any

from bioetl.clients.base.client_abc import ConfiguredExternalDataClient
from bioetl.clients.config.models import ClientConfig


class ChemblClient(ConfiguredExternalDataClient):
    """Клиент внешнего источника, работающий по единому контракту."""

    def __init__(self, *, config: ClientConfig, transport: Any) -> None:
        super().__init__(config=config, transport=transport)
