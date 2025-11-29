"""Базовые абстракции клиентской подсистемы."""

from bioetl.clients.base.client_abc import (
    BaseExternalDataClient,
    ClientRequest,
    ConfiguredExternalDataClient,
    ExternalDataClient,
    Page,
    Record,
    RequestContext,
)
from bioetl.clients.base.http_backend import BaseHttpBackend, RequestsHttpBackend
from bioetl.clients.base.types import LoggingTransportAdapter, PaginationConfig, DataProviderError, TransportError

__all__ = [
    "BaseExternalDataClient",
    "ClientRequest",
    "ConfiguredExternalDataClient",
    "DataProviderError",
    "ExternalDataClient",
    "LoggingTransportAdapter",
    "Page",
    "PaginationConfig",
    "Record",
    "RequestContext",
    "RequestsHttpBackend",
    "BaseHttpBackend",
    "TransportError",
]
