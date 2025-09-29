"""Client-facing abstractions for third-party integrations.

This package defines the public interface for service clients used across the
ChEMBL data acquisition stack. Downstream modules should import client
contracts from ``library.clients`` instead of reaching into concrete client
implementations. This keeps the dependency direction flowing outwards from the
pipelines toward infrastructure and simplifies substitution during testing.

Import rules for higher layers:
* Business logic in ``library`` may depend on these exported protocols and
  dataclasses.
* Concrete pipeline steps must not import private modules from future client
  packages directly; rely on the abstractions defined here instead.
* Test suites may import concrete clients for fixtures, but prefer the
  interfaces below for type-safe mocking.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Mapping, MutableMapping, Protocol, Sequence

from .crossref import fetch_crossref
from .openalex import fetch_openalex

__all__ = [
    "ClientConfig",
    "ClientError",
    "ClientPayload",
    "ClientProtocol",
    "fetch_crossref",
    "fetch_openalex",
]


@dataclass(slots=True)
class ClientConfig:
    """Common configuration options for HTTP-like service clients."""

    base_url: str
    headers: Mapping[str, str]
    timeout_seconds: float


@dataclass(slots=True)
class ClientPayload:
    """Normalized response payload shared by client implementations."""

    data: Sequence[MutableMapping[str, Any]]


class ClientError(RuntimeError):
    """Raised when a client encounters an unrecoverable error."""


class ClientProtocol(Protocol):
    """Skeleton protocol for strongly typed service clients."""

    config: ClientConfig

    def fetch(self, *, params: Mapping[str, Any] | None = None) -> ClientPayload:
        """Retrieve data from the upstream service."""

