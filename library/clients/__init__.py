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

from collections.abc import Mapping, MutableMapping, Sequence
from dataclasses import dataclass
from typing import Any, Protocol

from .base import BasePaginatedClient
from library.chembl.clients.chembl import ChemblClient, _chunked
from library.crossref.clients.crossref import fetch_crossref
from library.iuphar.clients.iuphar import (
    download_gtp_to_hgnc_mapping,
    download_gtp_to_uniprot_mapping,
)
from library.iuphar.clients.iuphar import (
    init_session as init_iuphar_session,
)
from library.iuphar.clients.iuphar import (
    load_families as load_iuphar_families,
)
from library.iuphar.clients.iuphar import (
    load_targets as load_iuphar_targets,
)
from library.iuphar.clients.iuphar import (
    query_gene_symbol as query_iuphar_gene_symbol,
)
from library.openalex.clients.openalex import fetch_openalex
from library.pubmed.clients.pubmed import PubMedClient
from library.semantic_scholar.clients.semantic_scholar import (
    fetch_semantic_scholar,
    fetch_semantic_scholar_batch,
)

__all__ = [
    "ChemblClient",
    "ClientConfig",
    "ClientError",
    "ClientPayload",
    "ClientProtocol",
    "BasePaginatedClient",
    "PubMedClient",
    "download_gtp_to_hgnc_mapping",
    "download_gtp_to_uniprot_mapping",
    "init_iuphar_session",
    "load_iuphar_families",
    "load_iuphar_targets",
    "query_iuphar_gene_symbol",
    "fetch_crossref",
    "fetch_openalex",
    "fetch_semantic_scholar",
    "fetch_semantic_scholar_batch",
    "_chunked",
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
