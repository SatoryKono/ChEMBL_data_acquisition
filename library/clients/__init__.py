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

from .base import (
    BasePaginatedClient,
    ClientConfig,
    ClientError,
    ClientPayload,
    ClientProtocol,
    RateLimitedRequestMixin,
)
from .chembl import ChemblClient, _chunked
from .crossref import fetch_crossref
from .iuphar import (
    download_gtp_to_hgnc_mapping,
    download_gtp_to_uniprot_mapping,
)
from .iuphar import (
    init_session as init_iuphar_session,
)
from .iuphar import (
    load_families as load_iuphar_families,
)
from .iuphar import (
    load_targets as load_iuphar_targets,
)
from .iuphar import (
    query_gene_symbol as query_iuphar_gene_symbol,
)
from .openalex import fetch_openalex
from .pubmed import PubMedClient
from .semantic_scholar import (
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
    "RateLimitedRequestMixin",
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


