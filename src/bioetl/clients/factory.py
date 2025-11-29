from __future__ import annotations

from pathlib import Path
from typing import Mapping, Type

from bioetl.clients.base.client_abc import BaseExternalDataClient
from bioetl.clients.base.http_backend import RequestsHttpBackend
from bioetl.clients.config.loader import load_client_config

from bioetl.clients.chembl.client import ChemblClient
from bioetl.clients.pubchem.client import PubChemClient
from bioetl.clients.pubmed.client import PubMedClient
from bioetl.clients.openalex.client import OpenAlexClient
from bioetl.clients.crossref.client import CrossrefClient
from bioetl.clients.uniprot.client import UniprotClient
from bioetl.clients.semantic_scholar.client import SemanticScholarClient
from bioetl.clients.iuphar.client import IupharClient

ClientClass = Type[BaseExternalDataClient]


CLIENTS_REGISTRY: Mapping[str, ClientClass] = {
    "chembl": ChemblClient,
    "pubchem": PubChemClient,
    "pubmed": PubMedClient,
    "openalex": OpenAlexClient,
    "crossref": CrossrefClient,
    "uniprot": UniprotClient,
    "semantic_scholar": SemanticScholarClient,
    "iuphar": IupharClient,
}


def create_client(source: str, *, config_path: str | Path | None = None) -> BaseExternalDataClient:
    """Создаёт клиента по его имени и YAML-конфигурации."""

    source_lower = source.lower()
    if source_lower not in CLIENTS_REGISTRY:
        raise KeyError(f"Неизвестный источник: {source}")

    config = load_client_config(source_lower, path=Path(config_path) if config_path else None)
    http_backend = RequestsHttpBackend(base_headers=config.headers, timeout=config.timeout)

    client_cls = CLIENTS_REGISTRY[source_lower]
    return client_cls(config=config, transport=http_backend)
