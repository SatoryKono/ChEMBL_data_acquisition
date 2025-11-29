from __future__ import annotations

from typing import Any, Mapping

from bioetl.clients.enricher_base import OptionsAwareApiClientImpl
from bioetl.clients.base.interfaces import LoggingTransportAdapter
from bioetl.clients.base.normalizers import INormalizer
from bioetl.clients.base.pagination import PaginationStrategyABC, TransportPaginationStrategyImpl
from bioetl.clients.factories.client_factory import ClientFactory
from bioetl.clients.factories.helpers import build_normalizer
from bioetl.clients.providers.pubchem.pubchem_data_client_impl import PubChemDataClientImpl
from bioetl.clients.providers.pubchem.pubchem_normalizer_impl import PubChemNormalizerImpl


class PubChemClientFactoryImpl(ClientFactory):
    """Фабрика клиентов PubChem."""

    def create(self, config: Mapping[str, Any]) -> PubChemDataClientImpl:
        api_client = config.get("api_client")
        if api_client is None:
            raise ValueError("api_client must be provided to build PubChem client")

        adapter: LoggingTransportAdapter = config.get("adapter") or LoggingTransportAdapter()
        http_client = OptionsAwareApiClientImpl(base_client=api_client, adapter=adapter)

        normalizer: INormalizer = config.get("normalizer") or build_normalizer(
            {"normalizer": PubChemNormalizerImpl()}
        )
        pagination_strategy: PaginationStrategyABC | None = config.get("pagination_strategy")

        return PubChemDataClientImpl(
            http_client,
            normalizer=normalizer,
            pagination_strategy=pagination_strategy or TransportPaginationStrategyImpl(),
            logger=config.get("logger"),
            options=config.get("options"),
        )
