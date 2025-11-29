from __future__ import annotations

from typing import Any, Mapping

from bioetl.clients import OptionsAwareApiClientImpl
from bioetl.clients.base.interfaces import LoggingTransportAdapter
from bioetl.clients.factories.client_factory import ClientFactory
from bioetl.clients.factories.helpers import build_normalizer
from bioetl.clients.providers.pubmed import PubMedDataClientImpl, PubMedNormalizerImpl


class PubMedClientFactoryImpl(ClientFactory):
    """Фабрика построения PubMed data client."""

    def create(self, config: Mapping[str, Any]) -> PubMedDataClientImpl:
        base_client = config.get("base_client")
        adapter = config.get("adapter")
        normalizer_value = config.get("normalizer") or PubMedNormalizerImpl()
        normalizer = build_normalizer({"normalizer": normalizer_value})

        if not isinstance(adapter, LoggingTransportAdapter):
            raise ValueError("Adapter must be a LoggingTransportAdapter instance")

        required_attrs = ("request", "get_json", "paginate_json", "close")
        if not all(hasattr(base_client, attr) for attr in required_attrs):
            raise ValueError("base_client must implement BaseApiClient protocol")

        transport = OptionsAwareApiClientImpl(base_client, adapter)
        return PubMedDataClientImpl(transport, normalizer=normalizer, logger=config.get("logger"))
