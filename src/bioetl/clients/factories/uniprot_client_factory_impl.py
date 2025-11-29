from __future__ import annotations

from typing import Any, Mapping

from bioetl.clients.base import BaseHttpClientABC, PaginationParams
from bioetl.clients.base.normalizers import INormalizer
from bioetl.clients.factories import ClientFactory, build_normalizer, build_transport
from bioetl.clients.providers.uniprot.uniprot_data_client_impl import UniprotDataClientImpl
from bioetl.clients.providers.uniprot.uniprot_normalizer_impl import UniprotNormalizerImpl


class UniprotClientFactoryImpl(ClientFactory):
    """Фабрика клиентов UniProt."""

    def create(self, config: Mapping[str, Any]) -> UniprotDataClientImpl:
        transport: BaseHttpClientABC = build_transport(config)

        normalizer: INormalizer
        try:
            normalizer = build_normalizer(config)
        except ValueError:
            normalizer = UniprotNormalizerImpl()

        pagination_config = config.get("pagination") if isinstance(config, Mapping) else None
        pagination = None
        if isinstance(pagination_config, Mapping):
            pagination = PaginationParams(
                page_key=pagination_config.get("page_key", "results"),
                next_key=pagination_config.get("next_key", "next"),
                page_param=pagination_config.get("page_param", "cursor"),
                page_size=pagination_config.get("page_size"),
            )

        client = UniprotDataClientImpl(
            transport,
            pagination_strategy=config.get("pagination_strategy"),
            normalizer=normalizer,
            options=config.get("options"),
        )

        if pagination is not None:
            client.configure(pagination=pagination)

        return client
