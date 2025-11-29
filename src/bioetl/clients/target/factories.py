from __future__ import annotations

from typing import Any, Mapping

from bioetl.clients.base import BaseHttpClientABC, DataProviderProtocol, PaginationParams, TransportPaginationStrategyImpl
from bioetl.clients.providers.iuphar.impl.target_client_impl import TargetIupharDataProviderImpl
from bioetl.clients.providers.iuphar.normalization import ProviderNormalizer


def default_target_iuphar(
    transport: BaseHttpClientABC,
    *,
    options: Mapping[str, Any] | None = None,
    logger: Any | None = None,
) -> DataProviderProtocol:
    """Фабрика по умолчанию для IUPHAR target (контракт DataProviderProtocol).

    Имплементация провайдера расположена в ``providers/iuphar/impl`` и
    регистрируется через эту фабрику.
    """

    pagination = PaginationParams(page_key=TargetIupharDataProviderImpl.DEFAULT_PAGINATION.page_key)

    return TargetIupharDataProviderImpl(
        transport,
        pagination_strategy=TransportPaginationStrategyImpl(),
        normalizer=ProviderNormalizer(),
        logger=logger,
        options={"page_key": pagination.page_key, **(options or {})},
    )
