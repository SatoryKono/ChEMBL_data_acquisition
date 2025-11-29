from __future__ import annotations

from unittest.mock import Mock

from src.bioetl.clients.base.http import BaseHttpClientABC
from src.bioetl.clients.base.pagination import TransportPaginationStrategyImpl


def test_transport_pagination_strategy_iterates_payloads() -> None:
    transport = Mock(spec=BaseHttpClientABC)
    transport.paginate_json.return_value = iter([[1], [2], [3]])
    strategy = TransportPaginationStrategyImpl()

    result = list(
        strategy.iter_pages(
            None,
            transport,
            endpoint="/items",
            params={"q": 1},
            page_key="items",
            next_key="next",
            page_param="page",
        )
    )

    assert result == [[1], [2], [3]]
    transport.paginate_json.assert_called_once_with(
        "/items",
        params={"q": 1},
        page_key="items",
        next_key="next",
        page_param="page",
    )
