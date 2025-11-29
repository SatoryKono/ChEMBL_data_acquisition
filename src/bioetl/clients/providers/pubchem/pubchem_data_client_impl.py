from __future__ import annotations

import warnings
from typing import Any, Iterator, Mapping

from bioetl.clients.base.interfaces import PaginationParams, RequestContext
from bioetl.clients.base.normalizers import INormalizer
from bioetl.clients.base.pagination import PaginationStrategyABC, TransportPaginationStrategyImpl
from bioetl.clients.providers._provider_template import (
    ProviderDataClientImpl,
    RouteConfig,
    normalize_payload,
)


class PubChemDataClientImpl(ProviderDataClientImpl):
    """HTTP клиент PubChem с поддержкой поиска и нормализации."""

    ROUTES = [
        RouteConfig("fetch", "/compound/cid/{value}"),
        RouteConfig("search", "/compound/name/{value}", query_param="name"),
        RouteConfig("search_smiles", "/compound/smiles", query_param="smiles"),
    ]

    def __init__(
        self,
        http: Any,
        *,
        normalizer: INormalizer | None = None,
        pagination_strategy: PaginationStrategyABC | None = None,
        logger: Any | None = None,
        options: Mapping[str, Any] | None = None,
    ) -> None:
        default_pagination = PaginationParams(page_key="PC_Compounds")
        super().__init__(
            http,
            routes=self.ROUTES,
            default_pagination=default_pagination,
            pagination_strategy=pagination_strategy or TransportPaginationStrategyImpl(),
            normalizer=normalizer,
            logger=logger,
            options=options,
        )

    def fetch(
        self,
        cid: str,
        *,
        params: Mapping[str, Any] | None = None,
        context: RequestContext | None = None,
    ) -> Iterator[Mapping[str, Any]]:
        """Запрос соединения по CID."""

        path, params_with_value = self._resolve_route("fetch", value=cid, params=params)
        payload = self._http.get_json(path, params=params_with_value, context=context)
        yield from self._normalize_records(payload)

    def search(
        self,
        value: str,
        *,
        params: Mapping[str, Any] | None = None,
        page_size: int | None = None,
        pagination: PaginationParams | None = None,
        context: RequestContext | None = None,
        query_param: str = "name",
    ) -> Iterator[Mapping[str, Any]]:
        """Поиск соединений по имени или SMILES."""

        route_name = "search_smiles" if query_param == "smiles" else "search"
        path, params_with_value = self._resolve_route(route_name, value=value, params=params)
        merged_params: dict[str, Any] = {"_path": path}
        if params_with_value:
            merged_params.update(params_with_value)

        yield from super().fetch_many(
            params=merged_params,
            page_size=page_size,
            pagination=pagination,
            context=context,
        )

    def fetch_by_cid(self, cid: str, **kwargs: Any) -> Iterator[Mapping[str, Any]]:
        warnings.warn(
            "fetch_by_cid is deprecated, use fetch('cid', ...) or fetch(...) instead",
            DeprecationWarning,
        )
        return self.fetch(cid, **kwargs)

    def search_by_smiles(self, smiles: str, **kwargs: Any) -> Iterator[Mapping[str, Any]]:
        warnings.warn(
            "search_by_smiles is deprecated, use search(..., query_param='smiles') instead",
            DeprecationWarning,
        )
        return self.search(smiles, query_param="smiles", **kwargs)

    def _normalize_records(self, payload: Any) -> Iterator[Mapping[str, Any]]:
        records = normalize_payload(payload, page_key=self._default_pagination.page_key)
        for record in records:
            yield self.normalizer.normalize(record)
