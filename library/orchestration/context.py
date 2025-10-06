"""Runtime helpers for orchestrating shared ETL resources."""

from __future__ import annotations

from collections.abc import Callable
from contextlib import ExitStack
from types import TracebackType

from ..clients import ChemblClient, PubMedClient
from ..common.rate_limiter import RateLimiter, get_global_limiter
from ..config import Config

ChemblClientFactory = Callable[[Config, RateLimiter | None], ChemblClient]
PubMedClientFactory = Callable[[Config], PubMedClient]

__all__ = ["ETLContext", "ChemblClientFactory", "PubMedClientFactory"]


class ETLContext:
    """Lazy container that provisions service clients for ETL pipelines."""

    def __init__(
        self,
        cfg: Config,
        *,
        chembl_client_factory: ChemblClientFactory | None = None,
        pubmed_client_factory: PubMedClientFactory | None = None,
    ) -> None:
        self.cfg = cfg
        self._stack = ExitStack()
        self._chembl_client_factory = (
            chembl_client_factory if chembl_client_factory is not None else self._build_chembl_client
        )
        self._pubmed_client_factory = (
            pubmed_client_factory if pubmed_client_factory is not None else self._build_pubmed_client
        )
        self._chembl_client: ChemblClient | None = None
        self._chembl_handle: ChemblClient | None = None
        self._pubmed_client: PubMedClient | None = None
        self._global_limiter = self._build_global_limiter(cfg)

    @staticmethod
    def _build_global_limiter(cfg: Config) -> RateLimiter | None:
        rate_cfg = cfg.rate
        rps = getattr(rate_cfg, "global_rps", None)
        if rps is None or rps <= 0:
            return None
        burst = getattr(rate_cfg, "global_burst", None)
        return get_global_limiter(rps, burst)

    def _build_chembl_client(
        self, cfg: Config, limiter: RateLimiter | None
    ) -> ChemblClient:
        return ChemblClient(cfg.api, cfg.retry, cfg.chembl, global_limiter=limiter)

    def _build_pubmed_client(self, cfg: Config) -> PubMedClient:
        return PubMedClient(cfg.pubmed)

    @property
    def global_limiter(self) -> RateLimiter | None:
        return self._global_limiter

    @property
    def chembl(self) -> ChemblClient:
        if self._chembl_handle is None:
            client = self._chembl_client_factory(self.cfg, self._global_limiter)
            self._chembl_client = client
            self._chembl_handle = self._stack.enter_context(client)
        return self._chembl_handle

    @property
    def pubmed(self) -> PubMedClient:
        if self._pubmed_client is None:
            self._pubmed_client = self._pubmed_client_factory(self.cfg)
        return self._pubmed_client

    def close(self) -> None:
        pubmed_client = self._pubmed_client
        self._pubmed_client = None
        if pubmed_client is not None:
            close = getattr(pubmed_client, "close", None)
            if callable(close):
                close()
        self._stack.close()
        self._chembl_handle = None
        self._chembl_client = None

    def __enter__(self) -> "ETLContext":
        return self

    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        exc: BaseException | None,
        tb: TracebackType | None,
    ) -> bool:
        self.close()
        return False
