"""Runtime context for ETL scripts and pipelines."""

from __future__ import annotations

from contextlib import ExitStack, AbstractContextManager
from types import TracebackType
from typing import Callable, ContextManager

from ..clients import ChemblClient
from ..common.rate_limiter import RateLimiter, get_global_limiter
from ..config import Config


class ETLContext(AbstractContextManager["ETLContext"]):
    """Create and manage shared service clients for ETL runs."""

    def __init__(
        self,
        cfg: Config,
        *,
        chembl_client_factory: Callable[..., ContextManager[ChemblClient]] | None = None,
        global_limiter_factory: Callable[[float | None, int | None], RateLimiter | None]
        | None = None,
    ) -> None:
        self._cfg = cfg
        self._chembl_factory = chembl_client_factory or self._default_chembl_factory
        self._limiter_factory = global_limiter_factory or get_global_limiter
        self._stack: ExitStack | None = None
        self._chembl_client: ChemblClient | None = None
        self._global_limiter: RateLimiter | None = None

    def _default_chembl_factory(
        self,
        api: object,
        retry: object,
        chembl: object,
        *,
        global_limiter: RateLimiter | None = None,
    ) -> ContextManager[ChemblClient]:
        return ChemblClient(api, retry, chembl, global_limiter=global_limiter)

    def __enter__(self) -> "ETLContext":
        if self._stack is not None:
            raise RuntimeError("ETLContext cannot be re-entered")
        self._stack = ExitStack()
        rate_cfg = self._cfg.rate
        self._global_limiter = self._limiter_factory(
            getattr(rate_cfg, "global_rps", None),
            getattr(rate_cfg, "global_burst", None),
        )
        self._chembl_client = self._stack.enter_context(
            self._chembl_factory(
                self._cfg.api,
                self._cfg.retry,
                self._cfg.chembl,
                global_limiter=self._global_limiter,
            )
        )
        return self

    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        exc: BaseException | None,
        tb: TracebackType | None,
    ) -> None:
        if self._stack is not None:
            self._stack.__exit__(exc_type, exc, tb)
        self._stack = None
        self._chembl_client = None
        self._global_limiter = None
        return None

    @property
    def chembl_client(self) -> ChemblClient:
        if self._chembl_client is None:
            raise RuntimeError("Chembl client is not initialised")
        return self._chembl_client

    @property
    def global_limiter(self) -> RateLimiter | None:
        return self._global_limiter

    def close(self) -> None:
        if self._stack is not None:
            self._stack.close()
        self._stack = None
        self._chembl_client = None
        self._global_limiter = None
