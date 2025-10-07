"""Runtime helpers for constructing shared ETL resources."""

from __future__ import annotations

import logging
from collections.abc import Callable
from dataclasses import dataclass, field

from library.clients import ChemblClient
from library.common.rate_limiter import RateLimiter, get_global_limiter
from library.config import Config

__all__ = ["ETLContext", "ChemblClientFactory"]

logger = logging.getLogger(__name__)

ChemblClientFactory = Callable[["ETLContext"], ChemblClient]


@dataclass
class ETLContext:
    """Context manager provisioning shared ETL resources."""

    cfg: Config
    chembl_client_factory: ChemblClientFactory | None = None
    _global_limiter: RateLimiter | None = field(init=False, default=None)
    _chembl_client: ChemblClient | None = field(init=False, default=None)
    _closers: list[Callable[[], None]] = field(init=False, default_factory=list)
    _closed: bool = field(init=False, default=False)

    def __post_init__(self) -> None:
        self._global_limiter = self._initialise_global_limiter()

    def _initialise_global_limiter(self) -> RateLimiter | None:
        rate_cfg = getattr(self.cfg, "rate", None)
        if rate_cfg is None:
            return None
        rps = getattr(rate_cfg, "global_rps", None) or 0
        if rps <= 0:
            return None
        burst = getattr(rate_cfg, "global_burst", None)
        return get_global_limiter(rps, burst)

    def __enter__(self) -> ETLContext:
        if self._closed:
            self._closed = False
        return self

    def __exit__(self, exc_type, exc, tb) -> None:
        self.close()
        return None

    def _register_cleanup(self, closer: Callable[[], None]) -> None:
        self._closers.append(closer)

    def register_cleanup(self, closer: Callable[[], None]) -> None:
        """Register ``closer`` to be executed when the context closes.

        The helper provides a public escape hatch for code that needs to tie
        additional resources to the lifecycle of the :class:`ETLContext`.
        Callers must ensure that ``closer`` is idempotent because it may be
        invoked multiple times when the context is reused across ``with``
        statements.
        """

        if self._closed:
            raise RuntimeError("Cannot register cleanup on closed context")
        self._register_cleanup(closer)

    def close(self) -> None:
        if self._closed:
            return
        while self._closers:
            closer = self._closers.pop()
            try:
                closer()
            except Exception:  # pragma: no cover - defensive cleanup
                logger.warning("etl_context_cleanup_failed", exc_info=True)
        self._chembl_client = None
        self._closed = True

    @property
    def global_limiter(self) -> RateLimiter | None:
        return self._global_limiter

    def _default_chembl_client_factory(self) -> ChemblClient:
        return ChemblClient(
            self.cfg.api,
            self.cfg.retry,
            self.cfg.chembl,
            global_limiter=self._global_limiter,
        )

    @property
    def chembl_client(self) -> ChemblClient:
        if self._closed:
            raise RuntimeError("ETLContext is closed")
        if self._chembl_client is None:
            if self.chembl_client_factory is None:
                client = self._default_chembl_client_factory()
            else:
                client = self.chembl_client_factory(self)
            close = getattr(client, "close", None)
            if callable(close):
                self._register_cleanup(close)
            self._chembl_client = client
        return self._chembl_client
