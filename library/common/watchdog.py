"""Utility helpers for logging stage SLA breaches."""

from __future__ import annotations

from contextlib import AbstractContextManager
from dataclasses import dataclass
from threading import Event, Timer
from time import perf_counter
from typing import Any, Mapping

from .log import logger


@dataclass(slots=True)
class StageWatchdog(AbstractContextManager["StageWatchdog"]):
    """Monitor execution time for long running stages.

    The watchdog starts a background timer that emits a warning when the
    configured SLA threshold is exceeded. The caller can inspect the elapsed
    time and whether the SLA was breached once the monitored block finishes.
    """

    stage: str
    event: str
    sla_seconds: float
    context: Mapping[str, Any] | None = None

    def __post_init__(self) -> None:
        if self.sla_seconds < 0:
            raise ValueError("sla_seconds must be non-negative")
        self._timer: Timer | None = None
        self._start = 0.0
        self._elapsed = 0.0
        self._breached = Event()

    def __enter__(self) -> "StageWatchdog":
        self._start = perf_counter()
        if self.sla_seconds:
            self._timer = Timer(self.sla_seconds, self._trigger)
            self._timer.daemon = True
            self._timer.start()
        return self

    def __exit__(self, exc_type, exc, tb) -> bool | None:
        self.cancel()
        self._elapsed = perf_counter() - self._start
        return None

    def cancel(self) -> None:
        """Cancel the watchdog timer if it is running."""

        if self._timer is not None:
            self._timer.cancel()
            self._timer = None

    def _trigger(self) -> None:
        self._breached.set()
        payload = {"stage": self.stage, "sla_seconds": self.sla_seconds}
        if self.context:
            payload.update(self.context)
        logger.warning(self.event, **payload)

    @property
    def breached(self) -> bool:
        """Return ``True`` when the SLA threshold was exceeded."""

        return self._breached.is_set()

    @property
    def elapsed(self) -> float:
        """Return the elapsed time in seconds for the monitored block."""

        return self._elapsed


__all__ = ["StageWatchdog"]
