"""Utility functions for measuring and logging execution durations."""

from __future__ import annotations

from time import perf_counter

from .log import logger


def log_duration(start: float) -> float:
    """Log the elapsed time since ``start``.

    Parameters
    ----------
    start : float
        Start time as returned by :func:`time.perf_counter`.

    Returns
    -------
    float
        The elapsed duration in seconds.
    """
    duration = perf_counter() - start
    logger.info("duration_sec", extra={"value": duration})
    return duration


__all__ = ["log_duration"]
