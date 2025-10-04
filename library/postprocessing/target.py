"""Lightweight hooks for final target export post-processing."""

from __future__ import annotations

from pathlib import Path

from ..common.log import logger


def process_targets(path: str | Path, *, verbose: bool = False) -> Path:
    """Emit a post-processing hook event for the final target export.

    Parameters
    ----------
    path:
        Location of the produced target table.
    verbose:
        When ``True`` the hook logs at ``INFO`` level, otherwise ``DEBUG``.

    Returns
    -------
    pathlib.Path
        The resolved path to the processed file.
    """

    resolved = Path(path)
    log_fn = logger.info if verbose else logger.debug
    log_fn("target_postprocess", path=str(resolved))
    return resolved


__all__ = ["process_targets"]

