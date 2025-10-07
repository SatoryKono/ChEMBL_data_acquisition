"""Logging utilities for postprocessing pipelines."""
from __future__ import annotations

import logging
from typing import Optional

_DEFAULT_LOGGER_NAME = "chembl.postprocess"


def get_logger(name: Optional[str] = None) -> logging.Logger:
    """Return a configured logger for postprocessing steps.

    The logger defaults to :mod:`chembl.postprocess` and ensures that
    double-configuration is avoided when used inside notebooks or CLI tools.
    """

    logger_name = name or _DEFAULT_LOGGER_NAME
    logger = logging.getLogger(logger_name)

    if not logger.handlers:
        handler = logging.StreamHandler()
        formatter = logging.Formatter(
            fmt="%(asctime)s | %(levelname)s | %(name)s | %(message)s",
            datefmt="%Y-%m-%dT%H:%M:%S",
        )
        handler.setFormatter(formatter)
        logger.addHandler(handler)
        logger.setLevel(logging.INFO)

    logger.propagate = False
    return logger


__all__ = ["get_logger"]
