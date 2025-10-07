"""Runtime checks for interpreter compatibility.

This module provides a helper ensuring the running Python interpreter meets
project requirements. The package targets Python 3.11 and newer (with CI
coverage on 3.11 and 3.12), so running the tools on an older version should
fail fast with a clear error message.
"""

from __future__ import annotations

import sys

from .common.log import logger


def require_python_version(min_version: tuple[int, int] = (3, 11)) -> None:
    """Validate that the interpreter meets ``min_version``.

    Parameters
    ----------
    min_version:
        Required ``(major, minor)`` Python version. Defaults to ``(3, 11)``.

    Raises
    ------
    RuntimeError
        If the current interpreter is older than ``min_version``.
    """
    current_info = sys.version_info
    current = (current_info.major, current_info.minor, current_info.micro)
    required = min_version + (0,)
    if current < required:
        current_str = f"{current_info.major}.{current_info.minor}"
        needed = f"{min_version[0]}.{min_version[1]}"
        logger.error(
            "python_version_unsupported",
            required=needed,
            found=current_str,
        )
        raise RuntimeError(f"Python {needed} or later is required; found {current_str}")
