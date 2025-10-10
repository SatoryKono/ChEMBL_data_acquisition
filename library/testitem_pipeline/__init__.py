"""Legacy shim for the test item pipeline package.

This module preserves backwards compatibility for consumers that still import
``library.testitem_pipeline`` while the implementation now lives in
:mod:`library.pipelines.testitem`.
"""

from __future__ import annotations

import warnings

from library.pipelines.testitem import *  # noqa: F401,F403
from library.pipelines.testitem import __all__ as _PUBLIC_API

warnings.warn(
    "'library.testitem_pipeline' is deprecated; please import from "
    "'library.pipelines.testitem' instead.",
    DeprecationWarning,
    stacklevel=2,
)

__all__ = list(_PUBLIC_API)
