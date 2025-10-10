"""Compatibility wrapper around :mod:`library.postprocess.common`."""

from __future__ import annotations

from library.postprocess import common as _common
from library.postprocess.common import *  # noqa: F401,F403

__all__ = list(getattr(_common, "__all__", ()))
