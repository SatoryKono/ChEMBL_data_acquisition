"""Compatibility helpers for importing :mod:`pandera` namespaces.

Pandera reorganised its public modules in the 0.20 series moving the
``pandera.pandas`` namespace to ``pandera.api.pandas``.  The rest of the code
base still imports :mod:`pandera.pandas` directly which raises a
``ModuleNotFoundError`` on new releases when the legacy shim is absent.  The
tests already emulate the old location; mirror the behaviour at runtime so the
CLI tools keep working without requiring an explicit dependency shim during
installation.
"""

from __future__ import annotations

import sys
import types
from importlib import import_module
from types import ModuleType

__all__ = ["load_pandera_pandas", "pa"]


def load_pandera_pandas() -> ModuleType:
    """Return the ``pandera.pandas`` module, installing a shim when required."""

    try:
        return import_module("pandera.pandas")
    except ModuleNotFoundError as exc:
        if exc.name != "pandera.pandas":
            raise

        pandas_model = import_module("pandera.api.pandas.model")
        shim = types.ModuleType("pandera.pandas")
        for attr in dir(pandas_model):
            if attr.startswith("__"):
                continue
            setattr(shim, attr, getattr(pandas_model, attr))

        sys.modules.setdefault("pandera.pandas", shim)
        return shim


pa = load_pandera_pandas()
