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
from typing import TYPE_CHECKING, Protocol, cast

if TYPE_CHECKING:
    from pandera.api.pandas.model import (  # type: ignore[attr-defined]
        Check as PanderaCheck,
        Column as PanderaColumn,
        DataFrameModel as PanderaDataFrameModel,
        DataFrameSchema as PanderaDataFrameSchema,
    )


class _PanderaModule(Protocol):
    Column: type["PanderaColumn"]
    Check: type["PanderaCheck"]
    DataFrameSchema: type["PanderaDataFrameSchema"]
    DataFrameModel: type["PanderaDataFrameModel"]

__all__ = [
    "load_pandera_pandas",
    "pa",
    "Column",
    "Check",
    "DataFrameSchema",
    "DataFrameModel",
]


def load_pandera_pandas() -> _PanderaModule:
    """Return the ``pandera.pandas`` module, installing a shim when required."""

    try:
        return cast(_PanderaModule, import_module("pandera.pandas"))
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
        return cast(_PanderaModule, shim)


pa: _PanderaModule = load_pandera_pandas()

Column = pa.Column
Check = pa.Check
DataFrameSchema = pa.DataFrameSchema
DataFrameModel = pa.DataFrameModel
