"""Deprecated PubChem helpers for :mod:`library.pipelines.testitem`."""

from __future__ import annotations

import warnings

from library.pipelines.testitem import pubchem as _pubchem
from library.pipelines.testitem.core import _PUBCHEM_EXPORTS

warnings.warn(
    "library.testitem_pipeline.pubchem is deprecated; import from "
    "library.pipelines.testitem.pubchem instead.",
    DeprecationWarning,
    stacklevel=2,
)

for _name in _PUBCHEM_EXPORTS:
    globals()[_name] = getattr(_pubchem, _name)
del _name

__all__ = list(_PUBCHEM_EXPORTS)
