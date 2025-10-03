"""Legacy compatibility wrapper for test item PubChem helpers.

Changelog
========
* Re-export implementation from :mod:`library.pipelines.testitem.pubchem` and
  warn about the deprecated import path.
"""

from __future__ import annotations

from warnings import warn

from library.pipelines.testitem.pubchem import *  # noqa: F401,F403
from library.pipelines.testitem.pubchem import __all__ as _PUBCHEM_ALL

warn(
    "library.testitem_pipeline.pubchem is deprecated; use "
    "library.pipelines.testitem.pubchem instead.",
    DeprecationWarning,
    stacklevel=2,
)

__all__ = list(_PUBCHEM_ALL)
