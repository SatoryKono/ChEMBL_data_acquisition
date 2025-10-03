"""Adapters and façade modules wrapping third-party integrations."""

from __future__ import annotations

from importlib import import_module
from typing import Any

_LAZY_MODULES = {
  "chembl_library",
  "input_initialisation_library",
  "iuphar_library",
  "mapper_batch_library",
  "mapper_library",
  "molecule_catalog",
  "openalex_crossref_library",
  "pubchem_library",
  "pubmed_library",
  "semantic_scholar_library",
  "uniprot_library",
}

__all__ = [
  "ChemblClient",
  *_LAZY_MODULES,
]


def __getattr__(name: str) -> Any:
  """Lazily import integration helpers on first access."""

  if name == "ChemblClient":
    from .chembl_client import ChemblClient as _ChemblClient

    globals()[name] = _ChemblClient
    return _ChemblClient
  if name in _LAZY_MODULES:
    module = import_module(f"{__name__}.{name}")
    globals()[name] = module
    return module
  raise AttributeError(f"module '{__name__}' has no attribute '{name}'")


def __dir__() -> list[str]:
  """Return available attributes including lazily loaded modules."""

  return sorted({*globals(), *__all__, *_LAZY_MODULES})
