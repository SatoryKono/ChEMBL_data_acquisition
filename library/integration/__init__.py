"""Adapters and façade modules wrapping third-party integrations."""

from __future__ import annotations

from importlib import import_module
from typing import Any

from .chembl_client import ChemblClient

_LAZY_MODULES = {
  "chembl_library": "library.integration.chembl_library",
  "input_initialisation_library": "library.integration.input_initialisation_library",
  "iuphar_library": "library.integration.iuphar_library",
  "mapper_batch_library": "library.integration.mapper_batch_library",
  "mapper_library": "library.integration.mapper_library",
  "molecule_catalog": "library.integration.molecule_catalog",
  "openalex_crossref_library": "library.integration.openalex_crossref_library",
  "pubchem_library": "library.integration.pubchem_library",
  "pubmed_library": "library.integration.pubmed_library",
  "semantic_scholar_library": "library.integration.semantic_scholar_library",
  "uniprot_library": "library.integration.uniprot_library",
}

__all__ = [
  "ChemblClient",
  *_LAZY_MODULES,
]


def __getattr__(name: str) -> Any:
  """Dynamically import integration submodules on first access."""

  if name in _LAZY_MODULES:
    module = import_module(_LAZY_MODULES[name])
    globals()[name] = module
    return module
  raise AttributeError(f"module '{__name__}' has no attribute '{name}'")


def __dir__() -> list[str]:
  """Return available attributes including lazily loaded submodules."""

  return sorted({*globals(), *__all__})
