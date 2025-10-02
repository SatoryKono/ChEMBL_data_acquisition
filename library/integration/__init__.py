"""Adapters and façade modules wrapping third-party integrations."""

from __future__ import annotations

from .chembl_client import ChemblClient
from . import (
  chembl_library,
  input_initialisation_library,
  iuphar_library,
  mapper_batch_library,
  mapper_library,
  molecule_catalog,
  openalex_crossref_library,
  pubchem_library,
  pubmed_library,
  semantic_scholar_library,
  uniprot_library,
)

__all__ = [
  "ChemblClient",
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
]
