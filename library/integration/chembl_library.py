"""Compatibility layer aggregating ChEMBL helpers."""

from __future__ import annotations

from importlib import import_module
from typing import Any

from library.clients import _chunked

_PIPELINE_EXPORTS = {
  "get_assay": ("library.pipelines.assay", "get_assay"),
  "get_assays": ("library.pipelines.assay", "get_assays"),
  "get_activities": ("library.pipelines.assay", "get_activities"),
  "get_testitem": ("library.pipelines.assay", "get_testitem"),
  "get_documents": ("library.pipelines.document", "get_documents"),
  "get_target": ("library.pipelines.target.chembl_target", "get_target"),
  "get_target_payload": ("library.pipelines.target.chembl_target", "get_target_payload"),
  "get_targets": ("library.pipelines.target.chembl_target", "get_targets"),
  "get_targets_payloads": ("library.pipelines.target.chembl_target", "get_targets_payloads"),
  "get_targets_raw_frame": ("library.pipelines.target.chembl_target", "get_targets_raw_frame"),
  "iter_target_batches": ("library.pipelines.target.chembl_target", "iter_target_batches"),
  "extend_target": ("library.pipelines.target.chembl_target", "extend_target"),
}

__all__ = [
  *_PIPELINE_EXPORTS,
  "_chunked",
]


def __getattr__(name: str) -> Any:
  """Return lazily imported pipeline helpers."""

  if name in _PIPELINE_EXPORTS:
    module_name, attribute = _PIPELINE_EXPORTS[name]
    module = import_module(module_name)
    value = getattr(module, attribute)
    globals()[name] = value
    return value
  raise AttributeError(f"module '{__name__}' has no attribute '{name}'")


def __dir__() -> list[str]:
  """Return available attributes for introspection tools."""

  return sorted({*globals(), *__all__})
