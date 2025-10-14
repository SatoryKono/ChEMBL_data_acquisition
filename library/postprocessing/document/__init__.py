"""Document postprocessing helpers and legacy compatibility exports."""

from __future__ import annotations

import importlib.util
from pathlib import Path

from . import steps

__all__ = [
    "DocumentFetchResult",
    "clean_doi_value",
    "fetch_chembl_documents",
    "fetch_crossref_metadata",
    "fetch_normalize_document",
    "fetch_openalex_metadata",
    "merge_document_metadata",
    "normalize_document_frame",
    "postprocess_export_file",
    "preprocess_document_export",
]

for name in __all__[:-2]:
    globals()[name] = getattr(steps, name)

_LEGACY_PATH = Path(__file__).resolve().parent.parent / "document.py"
_LEGACY_SPEC = importlib.util.spec_from_file_location(
    "library.postprocessing._document_legacy",
    _LEGACY_PATH,
)

if _LEGACY_SPEC is None or _LEGACY_SPEC.loader is None:  # pragma: no cover - defensive
    raise ImportError("Unable to load legacy document module")

_legacy_module = importlib.util.module_from_spec(_LEGACY_SPEC)
_LEGACY_SPEC.loader.exec_module(_legacy_module)

preprocess_document_export = getattr(_legacy_module, "preprocess_document_export")
postprocess_export_file = getattr(_legacy_module, "postprocess_export_file")
