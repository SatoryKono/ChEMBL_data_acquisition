"""Deprecated compatibility layer for :mod:`library.pipelines.testitem`."""

from __future__ import annotations

import warnings

from library.pipelines.testitem import core as _core
from library.pipelines.testitem import (
    ChemblClient,  # noqa: F401
    TESTITEM_PUBCHEM_COLUMNS,  # noqa: F401
    analyze_table_quality,  # noqa: F401
    file_sha256,  # noqa: F401
    logger,  # noqa: F401
    pc,  # noqa: F401
    pl,  # noqa: F401
    run_pipeline,  # noqa: F401
    run_testitem_pipeline,  # noqa: F401
    TestitemPipelineOptions,  # noqa: F401
    validate_testitems,  # noqa: F401
    write_csv_deterministic,  # noqa: F401
    write_meta_yaml,  # noqa: F401
)
from library.pipelines.testitem import enrichment as testitem_enrichment  # noqa: F401

_DEPRECATION_MESSAGE = (
    "library.testitem_pipeline is deprecated and will be removed in a future "
    "release; use library.pipelines.testitem instead."
)

warnings.warn(_DEPRECATION_MESSAGE, DeprecationWarning, stacklevel=2)

for _name in _core.__all__:
    globals()[_name] = getattr(_core, _name)
del _name

__all__ = list(
    dict.fromkeys(
        [
            *_core.__all__,
            "ChemblClient",
            "TESTITEM_PUBCHEM_COLUMNS",
            "analyze_table_quality",
            "file_sha256",
            "logger",
            "pc",
            "pl",
            "run_pipeline",
            "testitem_enrichment",
            "validate_testitems",
            "write_csv_deterministic",
            "write_meta_yaml",
        ]
    )
)
