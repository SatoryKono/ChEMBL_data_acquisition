from __future__ import annotations

import importlib

import library.pipelines.testitem as pipeline


def test_pipeline_exports_are_available() -> None:
    """The pipeline module exposes documented helpers."""

    expected_attrs = {
        "_log_missing_identifier_summary",
        "_merge_pubchem_properties",
    }

    # Sanity check for cache invalidation issues.
    importlib.reload(pipeline)

    for attr in expected_attrs:
        assert hasattr(pipeline, attr), attr

    # Public API should not contain duplicates and must be resolvable.
    assert len(pipeline.__all__) == len(set(pipeline.__all__))
    for attr in pipeline.__all__:
        assert hasattr(pipeline, attr), attr
