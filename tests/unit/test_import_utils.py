"""Tests for :mod:`library.postprocessing.pipeline.common.import_utils`."""

from __future__ import annotations

from pathlib import Path

import pytest

from library.postprocessing.pipeline.common.import_utils import (
    ImportResolutionError,
    import_by_path,
    resolve_dotted_path,
)


@pytest.mark.unit
def test_import_by_path__resolves_callable() -> None:
    sqrt = import_by_path("math:sqrt")

    assert sqrt(9) == 3


@pytest.mark.unit
def test_import_by_path__missing_module() -> None:
    with pytest.raises(ImportResolutionError) as excinfo:
        import_by_path("not_a_real_module:function")

    assert "Cannot import module" in str(excinfo.value)


@pytest.mark.unit
def test_import_by_path__missing_attribute() -> None:
    with pytest.raises(ImportResolutionError) as excinfo:
        import_by_path("math:not_real")

    assert "has no attribute" in str(excinfo.value)


@pytest.mark.unit
def test_import_by_path__non_callable_object_rejected() -> None:
    with pytest.raises(ImportResolutionError) as excinfo:
        import_by_path("math:pi")

    assert "must be callable" in str(excinfo.value)


@pytest.mark.unit
def test_import_by_path__custom_expected_type() -> None:
    value = import_by_path("math:pi", expected_type=(int, float))

    assert isinstance(value, float)


@pytest.mark.unit
def test_resolve_dotted_path__supports_expected_type_argument() -> None:
    result = resolve_dotted_path("pathlib.Path", expected_type=type)

    assert result is Path


@pytest.mark.unit
def test_resolve_dotted_path__raises_for_unmatched_expected_type() -> None:
    with pytest.raises(ImportResolutionError):
        resolve_dotted_path("pathlib.Path", expected_type=int)
