"""Tests for :mod:`library.document_defaults`."""

from __future__ import annotations

import pytest

from library.document_defaults import DocumentCfg


@pytest.mark.unit
def test_document_cfg__accepts_valid_values() -> None:
    cfg = DocumentCfg(
        column="document_chembl_id",
        batch_size=10,
        chunk_size=25,
        sleep=0.5,
        workers=2,
        timeout=30.0,
    )

    assert cfg.column == "document_chembl_id"
    assert cfg.batch_size == 10
    assert cfg.chunk_size == 25
    assert cfg.sleep == 0.5
    assert cfg.workers == 2
    assert cfg.timeout == 30.0


@pytest.mark.unit
@pytest.mark.parametrize(
    ("field", "value", "exception"),
    [
        ("column", "", ValueError),
        ("batch_size", 0, ValueError),
        ("batch_size", 1.5, TypeError),
        ("chunk_size", 0, ValueError),
        ("chunk_size", True, TypeError),
        ("sleep", -0.1, ValueError),
        ("sleep", True, TypeError),
        ("workers", 0, ValueError),
        ("workers", False, TypeError),
        ("timeout", 0.0, ValueError),
        ("timeout", True, TypeError),
    ],
)
def test_document_cfg__invalid_boundaries(
    field: str, value: object, exception: type[Exception]
) -> None:
    kwargs: dict[str, object] = {
        "column": "valid",
        "batch_size": 10,
        "chunk_size": 20,
        "sleep": 1.0,
        "workers": 1,
        "timeout": 5.0,
    }
    kwargs[field] = value

    with pytest.raises(exception):
        DocumentCfg(**kwargs)
