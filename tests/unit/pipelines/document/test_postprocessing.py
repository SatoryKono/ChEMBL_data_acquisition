"""Unit tests for document post-processing helpers."""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest

from library.pipelines.document import postprocessing as module


def _write_reference_csv(path: str, *, doc_id: str, classification: int) -> None:
    frame = pd.DataFrame(
        {
            "document_chembl_id": [doc_id],
            "classification": [classification],
            "document_contains_external_links": ["true"],
            "is_experimental_doc": ["false"],
        }
    )
    frame.to_csv(path, index=False, encoding=module.CP1252_ENCODING, sep=module.CSV_DELIMITER)


@pytest.mark.unit
def test_resolve_reference__prefers_cli_override(tmp_path: Path) -> None:
    resources_dir = tmp_path / "resources"
    (resources_dir / "_document").mkdir(parents=True)
    resource_path = resources_dir / "_document" / "document.csv"
    _write_reference_csv(resource_path, doc_id="RESOURCE_REF", classification=0)

    cli_path = tmp_path / "override.csv"
    _write_reference_csv(cli_path, doc_id="CLI_REF", classification=1)

    resources = SimpleNamespace(dictionary_dir=resources_dir)

    resolved = module._resolve_reference(None, cli_path, resources)

    assert list(resolved.columns) == list(module.REFERENCE_REQUIRED_COLUMNS)
    assert resolved.loc[0, "document_chembl_id"] == "CLI_REF"
    assert bool(resolved.loc[0, "doctype_review"]) is True


@pytest.mark.unit
def test_resolve_reference__uses_resources_dictionary(tmp_path: Path) -> None:
    resources_dir = tmp_path / "resources"
    (resources_dir / "_document").mkdir(parents=True)
    resource_path = resources_dir / "_document" / "document.csv"
    _write_reference_csv(resource_path, doc_id="RESOURCE_REF", classification=1)

    resources = SimpleNamespace(dictionary_dir=resources_dir)

    resolved = module._resolve_reference(None, None, resources)

    assert resolved.loc[0, "document_chembl_id"] == "RESOURCE_REF"
    assert bool(resolved.loc[0, "doctype_review"]) is True
    assert resolved.loc[0, "document_contains_external_links"] is True
    assert resolved.loc[0, "is_experimental_doc"] is False


@pytest.mark.unit
def test_resolve_reference__missing_reference_reports_checked_paths(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    resources_dir = tmp_path / "resources"
    resources_dir.mkdir(parents=True)
    resources = SimpleNamespace(dictionary_dir=resources_dir)

    fallback_path = tmp_path / "missing_default.csv"
    monkeypatch.setattr(module, "DEFAULT_REF_DOCUMENT_PATH", fallback_path)

    with pytest.raises(FileNotFoundError) as excinfo:
        module._resolve_reference(None, None, resources)

    message = str(excinfo.value)
    assert "data/input" not in message
    expected_resource_path = resources_dir / module.REFERENCE_RESOURCE_SUBPATH
    assert str(expected_resource_path) in message
    assert str(fallback_path) in message
