"""Tests for resilient loading of the document schema declaration."""

from __future__ import annotations

from pathlib import Path

import pytest

from library.schemas import document_spec


@pytest.mark.unit
def test_load_document_declaration__packaged_resource(monkeypatch, tmp_path) -> None:
    monkeypatch.setattr(document_spec, "SCHEMA_DIR", tmp_path)

    declaration = document_spec.load_document_declaration()

    assert declaration.ordered_columns
    assert any(group.columns for group in declaration.groups)


@pytest.mark.unit
def test_load_document_declaration__missing_default_file(monkeypatch, tmp_path) -> None:
    schema_dir = tmp_path / "schema"
    schema_dir.mkdir()
    expected_path = schema_dir / "document.yaml"

    monkeypatch.setattr(document_spec, "SCHEMA_DIR", schema_dir)

    original_exists = Path.exists
    original_open = Path.open

    def fake_exists(self: Path) -> bool:  # pragma: no cover - thin shim
        if self == expected_path:
            return True
        return original_exists(self)

    def fake_open(self: Path, *args, **kwargs):  # pragma: no cover - thin shim
        if self == expected_path:
            raise FileNotFoundError(self)
        return original_open(self, *args, **kwargs)

    monkeypatch.setattr(Path, "exists", fake_exists)
    monkeypatch.setattr(Path, "open", fake_open)

    declaration = document_spec.load_document_declaration()

    assert declaration.ordered_columns
    assert any(group.columns for group in declaration.groups)


@pytest.mark.unit
def test_load_document_declaration__explicit_path(tmp_path) -> None:
    schema_dir = tmp_path / "schema"
    schema_dir.mkdir()
    schema_path = schema_dir / "document.yaml"
    schema_path.write_text("groups: []\n", encoding="utf-8")

    declaration = document_spec.load_document_declaration(schema_dir)

    assert declaration.groups == tuple()
    assert declaration.ordered_columns == tuple()
