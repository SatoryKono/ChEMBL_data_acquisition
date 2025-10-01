"""Regression tests for molecule catalog cache parsing."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from library.config import MoleculeCatalogCfg
from library import molecule_catalog as mc
from library.utils import atomic as atomic_utils


def _cfg(cache_path: Path) -> MoleculeCatalogCfg:
    return MoleculeCatalogCfg().model_copy(
        update={
            "cache_path": cache_path,
            "sqlite_path": cache_path.with_suffix(".sqlite"),
        }
    )


def test_read_cache_csv_normalises_ids(tmp_path: Path) -> None:
    path = tmp_path / "catalog.csv"
    path.write_text(
        "\n".join(
            [
                "molecule_chembl_id,parent_molecule_chembl_id",
                " chembl1 , CHEMBL999 ",
                "CHEMBL2,",  # missing parent skipped
                ",CHEMBL3",  # missing child skipped
            ]
        ),
        encoding="utf-8",
    )

    result = mc._read_cache(path, _cfg(path))

    assert result == {"CHEMBL1": "CHEMBL999"}


def test_read_cache_invalid_json_logs_warning(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    path = tmp_path / "catalog.json"
    path.write_text("not-json", encoding="utf-8")

    events: list[tuple[str, dict[str, object]]] = []

    def capture(event: str, *args, extra: dict[str, object] | None = None, **kwargs):
        events.append((event, extra or {}))

    monkeypatch.setattr(mc.logger, "warning", capture)

    result = mc._read_cache(path, _cfg(path))

    assert result == {}
    assert any(event == "invalid_catalog_cache" for event, _ in events)


def test_write_cache_crash_leaves_existing_file(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    path = tmp_path / "catalog.json"
    original = {"CHEMBL1": "CHEMBL2"}
    path.write_text(json.dumps(original), encoding="utf-8")
    cfg = _cfg(path)

    def fail_replace(src: str | Path, dst: str | Path) -> None:
        raise OSError("boom")

    monkeypatch.setattr(atomic_utils.os, "replace", fail_replace)

    with pytest.raises(OSError):
        mc._write_cache_from_items([("CHEMBL3", "CHEMBL4")], cfg)

    assert json.loads(path.read_text(encoding="utf-8")) == original
    assert not list(path.parent.glob(f".{path.name}.*.tmp"))
