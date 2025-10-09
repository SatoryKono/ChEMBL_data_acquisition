from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock

import pandas as pd
import pytest

from library.config import MoleculeCatalogCfg
from library.pipelines.testitem import catalog


def _stub_parent_catalog_calls(monkeypatch: pytest.MonkeyPatch) -> None:
    """Replace parent catalog I/O helpers with deterministic stubs."""

    monkeypatch.setattr(catalog, "query_parent_catalog", lambda *args, **kwargs: {})
    monkeypatch.setattr(catalog, "load_parent_catalog", lambda *args, **kwargs: {})
    monkeypatch.setattr(catalog, "update_parent_catalog_cache", lambda *args, **kwargs: None)
    monkeypatch.setattr(catalog, "write_parent_catalog_cache", lambda *args, **kwargs: None)
    monkeypatch.setattr(
        catalog.molecule_catalog,
        "fetch_parent_catalog_for",
        lambda *args, **kwargs: {},
    )


def _build_catalog_cfg(tmp_path: Path) -> MoleculeCatalogCfg:
    return MoleculeCatalogCfg(
        cache_path=tmp_path / "parent_catalog.json",
        sqlite_path=tmp_path / "parent_catalog.sqlite",
    )


def _capture_events(
    monkeypatch: pytest.MonkeyPatch,
) -> list[tuple[str, str, dict[str, object]]]:
    events: list[tuple[str, str, dict[str, object]]] = []

    def _capture(level: str):
        def _record(event: str, **fields: object) -> None:
            events.append((level, event, fields))

        return _record

    monkeypatch.setattr(catalog.logger, "warning", _capture("warning"))
    monkeypatch.setattr(catalog.logger, "info", _capture("info"))
    return events


def _run_parent_lookup(
    child_ids: list[str],
    *,
    cfg: MoleculeCatalogCfg,
    api_cfg,
    monkeypatch: pytest.MonkeyPatch,
) -> list[tuple[str, dict[str, object]]]:
    _stub_parent_catalog_calls(monkeypatch)
    events = _capture_events(monkeypatch)

    client = MagicMock(spec=catalog.ChemblClient)
    frame = pd.DataFrame({"molecule_chembl_id": child_ids})

    catalog.attach_parent_molecule_ids(
        frame,
        client=client,
        api_cfg=api_cfg,
        catalog_cfg=cfg,
        timeout=1.0,
        catalog=None,
    )

    return events


def _extract_event(
    events: list[tuple[str, str, dict[str, object]]], name: str
) -> tuple[str, dict[str, object]]:
    level, fields = next(
        (level, fields) for level, event, fields in events if event == name
    )
    return level, fields


@pytest.mark.integration
@pytest.mark.parametrize(
    "child_ids, truncated",
    [
        (
            [
                "CHEMBL0005",
                "CHEMBL0002",
                "CHEMBL0003",
            ],
            False,
        ),
        (
            [
                f"CHEMBL{value:07d}"
                for value in range(1024, 998, -1)
            ],
            True,
        ),
    ],
    ids=["no-truncation", "truncated"],
)
def test_attach_parent_molecule_ids__summarises_identifier_payload(
    child_ids: list[str],
    truncated: bool,
    tmp_path,
    cfg,
    monkeypatch,
) -> None:
    catalog_cfg = _build_catalog_cfg(tmp_path)
    events = _run_parent_lookup(
        child_ids,
        cfg=catalog_cfg,
        api_cfg=cfg.api,
        monkeypatch=monkeypatch,
    )

    expected_identifiers = sorted(child_ids)
    if truncated:
        assert len(expected_identifiers) > 20
        expected_identifiers = expected_identifiers[:20]

    skip_level, skip_event = _extract_event(
        events, "parent_lookup_full_sync_skipped_parentless"
    )
    missing_level, missing_event = _extract_event(
        events, "parent_lookup_missing_parents"
    )

    assert skip_level == "info"
    assert missing_level == "info"

    warning_events = [event for event in events if event[0] == "warning"]
    assert not warning_events

    for payload in (skip_event, missing_event):
        assert payload["count"] == len(child_ids)
        assert payload["identifiers"] == expected_identifiers
        assert payload["truncated"] is truncated
