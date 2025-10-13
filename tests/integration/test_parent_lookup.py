from __future__ import annotations

import logging
import shutil
from pathlib import Path
from unittest.mock import MagicMock

import pandas as pd
import pytest

from library.config import MoleculeCatalogCfg
from library.pipelines.testitem import catalog
from scripts import get_testitem_data


@pytest.fixture()
def hierarchy_lookup_parity(
    cfg,
    tmp_path: Path,
    caplog: pytest.LogCaptureFixture,
):
    """Run CLI and pipeline lookup loaders and ensure parity."""

    from library.pipelines.testitem import catalog as pipeline_catalog
    from scripts import get_testitem_data

    default_source = Path(__file__).resolve().parents[1] / "resources" / "molecule_hierarchy.csv"

    def _run(
        *,
        source_path: Path | None = None,
        encoding: str | None = None,
        delimiter: str | None = None,
    ) -> dict[str, object]:
        csv_source = Path(source_path or default_source)
        resolved_path = tmp_path / csv_source.name
        shutil.copy(csv_source, resolved_path)

        caplog.set_level(logging.INFO, logger="chembl")

        caplog.clear()
        cli_mapping = get_testitem_data.load_molecule_hierarchy_lookup(
            resolved_path,
            io_cfg=cfg.io,
            encoding=encoding,
            delimiter=delimiter,
        )
        cli_messages = [
            record.getMessage()
            for record in caplog.records
            if "molecule_hierarchy_lookup_loaded" in record.getMessage()
        ]

        caplog.clear()
        pipeline_mapping = catalog.load_molecule_hierarchy_lookup(
            resolved_path,
            io_cfg=cfg.io,
            encoding=encoding,
            delimiter=delimiter,
        )
        pipeline_messages = [
            record.getMessage()
            for record in caplog.records
            if "molecule_hierarchy_lookup_loaded" in record.getMessage()
        ]

        assert cli_mapping == pipeline_mapping
        assert cli_messages == pipeline_messages

        return {
            "mapping": cli_mapping,
            "messages": cli_messages,
            "path": str(resolved_path),
        }

    return _run


@pytest.mark.integration
def test_load_molecule_hierarchy_lookup__cli_pipeline_parity(
    hierarchy_lookup_parity,
):
    """CLI compatibility wrapper should mirror pipeline behaviour."""

    result = hierarchy_lookup_parity()

    expected_message = (
        "molecule_hierarchy_lookup_loaded "
        f"path='{result['path']}' rows=2 rps=None status=None"
    )

    assert result["messages"] == [expected_message]


def _stub_parent_catalog_calls(monkeypatch: pytest.MonkeyPatch) -> None:
    """Replace parent catalog I/O helpers with deterministic stubs."""

    monkeypatch.setattr(catalog, "query_parent_catalog", lambda *args, **kwargs: {})
    monkeypatch.setattr(catalog, "load_parent_catalog", lambda *args, **kwargs: {})
    monkeypatch.setattr(
        catalog, "update_parent_catalog_cache", lambda *args, **kwargs: None
    )
    monkeypatch.setattr(
        catalog, "write_parent_catalog_cache", lambda *args, **kwargs: None
    )
    monkeypatch.setattr(
        catalog.molecule_catalog,
        "fetch_parent_catalog_for",
        lambda *args, **kwargs: {},
    )
    monkeypatch.setattr(
        catalog.molecule_catalog,
        "fetch_parent_catalog",
        lambda *args, **kwargs: {},
    )
    monkeypatch.setattr(
        catalog.molecule_catalog,
        "load_parent_catalog",
        lambda *args, **kwargs: {},
    )
    monkeypatch.setattr(
        catalog.molecule_catalog,
        "query_parent_catalog",
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


@pytest.mark.integration
@pytest.mark.parametrize(
    "child_ids, truncated, catalog_filters, expected_level, expected_severity",
    [
        (
            [
                "CHEMBL0005",
                "CHEMBL0002",
                "CHEMBL0003",
            ],
            False,
            None,
            "info",
            "info",
        ),
        (
            [f"CHEMBL{value:07d}" for value in range(1024, 998, -1)],
            True,
            None,
            "info",
            "info",
        ),
        (
            [
                "CHEMBL0005",
                "CHEMBL0002",
                "CHEMBL0003",
            ],
            False,
            {"parent_molecule_chembl_id__isnull": "true"},
            "warning",
            "warning",
        ),
    ],
    ids=["no-truncation", "truncated", "include-parentless"],
)
def test_attach_parent_molecule_ids__summarises_identifier_payload(
    child_ids: list[str],
    truncated: bool,
    catalog_filters: dict[str, str] | None,
    expected_level: str,
    expected_severity: str,
    tmp_path,
    cfg,
    monkeypatch,
) -> None:
    catalog_cfg = _build_catalog_cfg(tmp_path)
    if catalog_filters is not None:
        catalog_cfg.filters = catalog_filters
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

    skip_events = [
        (level, fields)
        for level, event, fields in events
        if event == "parent_lookup_full_sync_skipped_parentless"
    ]
    missing_events = [
        (level, fields)
        for level, event, fields in events
        if event == "parent_lookup_missing_parents"
    ]

    assert len(missing_events) == 1
    missing_level, missing_event = missing_events[0]

    assert missing_level == expected_level
    assert missing_event["severity"] == expected_severity
    assert missing_event["count"] == len(child_ids)
    assert missing_event["identifiers"] == expected_identifiers
    assert missing_event["truncated"] is truncated

    if expected_level == "info":
        assert len(skip_events) == 1
        skip_level, skip_event = skip_events[0]
        assert skip_level == "info"
        assert skip_event["count"] == len(child_ids)
        assert skip_event["identifiers"] == expected_identifiers
        assert skip_event["truncated"] is truncated
    else:
        assert not skip_events

    warning_events = [event for event in events if event[0] == "warning"]
    if expected_level == "warning":
        assert warning_events == [
            ("warning", "parent_lookup_missing_parents", missing_event)
        ]
    else:
        assert not warning_events
