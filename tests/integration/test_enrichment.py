from pathlib import Path

import sys
from pathlib import Path

import pandas as pd
import pytest

from library.pipelines.testitem import enrichment


@pytest.mark.integration
def test_enrich__attaches_flags_and_parent(
    tmp_path: Path, cfg, snapshot_resource, monkeypatch
) -> None:
    cfg.testitem_molecule_enrichment.enable = True
    cfg.testitem_molecule_enrichment.sources.molecule_catalog_path = (
        snapshot_resource / "molecule_catalog.csv"
    )
    cfg.testitem_molecule_enrichment.sources.molecule_hierarchy_path = (
        snapshot_resource / "molecule_hierarchy.csv"
    )

    events: list[tuple[str, dict[str, object]]] = []

    def _capture(event: str, **fields: object) -> None:
        events.append((event, fields))

    monkeypatch.setattr(enrichment.logger, "warning", _capture)

    frame = pd.DataFrame(
        {
            "molecule_chembl_id": ["CHEMBL1", "CHEMBL2", "CHEMBL3"],
            "parent_molecule_chembl_id": [pd.NA, pd.NA, pd.NA],
        }
    )

    enriched = enrichment.enrich(
        frame,
        cfg=cfg.testitem_molecule_enrichment,
        io_cfg=cfg.io,
    )

    assert list(enriched["parent_molecule_chembl_id"]) == ["CHEMBL10", pd.NA, "CHEMBL30"]
    assert enriched["salt_chembl_id"].tolist() == ["CHEMBL1", pd.NA, "CHEMBL3"]
    assert enriched["natural_product"].dtype == "boolean"
    assert enriched["natural_product"].tolist() == [True, pd.NA, False]
    missing_child_events = [
        fields for event, fields in events if event == "testitem_enrichment_missing_child_flags"
    ]
    assert missing_child_events == [
        {"count": 2, "identifiers": ["CHEMBL2", "CHEMBL3"]}
    ]


@pytest.mark.integration
def test_enrich__handles_unknown_flag_values(tmp_path: Path, cfg, snapshot_resource, caplog) -> None:
    cfg.testitem_molecule_enrichment.enable = True
    cfg.testitem_molecule_enrichment.sources.molecule_catalog_path = snapshot_resource / "molecule_catalog.csv"
    cfg.testitem_molecule_enrichment.sources.molecule_hierarchy_path = snapshot_resource / "molecule_hierarchy.csv"
    cfg.testitem_molecule_enrichment.flags.coerce_to_bool = True

    frame = pd.DataFrame(
        {
            "molecule_chembl_id": ["CHEMBL1"],
            "parent_molecule_chembl_id": [pd.NA],
        }
    )

    caplog.set_level("WARNING")

    enriched = enrichment.enrich(
        frame,
        cfg=cfg.testitem_molecule_enrichment,
        io_cfg=cfg.io,
    )

    assert enriched["natural_product"].tolist() == [True]
    assert any("unknown_flag_values" in record.message for record in caplog.records) is False


@pytest.mark.integration
def test_enrich__missing_sources_gracefully_fills_columns(
    tmp_path: Path, cfg, monkeypatch
) -> None:
    cfg.testitem_molecule_enrichment.enable = True
    cfg.testitem_molecule_enrichment.sources.molecule_catalog_path = tmp_path / "missing_catalog.csv"
    cfg.testitem_molecule_enrichment.sources.molecule_hierarchy_path = tmp_path / "missing_hierarchy.csv"

    frame = pd.DataFrame({"molecule_chembl_id": ["CHEMBL1"]})
    events: list[str] = []

    def fake_read_csv(*_args: object, **_kwargs: object) -> pd.DataFrame:
        raise FileNotFoundError("missing file")

    def fake_warning(event: str, **_fields: object) -> None:
        events.append(event)

    monkeypatch.setattr(enrichment.io, "read_csv", fake_read_csv)
    monkeypatch.setattr(enrichment.logger, "warning", fake_warning)
    enriched = enrichment.enrich(
        frame,
        cfg=cfg.testitem_molecule_enrichment,
        io_cfg=cfg.io,
    )

    assert list(enriched.columns) == [
        "molecule_chembl_id",
        "parent_molecule_chembl_id",
        "salt_chembl_id",
        "natural_product",
        "prodrug",
        "polymer_flag",
    ]
    assert enriched["parent_molecule_chembl_id"].tolist() == [pd.NA]
    assert enriched["natural_product"].dtype == "boolean"
    assert events.count("testitem_enrichment_missing_hierarchy") == 1
    assert events.count("testitem_enrichment_missing_catalog") == 1
