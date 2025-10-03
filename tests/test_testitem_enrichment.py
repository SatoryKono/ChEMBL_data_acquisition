from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from library.pipelines.testitem import enrichment as testitem_enrichment
from library.config import Config


def _write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    frame = pd.DataFrame(rows)
    frame.to_csv(path, index=False)


def test_enrich_adds_flags_and_salt(tmp_path: Path, cfg: Config) -> None:
    df = pd.DataFrame(
        [
            {"molecule_chembl_id": "CHEMBL100"},
            {
                "molecule_chembl_id": "CHEMBL200",
                "parent_molecule_chembl_id": "CHEMBL0200P",
            },
            {
                "molecule_chembl_id": "CHEMBL300",
                "parent_molecule_chembl_id": "CHEMBL300",
            },
        ]
    )

    hierarchy_path = tmp_path / "hierarchy.csv"
    _write_csv(
        hierarchy_path,
        [
            {
                "molecule_chembl_id": "CHEMBL100",
                "parent_molecule_chembl_id": "CHEMBL0100P",
            },
            {
                "molecule_chembl_id": "CHEMBL200",
                "parent_molecule_chembl_id": "CHEMBL0200P",
            },
        ],
    )

    catalog_path = tmp_path / "catalog.csv"
    _write_csv(
        catalog_path,
        [
            {
                "molecule_chembl_id": "CHEMBL100",
                "natural_product": "Y",
                "prodrug": "",
                "polymer_flag": "0",
            },
            {
                "molecule_chembl_id": "CHEMBL0100P",
                "natural_product": "",
                "prodrug": "N",
                "polymer_flag": "",
            },
            {
                "molecule_chembl_id": "CHEMBL200",
                "natural_product": "N",
                "prodrug": "Y",
                "polymer_flag": "1",
            },
            {
                "molecule_chembl_id": "CHEMBL0200P",
                "natural_product": "",
                "prodrug": "Y",
                "polymer_flag": "0",
            },
        ],
    )

    cfg.testitem_molecule_enrichment.sources.molecule_hierarchy_path = hierarchy_path
    cfg.testitem_molecule_enrichment.sources.molecule_catalog_path = catalog_path

    result = testitem_enrichment.enrich(
        df,
        cfg=cfg.testitem_molecule_enrichment,
        io_cfg=cfg.io,
    )

    assert list(result["parent_molecule_chembl_id"]) == [
        "CHEMBL0100P",
        "CHEMBL0200P",
        "CHEMBL300",
    ]
    assert result["salt_chembl_id"].tolist()[:2] == ["CHEMBL100", "CHEMBL200"]
    assert pd.isna(result.loc[2, "salt_chembl_id"])
    assert result["natural_product"].dtype == "boolean"
    assert result["prodrug"].dtype == "boolean"
    assert bool(result.loc[0, "natural_product"]) is True
    assert bool(result.loc[0, "prodrug"]) is False
    assert bool(result.loc[1, "polymer_flag"]) is True
    assert pd.isna(result.loc[2, "natural_product"])


def test_enrich_logs_missing_and_inconsistent(
    tmp_path: Path, cfg: Config, monkeypatch: pytest.MonkeyPatch
) -> None:
    df = pd.DataFrame(
        [
            {
                "molecule_chembl_id": "CHEMBL500",
                "parent_molecule_chembl_id": "CHEMBL0500P",
            },
            {
                "molecule_chembl_id": "CHEMBL600",
                "parent_molecule_chembl_id": "CHEMBL0600P",
            },
        ]
    )

    hierarchy_path = tmp_path / "hierarchy.csv"
    _write_csv(
        hierarchy_path,
        [
            {
                "molecule_chembl_id": "CHEMBL500",
                "parent_molecule_chembl_id": "CHEMBL0500P",
            },
        ],
    )

    catalog_path = tmp_path / "catalog.csv"
    _write_csv(
        catalog_path,
        [
            {
                "molecule_chembl_id": "CHEMBL500",
                "natural_product": "Y",
                "prodrug": "N",
                "polymer_flag": "0",
            },
            {
                "molecule_chembl_id": "CHEMBL0500P",
                "natural_product": "N",
                "prodrug": "Y",
                "polymer_flag": "1",
            },
        ],
    )

    cfg.testitem_molecule_enrichment.sources.molecule_hierarchy_path = hierarchy_path
    cfg.testitem_molecule_enrichment.sources.molecule_catalog_path = catalog_path

    events: list[tuple[str, dict[str, object]]] = []

    def capture(event: str, *args: object, **kwargs: object) -> None:
        events.append((event, kwargs))

    monkeypatch.setattr(testitem_enrichment.logger, "warning", capture)

    testitem_enrichment.enrich(
        df,
        cfg=cfg.testitem_molecule_enrichment,
        io_cfg=cfg.io,
    )

    event_names = [event for event, _ in events]
    assert "testitem_enrichment_inconsistent_flag" in event_names
    assert "testitem_enrichment_missing_child_flags" in event_names
    assert "testitem_enrichment_missing_parent_flags" in event_names


def test_enrich_respects_configuration_options(tmp_path: Path, cfg: Config) -> None:
    df = pd.DataFrame(
        [
            {
                "molecule_chembl_id": "CHEMBL800",
                "parent_molecule_chembl_id": "CHEMBL0800P",
            },
            {
                "molecule_chembl_id": "CHEMBL810",
                "parent_molecule_chembl_id": "CHEMBL0810P",
            },
            {
                "molecule_chembl_id": "CHEMBL900",
                "parent_molecule_chembl_id": "CHEMBL900",
            },
        ]
    )

    hierarchy_path = tmp_path / "hierarchy.csv"
    _write_csv(
        hierarchy_path,
        [
            {
                "molecule_chembl_id": "CHEMBL800",
                "parent_molecule_chembl_id": "CHEMBL0800P",
            },
            {
                "molecule_chembl_id": "CHEMBL810",
                "parent_molecule_chembl_id": "CHEMBL0810P",
            },
        ],
    )

    catalog_path = tmp_path / "catalog.csv"
    _write_csv(
        catalog_path,
        [
            {
                "molecule_chembl_id": "CHEMBL800",
                "natural_product": "Y",
                "prodrug": "N",
                "polymer_flag": "0",
            },
            {
                "molecule_chembl_id": "CHEMBL0800P",
                "natural_product": "N",
                "prodrug": "Y",
                "polymer_flag": "1",
            },
            {
                "molecule_chembl_id": "CHEMBL0810P",
                "natural_product": "N",
                "prodrug": "N",
                "polymer_flag": "0",
            },
        ],
    )

    cfg.testitem_molecule_enrichment.output.salt_as_null_when_absent = False
    cfg.testitem_molecule_enrichment.flags.parent_fallback = False
    cfg.testitem_molecule_enrichment.flags.coerce_to_bool = False
    cfg.testitem_molecule_enrichment.sources.molecule_hierarchy_path = hierarchy_path
    cfg.testitem_molecule_enrichment.sources.molecule_catalog_path = catalog_path

    result = testitem_enrichment.enrich(
        df,
        cfg=cfg.testitem_molecule_enrichment,
        io_cfg=cfg.io,
    )

    assert result["natural_product"].dtype == "string"
    assert result.loc[0, "natural_product"] == "Y"
    assert pd.isna(result.loc[1, "natural_product"])
    assert result.loc[2, "salt_chembl_id"] == "-"
