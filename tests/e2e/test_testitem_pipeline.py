from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from library.common.csv_utils import sha256_file, write_csv_deterministic
from library.pipelines.testitem import cli, enrichment
from library.schemas.normalize import normalize_testitems


@pytest.mark.e2e
def test_testitem_pipeline_e2e__deterministic_output(
    tmp_path: Path,
    sample_input_csv: Path,
    cfg,
    snapshot_resource: Path,
) -> None:
    cfg.testitem_molecule_enrichment.enable = True
    cfg.testitem_molecule_enrichment.sources.molecule_catalog_path = (
        snapshot_resource / "molecule_catalog.csv"
    )
    cfg.testitem_molecule_enrichment.sources.molecule_hierarchy_path = (
        snapshot_resource / "molecule_hierarchy.csv"
    )

    rc, read_result = cli.read_input_ids(
        sample_input_csv,
        column="molecule_chembl_id",
        io_cfg=cfg.io,
        limit=None,
        offset=0,
    )
    assert rc == 0
    assert read_result is not None

    requested_ids = list(read_result.ids_iter)

    raw_frame = pd.DataFrame(
        {
            "molecule_chembl_id": requested_ids,
            "parent_molecule_chembl_id": [pd.NA, "CHEMBL2"],
            "relation": ["<", ">"],
            "units": ["1 μM", ""],
            "natural_product": ["yes", ""],
        }
    )

    normalised = normalize_testitems(raw_frame)
    enriched = enrichment.enrich(normalised, cfg=cfg.testitem_molecule_enrichment, io_cfg=cfg.io)
    final = cli.integrate_missing_identifiers(
        enriched,
        missing_ids=["CHEMBL999"],
        requested_ids=requested_ids + ["CHEMBL999"],
    )

    output_path = tmp_path / "output.csv"
    write_csv_deterministic(final, output_path, key_cols=sorted(final.columns))
    first_hash = sha256_file(output_path)

    repeat_path = tmp_path / "output_repeat.csv"
    write_csv_deterministic(final, repeat_path, key_cols=sorted(final.columns))
    second_hash = sha256_file(repeat_path)

    assert first_hash == second_hash
    reloaded = pd.read_csv(output_path)
    assert list(reloaded["molecule_chembl_id"]) == ["CHEMBL1", "CHEMBL2", "CHEMBL999"]
    assert pd.isna(reloaded.loc[2, "parent_molecule_chembl_id"])
    assert pd.isna(reloaded.loc[2, "natural_product"])
