from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from library.common.csv_utils import sha256_file, write_csv_deterministic
from library.pipelines.testitem import cli, enrichment
from library.pipelines.testitem.catalog import ParentLookupStats
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
    enriched = enrichment.enrich(
        normalised, cfg=cfg.testitem_molecule_enrichment, io_cfg=cfg.io
    )
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


@pytest.mark.e2e
@pytest.mark.pipeline_scenario("idempotence")
def test_testitem_pipeline_e2e__finalize_stage_idempotent(
    tmp_path: Path, sample_input_csv: Path, cfg
) -> None:
    cfg.system.doc_quality.enable = False

    chunk = pd.DataFrame(
        {
            "molecule_chembl_id": pd.Series(["CHEMBL1", "CHEMBL2"], dtype="string"),
            "parent_molecule_chembl_id": pd.Series(["CHEMBL10", pd.NA], dtype="string"),
            "natural_product": pd.Series([True, False], dtype="boolean"),
            "salt_chembl_id": pd.Series(["CHEMBL1", "CHEMBL2"], dtype="string"),
        }
    )

    stats = ParentLookupStats(
        source="lookup",
        missing=0,
        unique=2,
        attached=2,
        uncovered=0,
    )

    def supplier() -> ParentLookupStats:
        return stats

    output_path = tmp_path / "finalized.csv"

    first_result = cli.finalize_output(
        [chunk],
        cfg=cfg,
        output=output_path,
        parent_stats_supplier=supplier,
        input_csv=sample_input_csv,
    )
    if isinstance(first_result, tuple):
        first_exit, artifacts = first_result
        dataset_path = Path(artifacts.dataset)
    else:
        first_exit = first_result
        dataset_path = output_path
    assert first_exit == 0
    first_hash = sha256_file(dataset_path)

    second_result = cli.finalize_output(
        [chunk.copy()],
        cfg=cfg,
        output=output_path,
        parent_stats_supplier=supplier,
        input_csv=sample_input_csv,
    )
    if isinstance(second_result, tuple):
        second_exit, second_artifacts = second_result
        assert Path(second_artifacts.dataset) == dataset_path
    else:
        second_exit = second_result
    assert second_exit == 0
    second_hash = sha256_file(dataset_path)

    assert first_hash == second_hash

    final_frame = pd.read_csv(dataset_path)
    assert list(final_frame["molecule_chembl_id"]) == ["CHEMBL1", "CHEMBL2"]
