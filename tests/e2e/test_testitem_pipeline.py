from __future__ import annotations

from pathlib import Path
import tracemalloc

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


@pytest.mark.e2e
def test_testitem_pipeline_e2e__missing_identifier_stats_fields(
    tmp_path: Path,
    sample_input_csv: Path,
    cfg,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    cfg.system.doc_quality.enable = False

    chunk = pd.DataFrame(
        {
            "molecule_chembl_id": pd.Series(["CHEMBL1"], dtype="string"),
            "parent_molecule_chembl_id": pd.Series([pd.NA], dtype="string"),
        }
    )

    stats = ParentLookupStats(
        source="lookup",
        missing=0,
        unique=1,
        attached=1,
        uncovered=0,
    )

    def supplier() -> ParentLookupStats:
        return stats

    info_events: list[tuple[str, dict[str, object]]] = []

    def capture_info(event: str, *args: object, **fields: object) -> None:
        if args:
            return
        info_events.append((event, fields))

    monkeypatch.setattr(cli.logger, "info", capture_info)

    missing_ids = ["CHEMBL999"]

    exit_code, _ = cli.finalize_output(
        [chunk],
        cfg=cfg,
        output=tmp_path / "missing.csv",
        parent_stats_supplier=supplier,
        input_csv=sample_input_csv,
        missing_ids=missing_ids,
        emit_legacy_artifacts=False,
    )

    assert exit_code == 0

    stats_events = [fields for event, fields in info_events if event == "testitem_stats"]
    assert stats_events, "expected testitem_stats event to be emitted"
    payload = stats_events[-1]

    assert payload["missing_molecule_ids"] == missing_ids
    assert payload["missing_ids_sample"] == missing_ids
    assert payload["missing_molecule_ids_total"] == len(missing_ids)
    assert payload["missing_molecule_ids_count"] == len(missing_ids)
    assert payload["missing_molecule_ids_truncated"] is False


@pytest.mark.e2e
@pytest.mark.pipeline_scenario("missing_identifiers_streaming")
def test_testitem_pipeline_e2e__placeholder_generation_streaming_memory() -> None:
    missing_total = 6500
    missing_ids = [f"CHEMBL{idx:06d}" for idx in range(10_000, 10_000 + missing_total)]

    chunk_size = cli.MISSING_IDENTIFIER_PLACEHOLDER_CHUNK_SIZE
    tracemalloc.start()
    chunk_lengths: list[int] = []
    max_chunk_memory = 0
    peak = 0
    try:
        for chunk in cli.generate_missing_identifier_placeholders(
            missing_ids, columns=("molecule_chembl_id",)
        ):
            chunk_lengths.append(len(chunk))
            max_chunk_memory = max(
                max_chunk_memory, int(chunk.memory_usage(deep=True).sum())
            )
        _current, peak = tracemalloc.get_traced_memory()
    finally:
        tracemalloc.stop()

    expected_full_chunks, remainder = divmod(missing_total, chunk_size)
    expected_lengths = [chunk_size] * expected_full_chunks
    if remainder:
        expected_lengths.append(remainder)

    assert chunk_lengths == expected_lengths
    # Individual placeholder chunks should remain well below 2 MB even for
    # identifier-heavy datasets.
    assert max_chunk_memory <= 2_000_000
    # ``tracemalloc`` peak memory reflects the largest allocation during the
    # streaming loop.  The chunk-based generator should keep it within a few
    # megabytes for the tested payload size.
    assert peak <= 5_000_000
