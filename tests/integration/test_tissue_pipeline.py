"""Integration tests for the tissue pipeline using recorded fixtures."""

from __future__ import annotations

import hashlib
import json
import shutil
from pathlib import Path

import pandas as pd
import pytest

from library.pipelines.tissue import TISSUE_COLUMN_ORDER
from library.pipelines.tissue.pipeline import (
    TissuePipelineOptions,
    run_tissue_pipeline,
)


class _SequentialChemblStub:
    """Return pre-recorded JSON payloads in the order they were requested."""

    def __init__(self, payloads: list[dict[str, object]]) -> None:
        self._payloads = payloads
        self.calls: list[str] = []

    def request_json(
        self, url: str, *, cfg, timeout: float | None
    ) -> dict[str, object]:
        self.calls.append(url)
        if not self._payloads:
            raise AssertionError("No recorded payload available for request")
        return self._payloads.pop(0)

    def close(self) -> None:  # pragma: no cover - interface compatibility
        return None


@pytest.mark.integration
def test_run_tissue_pipeline__replays_recorded_chembl_payloads(
    tmp_path: Path, cfg, snapshot_resource: Path
) -> None:
    """The pipeline processes recorded Chembl responses deterministically."""

    resource_dir = snapshot_resource / "tissue_pipeline"
    input_csv = resource_dir / "tissue_ids.csv"
    working_input = tmp_path / input_csv.name
    shutil.copyfile(input_csv, working_input)
    output_csv = tmp_path / "tissue_output.csv"

    def _load_payloads() -> list[dict[str, object]]:
        filenames = ["chembl_tissue_page1.json", "chembl_tissue_page2.json"]
        payloads: list[dict[str, object]] = []
        for name in filenames:
            with (resource_dir / name).open("r", encoding="utf-8") as handle:
                payloads.append(json.load(handle))
        return payloads

    recorded_payloads = _load_payloads()
    available_ids = {
        str(item["tissue_chembl_id"])
        for payload in recorded_payloads
        for item in payload.get("tissues", [])
        if isinstance(item, dict) and item.get("tissue_chembl_id")
    }

    requested_ids = pd.read_csv(working_input, dtype="string")[
        "tissue_chembl_id"
    ].tolist()
    expected_missing = tuple(
        identifier for identifier in requested_ids if identifier not in available_ids
    )

    options = TissuePipelineOptions(
        input_csv=working_input,
        output_csv=output_csv,
        column="tissue_chembl_id",
        batch_size=8,
        limit=None,
        offset=0,
        timeout=None,
    )

    stub = _SequentialChemblStub(_load_payloads())
    result = run_tissue_pipeline(cfg, options, client=stub)

    assert result.exit_code == 0
    output_df = pd.read_csv(output_csv, sep=cfg.io.csv_sep, dtype="string")
    assert result.records == len(output_df) == len(requested_ids)
    assert result.output_path == output_csv
    assert result.missing_ids == expected_missing
    assert stub.calls, "the stub should have been invoked at least once"
    assert not stub._payloads, "all recorded payloads should be consumed"

    expected_order = sorted(requested_ids)
    assert output_df["tissue_chembl_id"].tolist() == expected_order
    assert list(output_df.columns) == TISSUE_COLUMN_ORDER

    missing_id = expected_missing[0] if expected_missing else None
    if missing_id is not None:
        missing_row = output_df.loc[output_df["tissue_chembl_id"] == missing_id].iloc[0]
        for column in (
            col
            for col in TISSUE_COLUMN_ORDER
            if col not in {"tissue_chembl_id", "pipeline_version", "timestamp_utc"}
        ):
            assert pd.isna(missing_row[column])

    failure_path = output_csv.with_name(f"{output_csv.stem}_validation_failures.csv")
    assert not failure_path.exists()

    hash_before = hashlib.sha256(output_csv.read_bytes()).hexdigest()

    rerun_stub = _SequentialChemblStub(_load_payloads())
    rerun_result = run_tissue_pipeline(cfg, options, client=rerun_stub)
    assert rerun_result.exit_code == 0
    assert rerun_result.records == result.records
    assert rerun_result.missing_ids == expected_missing

    hash_after = hashlib.sha256(output_csv.read_bytes()).hexdigest()
    assert hash_before == hash_after

    rerun_df = pd.read_csv(output_csv, sep=cfg.io.csv_sep, dtype="string")
    pd.testing.assert_frame_equal(output_df, rerun_df, check_dtype=False)
    assert not failure_path.exists()
