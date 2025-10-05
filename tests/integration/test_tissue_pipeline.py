from __future__ import annotations

import json
from copy import deepcopy
from datetime import datetime, timezone
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest
import yaml

from library.config import Config
from library.pipelines.common import metadata as pipeline_metadata
from library.pipelines.tissue.pipeline import (
    TissuePipelineOptions,
    run_tissue_pipeline,
)
from library.pipelines.tissue.chembl import TISSUE_COLUMN_ORDER


class _StubChemblClient:
    """Deterministic :class:`ChemblClient` replacement for tests."""

    def __init__(self, responses: dict[str, dict[str, object]]) -> None:
        self._responses = {url: deepcopy(payload) for url, payload in responses.items()}
        self.calls: list[dict[str, object]] = []

    def request_json(self, url: str, *, cfg: object, timeout: float | None = None) -> dict:
        self.calls.append({"url": url, "timeout": timeout})
        if url not in self._responses:
            raise AssertionError(f"Unexpected request for {url}")
        return json.loads(json.dumps(self._responses[url]))


def _load_responses(base: str, resource_dir: Path) -> dict[str, dict[str, object]]:
    mappings: dict[str, dict[str, object]] = {}
    chunk0_page0 = json.loads(
        (resource_dir / "chembl_response_chunk0_page0.json").read_text(encoding="utf-8")
    )
    chunk0_page1 = json.loads(
        (resource_dir / "chembl_response_chunk0_page1.json").read_text(encoding="utf-8")
    )
    chunk1_page0 = json.loads(
        (resource_dir / "chembl_response_chunk1_page0.json").read_text(encoding="utf-8")
    )
    mappings[
        f"{base}/tissue.json?format=json&tissue_chembl_id__in=CHEMBL613507,CHEMBL2109249&limit=2"
    ] = chunk0_page0
    mappings[
        (
            f"{base}/tissue.json?format=json"
            "&tissue_chembl_id__in=CHEMBL613507,CHEMBL2109249&limit=2&offset=1"
        )
    ] = chunk0_page1
    mappings[f"{base}/tissue.json?format=json&tissue_chembl_id__in=CHEMBL999999&limit=1"] = (
        chunk1_page0
    )
    return mappings


@pytest.mark.integration
def test_run_tissue_pipeline__writes_sorted_deterministic_csv(
    tmp_path: Path, snapshot_resource: Path, cfg: Config, monkeypatch: pytest.MonkeyPatch
) -> None:
    resource_dir = snapshot_resource / "tissue"
    input_csv = tmp_path / "input" / "tissue_ids.csv"
    input_csv.parent.mkdir(parents=True, exist_ok=True)
    input_csv.write_text(
        (resource_dir / "input_ids.csv").read_text(encoding="utf-8"),
        encoding="utf-8",
    )
    output_csv = tmp_path / "output" / "tissue_export.csv"
    output_csv.parent.mkdir(parents=True, exist_ok=True)

    cfg.io.output_dir = output_csv.parent
    cfg.io.cache_dir = tmp_path / "cache"
    cfg.io.cache_dir.mkdir()
    cfg.io.exist_ok = True
    cfg.io.csv_encoding = "utf-8"
    cfg.tissue.batch_size = 2

    pipeline_metadata.get_timestamp_utc.cache_clear()
    pipeline_metadata.pipeline_metadata.cache_clear()
    monkeypatch.setattr(
        pipeline_metadata,
        "get_timestamp_utc",
        lambda: "2020-01-01T00:00:00+00:00",
    )
    pipeline_metadata.pipeline_metadata.cache_clear()

    frozen_now = datetime(2020, 1, 1, tzinfo=timezone.utc)

    def _frozen_datetime_now(tz: timezone | None = None) -> datetime:
        if tz is None:
            return frozen_now.replace(tzinfo=None)
        return frozen_now.astimezone(tz)

    monkeypatch.setattr(
        "library.io.metadata.datetime",
        SimpleNamespace(now=_frozen_datetime_now),
    )

    base = cfg.api.chembl_base.rstrip("/")
    responses = _load_responses(base, resource_dir)
    client = _StubChemblClient(responses)

    options = TissuePipelineOptions(
        input_csv=input_csv,
        output_csv=output_csv,
        column="tissue_chembl_id",
        batch_size=2,
        limit=None,
        offset=0,
        timeout=15.0,
    )

    result = run_tissue_pipeline(cfg, options, client=client)

    assert result.exit_code == 0
    assert result.records == 3
    assert result.written is True
    assert result.output_path == output_csv
    assert result.missing_ids == ("CHEMBL999999",)
    assert result.failure_path is None

    output_df = pd.read_csv(output_csv, dtype="string")
    expected_df = pd.read_csv(resource_dir / "expected_output.csv", dtype="string")
    assert list(output_df.columns) == TISSUE_COLUMN_ORDER
    assert list(output_df["tissue_chembl_id"]) == sorted(output_df["tissue_chembl_id"].tolist())
    pd.testing.assert_frame_equal(output_df, expected_df, check_dtype=False)

    meta_path = output_csv.with_suffix(".csv.meta.yaml")
    assert meta_path.exists()
    meta = yaml.safe_load(meta_path.read_text(encoding="utf-8"))
    assert meta["columns"] == TISSUE_COLUMN_ORDER
    assert meta["dtypes"]["tissue_chembl_id"].lower().startswith("string")

    expected_calls = list(responses.keys())
    observed_calls = [entry["url"] for entry in client.calls]
    assert observed_calls == expected_calls
