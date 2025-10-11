from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest

from library.cli.pipeline_definition import PipelineDefinition
from library.cli_utils import RunPipelineResult, run_pipeline
from library.config import IoCfg


def _make_cfg(tmp_path: Path) -> SimpleNamespace:
    io_cfg = IoCfg(
        output_dir=tmp_path,
        exist_ok=True,
        csv_sep=",",
        csv_encoding="utf-8",
        default_date_prefix="20240101",
    )
    system_cfg = SimpleNamespace(doc_quality=SimpleNamespace(fatal_on_error=False))
    return SimpleNamespace(io=io_cfg, system=system_cfg)


@pytest.mark.integration
def test_run_pipeline__standard_outputs_sorted_by_multiple_keys(tmp_path: Path) -> None:
    frames = [
        pd.DataFrame(
            {
                "chembl_id": ["CHEMBL1", "CHEMBL2"],
                "batch": [2, 1],
                "value": [10, 20],
            }
        ),
        pd.DataFrame(
            {
                "chembl_id": ["CHEMBL1", "CHEMBL3"],
                "batch": [1, 0],
                "value": [5, 30],
            }
        ),
    ]

    def fetcher() -> list[pd.DataFrame]:
        return [chunk.copy() for chunk in frames]

    def writer(
        chunks: list[pd.DataFrame],
        destination: Path,
        col_order: list[str] | None,
        key_cols: list[str],
    ) -> Path:
        raise AssertionError(
            "legacy writer must not be invoked when standard outputs are enabled"
        )

    definition = PipelineDefinition(
        schema=None,
        schema_name="TestSchema",
        writer=writer,
        validators=(),
        metadata_hooks=(),
        command="test",
        config_snapshot={},
        inputs={},
        key_columns=("chembl_id", "batch"),
    )

    output_path = tmp_path / "output.documents_20240101.csv"
    failure_path = tmp_path / "failures.csv"
    cfg = _make_cfg(tmp_path)

    result = run_pipeline(
        definition=definition,
        fetcher=fetcher,
        output_path=output_path,
        failure_path=failure_path,
        cfg=cfg,
        emit_standard_outputs=True,
        emit_legacy_artifacts=False,
    )

    assert isinstance(result, RunPipelineResult)
    assert int(result) == 0
    artifacts = result.artifacts
    assert artifacts is not None, "standard outputs must be produced"
    assert artifacts.dataset.exists()

    dataset = pd.read_csv(artifacts.dataset)
    expected = (
        pd.concat(frames, ignore_index=True)
        .sort_values(by=["chembl_id", "batch"], kind="mergesort")
        .reset_index(drop=True)
    )
    expected = expected.reindex(columns=list(dataset.columns))

    pd.testing.assert_frame_equal(dataset.reset_index(drop=True), expected)

    assert list(dataset["batch"]) == [1, 2, 1, 0]
