import os
import random
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from library.cli.pipeline_definition import PipelineDefinition
from library.cli_utils import run_pipeline
from library.config import IoCfg


def _fix_seed(seed: int = 42) -> None:
    os.environ["PYTHONHASHSEED"] = str(seed)
    random.seed(seed)
    np.random.seed(seed)


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


@pytest.mark.slow
def test_run_pipeline_streaming__large_dataset(tmp_path: Path) -> None:
    _fix_seed()
    total_rows = 50_000
    chunk_size = 2_048

    def fetcher():
        for start in range(0, total_rows, chunk_size):
            stop = min(start + chunk_size, total_rows)
            yield pd.DataFrame(
                {
                    "identifier": [f"id-{idx}" for idx in range(start, stop)],
                    "value": list(range(start, stop)),
                }
            )

    def writer(
        chunks: list[pd.DataFrame],
        destination: Path,
        col_order: list[str] | None,
        key_cols: list[str],
    ) -> Path:
        raise AssertionError("legacy writer must not be used in streaming mode")

    definition = PipelineDefinition(
        schema=None,
        schema_name="PerfSchema",
        writer=writer,
        validators=(),
        metadata_hooks=(),
        command="perf-test",
        config_snapshot={},
        inputs={},
        key_columns=("identifier",),
    )

    output_path = tmp_path / "output.perf_20240101.csv"
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

    assert int(result) == 0
    artifacts = result.artifacts
    assert artifacts is not None
    dataset = pd.read_csv(artifacts.dataset)
    assert len(dataset) == total_rows
    values = dataset["identifier"].tolist()
    assert values[0] == "id-0"
    assert sorted(values) == values
    assert f"id-{total_rows - 1}" in values
