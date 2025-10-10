from __future__ import annotations

from pathlib import Path

import pandas as pd
import pandas.testing as pdt
import pytest

from library.config import Config
from library.pipelines.target import helpers

INPUT_FILE = "target_postprocess_power_query_input.csv"
EXPECTED_FILE = "target_postprocess_power_query_expected.csv"


@pytest.mark.integration
def test_postprocess_target_table__matches_golden(
    tmp_path: Path, cfg: Config, snapshot_resource: Path
) -> None:
    input_path = snapshot_resource / INPUT_FILE
    expected_path = snapshot_resource / EXPECTED_FILE

    working_input = tmp_path / INPUT_FILE
    working_output = tmp_path / "targets_final.csv"
    working_input.write_bytes(input_path.read_bytes())

    helpers.postprocess_target_table_file(working_input, working_output, cfg=cfg)

    # Normalise EOL across platforms
    assert working_output.read_bytes().replace(
        b"\r\n", b"\n"
    ) == expected_path.read_bytes().replace(b"\r\n", b"\n")

    result = pd.read_csv(working_output, dtype=str, keep_default_na=False).astype(
        "string"
    )
    expected = pd.read_csv(expected_path, dtype=str, keep_default_na=False).astype(
        "string"
    )
    pdt.assert_frame_equal(result, expected)


@pytest.mark.integration
def test_postprocess_target_table__is_idempotent(
    tmp_path: Path, cfg: Config, snapshot_resource: Path
) -> None:
    input_path = snapshot_resource / INPUT_FILE
    expected_path = snapshot_resource / EXPECTED_FILE

    working_input = tmp_path / INPUT_FILE
    working_output = tmp_path / "targets_final.csv"
    working_input.write_bytes(input_path.read_bytes())

    helpers.postprocess_target_table_file(working_input, working_output, cfg=cfg)
    first_run = working_output.read_bytes()

    second_output = tmp_path / "targets_final_second.csv"
    helpers.postprocess_target_table_file(working_input, second_output, cfg=cfg)
    second_run = second_output.read_bytes()

    assert (
        first_run.replace(b"\r\n", b"\n")
        == second_run.replace(b"\r\n", b"\n")
        == expected_path.read_bytes().replace(b"\r\n", b"\n")
    )
