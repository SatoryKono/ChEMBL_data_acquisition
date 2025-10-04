from __future__ import annotations

from pathlib import Path

import pandas as pd
import pandas.testing as pdt
import pytest

from library.postprocessing.target.main import postprocess_target_table


@pytest.mark.unit
def test_postprocess_target_table__produces_power_query_equivalent_csv(
    tmp_path: Path, snapshot_resource: Path
) -> None:
    input_path = snapshot_resource / "target_postprocess_input.csv"
    expected_path = snapshot_resource / "target_postprocess_expected.csv"

    working_input = tmp_path / input_path.name
    working_input.write_bytes(input_path.read_bytes())

    output_location = postprocess_target_table(working_input)
    output_path = Path(output_location)

    assert output_path == tmp_path / f"organism.{input_path.name}"
    assert output_path.read_bytes() == expected_path.read_bytes()

    result_frame = pd.read_csv(output_path, dtype=str, keep_default_na=False)
    expected_frame = pd.read_csv(expected_path, dtype=str, keep_default_na=False)
    pdt.assert_frame_equal(result_frame, expected_frame)

    bool_frame = pd.read_csv(
        output_path,
        dtype={"multifunctional_enzyme": "boolean"},
        keep_default_na=False,
    )
    assert str(bool_frame["multifunctional_enzyme"].dtype) == "boolean"

