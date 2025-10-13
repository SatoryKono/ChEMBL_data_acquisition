"""Postprocessing tests covering activity pipeline output pruning."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from library.cli.entrypoints import activity


def test_activity_final_output__drops_declared_empty_columns(tmp_path: Path) -> None:
    """Ensure empty activity columns are removed prior to CSV export."""

    frame = pd.DataFrame({"activity_chembl_id": ["ACT-1"], "value": [1.0]})
    for column in activity.EMPTY_COLUMNS:
        frame[column] = pd.NA

    pruned = activity._drop_empty_activity_columns(frame)

    output_path = tmp_path / "activity.csv"
    pruned.to_csv(output_path, index=False)

    reloaded = pd.read_csv(output_path)

    assert "activity_chembl_id" in reloaded.columns
    assert "value" in reloaded.columns
    for column in activity.EMPTY_COLUMNS:
        assert column not in reloaded.columns


def test_activity_final_output__keeps_non_empty_columns(tmp_path: Path) -> None:
    """Columns with real data are preserved when pruning empty columns."""

    populated_columns = {
        "compound_key": ["CMP-1"],
        "salt_chembl_id": ["CHEMBL123"],
        "nstereo": [0],
    }
    frame = pd.DataFrame({"activity_chembl_id": ["ACT-1"], **populated_columns})

    pruned = activity._drop_empty_activity_columns(frame)

    output_path = tmp_path / "activity.csv"
    pruned.to_csv(output_path, index=False)

    reloaded = pd.read_csv(output_path)

    for column in populated_columns:
        assert column in reloaded.columns
