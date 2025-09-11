from __future__ import annotations

from pathlib import Path

import pandas as pd
from pandera.errors import SchemaErrors

from library.sidecar import SidecarErrors
from schemas import ActivitiesSchema


def _validate(df: pd.DataFrame, out_csv: Path, failure_csv: Path) -> pd.DataFrame:
    """Validate *df* against :data:`ActivitiesSchema` and collect errors."""
    sc = SidecarErrors()
    try:
        validated = ActivitiesSchema.validate(df, lazy=True)
    except SchemaErrors as exc:
        for row in exc.failure_cases.to_dict(orient="records"):
            sc.add_error(row)
        validated = df.drop(index=exc.failure_cases["index"].unique())
    sc.save(failure_csv)
    validated.to_csv(out_csv, index=False)
    return validated


def test_activity_validation_pass(tmp_path: Path) -> None:
    df = pd.read_csv(Path("tests/data/activities_valid.csv"))
    out_csv = tmp_path / "out.csv"
    failure_csv = tmp_path / "failure_cases.csv"
    result = _validate(df, out_csv, failure_csv)
    assert result.equals(df)
    assert not failure_csv.exists()


def test_activity_validation_failures(tmp_path: Path) -> None:
    df = pd.read_csv(Path("tests/data/activities_mixed.csv"))
    out_csv = tmp_path / "out.csv"
    failure_csv = tmp_path / "failure_cases.csv"
    result = _validate(df, out_csv, failure_csv)
    assert len(result) == 1
    assert (result["standard_value"] >= 0).all()
    assert failure_csv.exists()
    sidecar_df = pd.read_csv(failure_csv)
    assert (sidecar_df["column"] == "standard_value").any()
    assert (sidecar_df["index"] == 1).any()
