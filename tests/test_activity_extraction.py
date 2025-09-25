"""Tests for :mod:`library.activity_extraction`."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from library.activity_extraction import extract_activities
from library.config import Config


class DummyChemblClient:
    """Test double replacing :class:`library.chembl_client.ChemblClient`."""

    def __init__(self, *args, **kwargs) -> None:  # noqa: D401 - signature matches context
        self.closed = False

    def __enter__(self) -> DummyChemblClient:
        return self

    def __exit__(self, exc_type, exc, tb) -> None:
        self.closed = True
        return None


@pytest.fixture()
def config(tmp_path: Path) -> Config:
    """Return a :class:`Config` instance writing outputs inside ``tmp_path``."""

    cfg = Config()
    cfg.io.output_dir = tmp_path
    return cfg


@pytest.fixture(autouse=True)
def stub_client(monkeypatch: pytest.MonkeyPatch) -> None:
    """Automatically replace :class:`ChemblClient` with a dummy implementation."""

    monkeypatch.setattr("library.activity_extraction.ChemblClient", DummyChemblClient)


@pytest.fixture(autouse=True)
def stub_quality(monkeypatch: pytest.MonkeyPatch) -> None:
    """Avoid generating table quality artefacts during tests."""

    def fake_quality(*args, **kwargs):  # type: ignore[no-untyped-def]
        return pd.DataFrame(), pd.DataFrame()

    monkeypatch.setattr(
        "library.activity_extraction.analyze_table_quality", fake_quality
    )


def _write_input(path: Path, identifiers: list[str]) -> None:
    rows = "\n".join(["activity_id"] + identifiers)
    path.write_text(rows, encoding="utf-8")


def test_extract_activities_success(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, config: Config
) -> None:
    """Successful run writes CSV and returns exit code ``0``."""

    input_csv = tmp_path / "activities.csv"
    output_csv = tmp_path / "result.csv"
    _write_input(input_csv, ["A1", "A2"])

    df = pd.DataFrame(
        {
            "activity_id": ["A1", "A2"],
            "molecule_chembl_id": ["M1", "M2"],
            "assay_chembl_id": ["AS1", "AS2"],
            "standard_value": [1.0, 2.0],
        }
    )

    def fake_get(ids, *, cfg, client, chunk_size, timeout):  # type: ignore[no-untyped-def]
        assert list(ids) == ["A1", "A2"]
        return df

    monkeypatch.setattr("library.activity_extraction.cl.get_activities", fake_get)

    exit_code = extract_activities(
        input_csv=input_csv,
        output_csv=output_csv,
        cfg=config,
        command="pytest",
    )

    assert exit_code == 0
    assert output_csv.exists()
    result = pd.read_csv(output_csv)
    expected_columns = set(df.columns) | {"pipeline_version", "timestamp_utc"}
    assert set(result.columns) == expected_columns
    assert result["activity_id"].tolist() == ["A1", "A2"]


def test_extract_activities_validation_failure(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, config: Config
) -> None:
    """Schema failures are exported to a sidecar CSV and signal an error."""

    input_csv = tmp_path / "activities.csv"
    output_csv = tmp_path / "result.csv"
    _write_input(input_csv, ["A1"])

    df = pd.DataFrame(
        {
            "activity_id": ["A1"],
            "molecule_chembl_id": ["M1"],
            "assay_chembl_id": ["AS1"],
            "standard_value": [1.0],
            "standard_type": ["invalid"],
        }
    )

    def fake_get(ids, *, cfg, client, chunk_size, timeout):  # type: ignore[no-untyped-def]
        return df

    monkeypatch.setattr("library.activity_extraction.cl.get_activities", fake_get)

    exit_code = extract_activities(
        input_csv=input_csv,
        output_csv=output_csv,
        cfg=config,
        command="pytest",
    )

    assert exit_code == 1
    failure_csv = output_csv.with_name("result_failure_cases.csv")
    assert failure_csv.exists()
    failure_df = pd.read_csv(failure_csv)
    assert not failure_df.empty
