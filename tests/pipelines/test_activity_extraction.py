"""Integration tests for :mod:`scripts.get_activity_data`."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest
import yaml

from library.config import Config


class DummyChemblClient:
    """Minimal stub replacing :class:`library.clients.ChemblClient`."""

    def __init__(
        self, *args, **kwargs
    ) -> None:  # noqa: D401 - signature mirrors context
        self.closed = False

    def __enter__(self) -> DummyChemblClient:
        return self

    def __exit__(self, exc_type, exc, tb) -> None:
        self.closed = True
        return None


@pytest.fixture()
def config_file(tmp_path: Path) -> Path:
    """Write a temporary configuration rooted in ``tmp_path``."""

    cfg = Config()
    cfg.io.output_dir = tmp_path / "output"
    cfg.io.cache_dir = tmp_path / "cache"
    cfg.activity.column = "activity_id"
    config_path = tmp_path / "config.yaml"
    config_path.write_text(
        yaml.safe_dump(cfg.model_dump(mode="json"), sort_keys=False),
        encoding="utf-8",
    )
    return config_path


@pytest.fixture(autouse=True)
def stub_quality(monkeypatch: pytest.MonkeyPatch) -> None:
    """Avoid generating table quality artefacts during tests."""

    def fake_quality(*args, **kwargs):  # type: ignore[no-untyped-def]
        return pd.DataFrame(), pd.DataFrame()

    monkeypatch.setattr(
        "scripts.get_activity_data.analyze_table_quality",
        fake_quality,
    )


def _write_input(path: Path, identifiers: list[str]) -> None:
    rows = "\n".join(["activity_id"] + identifiers)
    path.write_text(rows, encoding="utf-8")


def test_cli_success(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, config_file: Path
) -> None:
    """Successful CLI run writes CSV output and exits with code ``0``."""

    from scripts import get_activity_data

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

    monkeypatch.setattr("scripts.get_activity_data.ChemblClient", DummyChemblClient)

    def fake_get(ids, *, cfg, client, chunk_size, timeout, **kwargs):  # type: ignore[no-untyped-def]
        assert list(ids) == ["A1", "A2"]
        return df

    monkeypatch.setattr("scripts.get_activity_data.cl.get_activities", fake_get)

    exit_code = get_activity_data.main(
        [
            "--config",
            str(config_file),
            "--input",
            str(input_csv),
            "--output",
            str(output_csv),
        ]
    )

    assert exit_code == 0
    assert output_csv.exists()
    result = pd.read_csv(output_csv)
    expected_columns = set(df.columns) | {
        "pipeline_version",
        "timestamp_utc",
        "lower_value",
        "upper_value",
        "action_type",
        "activity_properties",
        "properties_hash",
    }
    assert set(result.columns) == expected_columns
    assert result["activity_id"].tolist() == ["A1", "A2"]


def test_cli_validation_failure(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, config_file: Path
) -> None:
    """Schema failures are exported to a sidecar CSV and signal an error."""

    from scripts import get_activity_data

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

    monkeypatch.setattr("scripts.get_activity_data.ChemblClient", DummyChemblClient)

    def fake_get(ids, *, cfg, client, chunk_size, timeout, **kwargs):  # type: ignore[no-untyped-def]
        return df

    monkeypatch.setattr("scripts.get_activity_data.cl.get_activities", fake_get)

    exit_code = get_activity_data.main(
        [
            "--config",
            str(config_file),
            "--input",
            str(input_csv),
            "--output",
            str(output_csv),
        ]
    )

    assert exit_code == 1
    failure_csv = output_csv.with_name("result_failure_cases.csv")
    assert failure_csv.exists()
    failure_df = pd.read_csv(failure_csv)
    assert not failure_df.empty
