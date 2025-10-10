"""End-to-end checks for :mod:`library.utils.cli_tools.chunk_io_main`."""

from __future__ import annotations

import csv
from pathlib import Path
from types import SimpleNamespace

import pytest

from library.utils.cli_tools import chunk_io_main


@pytest.mark.e2e
def test_chunk_io_main__copies_csv(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    input_csv = tmp_path / "input.csv"
    rows = [
        ["molecule_chembl_id", "name", "smiles"],
        ["CHEMBL1", "Aspirin", "CC(=O)OC1=CC=CC=C1C(=O)O"],
        ["CHEMBL2", "Unknown", ""],
    ]
    with input_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerows(rows)

    output_csv = tmp_path / "output" / "copy.csv"
    checkpoint = tmp_path / "state.json"

    monkeypatch.setenv("CHEMBL_DA_BASE_PATH", str(tmp_path / "chembl-data"))

    cfg_stub = SimpleNamespace(
        io=SimpleNamespace(
            exist_ok=True,
            output_dir=tmp_path,
            cache_dir=tmp_path / "cache",
            csv_sep=",",
            csv_encoding="utf-8",
        )
    )

    def _cfg_to_dict() -> dict[str, object]:
        return {
            "io": {
                "exist_ok": True,
                "output_dir": str(tmp_path),
                "cache_dir": str(tmp_path / "cache"),
                "csv_sep": ",",
                "csv_encoding": "utf-8",
            }
        }

    cfg_stub.to_dict = _cfg_to_dict  # type: ignore[attr-defined]

    captured: dict[str, object] = {}

    def _fake_apply_config(args, parser, config_path):
        captured["config_path"] = config_path
        return cfg_stub

    monkeypatch.setattr(chunk_io_main.cli, "apply_config_overrides", _fake_apply_config)

    exit_code = chunk_io_main.main(
        [
            "--input",
            str(input_csv),
            "--final-out",
            str(output_csv),
            "--checkpoint",
            str(checkpoint),
            "--chunk-size",
            "2",
        ]
    )

    assert exit_code == 0
    assert captured["config_path"] == chunk_io_main.DEFAULT_CONFIG_PATH
    assert output_csv.exists()
    assert checkpoint.exists()

    with output_csv.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.reader(handle)
        written_rows = list(reader)

    assert written_rows == rows


@pytest.mark.e2e
def test_chunk_io_main__resolves_checkpoint_with_base_path(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    input_csv = tmp_path / "input.csv"
    rows = [
        ["molecule_chembl_id", "name", "smiles"],
        ["CHEMBL1", "Aspirin", "CC(=O)OC1=CC=CC=C1C(=O)O"],
        ["CHEMBL2", "Unknown", ""],
    ]
    with input_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerows(rows)

    base_path = tmp_path / "chembl-data"
    output_dir_name = "outputs"
    final_output_name = "copy.csv"
    checkpoint_name = "state.json"

    cfg_stub = SimpleNamespace(
        io=SimpleNamespace(
            exist_ok=True,
            output_dir=base_path / output_dir_name,
            cache_dir=base_path / "cache",
            csv_sep=",",
            csv_encoding="utf-8",
        )
    )

    def _cfg_to_dict() -> dict[str, object]:
        return {
            "io": {
                "exist_ok": True,
                "output_dir": str(base_path / output_dir_name),
                "cache_dir": str(base_path / "cache"),
                "csv_sep": ",",
                "csv_encoding": "utf-8",
            }
        }

    cfg_stub.to_dict = _cfg_to_dict  # type: ignore[attr-defined]

    captured: dict[str, object] = {}

    def _fake_apply_config(args, parser, config_path):
        captured["config_path"] = config_path
        return cfg_stub

    monkeypatch.setattr(chunk_io_main.cli, "apply_config_overrides", _fake_apply_config)

    exit_code = chunk_io_main.main(
        [
            "--input",
            str(input_csv),
            "--final-out",
            final_output_name,
            "--checkpoint",
            checkpoint_name,
            "--chunk-size",
            "2",
            "--base-path",
            str(base_path),
            "--output-dir",
            output_dir_name,
        ]
    )

    expected_output = base_path / output_dir_name / final_output_name
    expected_checkpoint = base_path / output_dir_name / checkpoint_name

    assert exit_code == 0
    assert captured["config_path"] == chunk_io_main.DEFAULT_CONFIG_PATH
    assert expected_output.exists()
    assert expected_checkpoint.exists()

    with expected_output.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.reader(handle)
        written_rows = list(reader)

    assert written_rows == rows
