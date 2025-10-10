"""End-to-end checks for :mod:`library.utils.cli_tools.chunk_io_main`."""

from __future__ import annotations

import csv
from pathlib import Path
from types import SimpleNamespace

import pytest

from library.utils.cli_tools import chunk_io_main


@pytest.mark.e2e
def test_chunk_io_main__date_overrides_config_prefix(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    input_csv = tmp_path / "input.csv"
    with input_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["id", "value"])
        writer.writerow(["row-1", "42"])

    output_dir = tmp_path / "exports"
    checkpoint = tmp_path / "state.json"
    cli_date = "20240707"
    cfg_date = "19990101"

    cfg_stub = SimpleNamespace(
        io=SimpleNamespace(
            exist_ok=True,
            output_dir=output_dir,
            cache_dir=tmp_path / "cache",
            csv_sep=",",
            csv_encoding="utf-8",
            default_date_prefix=cfg_date,
        )
    )

    def _fake_apply_config(args, parser, config_path):
        return cfg_stub

    original_prepare = chunk_io_main.cli.prepare_io_paths

    def _prepare_and_clear(args, *positional, **keyword):
        original_prepare(args, *positional, **keyword)
        args.output_csv = None
        args.final_out = None

    observed: dict[str, Path] = {}

    def _fake_process(input_path, output_path, **_kwargs):
        output = Path(output_path)
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text("id,value\nrow-1,42\n", encoding="utf-8")
        observed["output"] = output
        return 1

    monkeypatch.setattr(chunk_io_main.cli, "apply_config_overrides", _fake_apply_config)
    monkeypatch.setattr(chunk_io_main.cli, "prepare_io_paths", _prepare_and_clear)
    monkeypatch.setattr(chunk_io_main, "process_csv_chunks", _fake_process)
    monkeypatch.setattr(chunk_io_main, "ensure_dirs", lambda _cfg: None)

    exit_code = chunk_io_main.main(
        [
            "--input",
            str(input_csv),
            "--chunk-size",
            "1",
            "--checkpoint",
            str(checkpoint),
            "--date",
            cli_date,
        ]
    )

    expected_output = output_dir / f"output.{input_csv.stem}_{cli_date}.csv"

    assert exit_code == 0
    assert observed["output"] == expected_output
    assert expected_output.exists()


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
