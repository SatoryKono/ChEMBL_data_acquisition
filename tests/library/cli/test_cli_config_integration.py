from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

from library.utils.config import DEFAULT_CONFIG_PATH
from library.utils.cli_tools import mapper_batch_main
from library.utils.cli_tools import chunk_io_main
from library.utils.cli_tools import csv_utils_main
from library.utils.cli_tools.table_quality_main import main


def test_table_quality_cli_with_config(tmp_path: Path) -> None:
    csv_path = tmp_path / "data.csv"
    csv_path.write_text("a\n1\n")
    config = Path("tests/fixtures/config.min.yaml")
    rc = main(
        [
            "--config",
            str(config),
            "--input",
            str(csv_path),
            "--table-name",
            "demo",
            "--output",
            str(tmp_path),
        ]
    )
    assert rc == 0


def test_mapper_batch_default_config_outside_repo(tmp_path: Path, monkeypatch) -> None:
    recorded: dict[str, Path] = {}

    def fake_apply_config_overrides(
        args,
        parser,
        config_path,
        mapping=None,
        **kwargs,
    ):
        recorded["config_path"] = Path(config_path)
        return SimpleNamespace(io=SimpleNamespace(exist_ok=True))

    monkeypatch.setattr(
        mapper_batch_main.cli,
        "apply_config_overrides",
        fake_apply_config_overrides,
    )
    monkeypatch.setattr(mapper_batch_main, "ensure_dirs", lambda cfg: None)
    monkeypatch.setattr(
        mapper_batch_main,
        "run",
        lambda cfg, args: 0,
    )
    monkeypatch.chdir(tmp_path)

    rc = mapper_batch_main.main([])

    assert rc == 0
    assert recorded["config_path"] == DEFAULT_CONFIG_PATH


def test_chunk_io_default_config_outside_repo(tmp_path: Path, monkeypatch) -> None:
    recorded: dict[str, Path] = {}

    def fake_apply_config_overrides(
        args,
        parser,
        config_path,
        mapping=None,
        **kwargs,
    ):
        recorded["config_path"] = Path(config_path)
        return SimpleNamespace(io=SimpleNamespace(exist_ok=True))

    monkeypatch.setattr(
        chunk_io_main.cli,
        "apply_config_overrides",
        fake_apply_config_overrides,
    )
    monkeypatch.setattr(chunk_io_main, "ensure_dirs", lambda cfg: None)
    monkeypatch.setattr(chunk_io_main, "process_csv_chunks", lambda *args, **kwargs: 0)

    input_csv = tmp_path / "input.csv"
    input_csv.write_text("id\n1\n")
    output_csv = tmp_path / "output.csv"

    monkeypatch.chdir(tmp_path)

    rc = chunk_io_main.main([
        "--input",
        str(input_csv),
        "--output",
        str(output_csv),
        "--log-level",
        "ERROR",
    ])

    assert rc == 0
    assert recorded["config_path"] == DEFAULT_CONFIG_PATH


def test_csv_utils_default_config_outside_repo(tmp_path: Path, monkeypatch) -> None:
    recorded: dict[str, Path] = {}

    io_cfg = SimpleNamespace(
        exist_ok=True,
        csv_sep=",",
        csv_encoding="utf8",
        csv_chunksize=10,
        csv_dtype=str,
    )

    def fake_apply_config_overrides(
        args,
        parser,
        config_path,
        mapping=None,
        **kwargs,
    ):
        recorded["config_path"] = Path(config_path)
        return SimpleNamespace(io=io_cfg)

    monkeypatch.setattr(
        csv_utils_main.cli,
        "apply_config_overrides",
        fake_apply_config_overrides,
    )
    monkeypatch.setattr(csv_utils_main, "ensure_dirs", lambda cfg: None)

    def fake_writer(reader, output, **kwargs):
        list(reader)
        output.write_text("done\n")

    monkeypatch.setattr(
        csv_utils_main,
        "write_csv_chunks_deterministic",
        fake_writer,
    )

    input_csv = tmp_path / "input.csv"
    input_csv.write_text("id,value\n1,2\n")
    output_csv = tmp_path / "output.csv"

    monkeypatch.chdir(tmp_path)

    rc = csv_utils_main.main([
        "--input",
        str(input_csv),
        "--output",
        str(output_csv),
        "--key-cols",
        "id",
        "--log-level",
        "ERROR",
    ])

    assert rc == 0
    assert recorded["config_path"] == DEFAULT_CONFIG_PATH
