from collections.abc import Iterable, Iterator
from datetime import date
from pathlib import Path
from typing import Any

import pandas as pd
import pytest
import yaml

from library.utils.cli_tools import csv_utils_main as cli


def _fd_count() -> int:
    fd_dir = Path("/proc/self/fd")
    if not fd_dir.exists():
        pytest.skip("/proc filesystem required to inspect file descriptors")
    return sum(1 for _ in fd_dir.iterdir())


class DummyReader(Iterator[pd.DataFrame]):
    def __init__(self, frames: Iterable[pd.DataFrame]):
        self._frames = iter(frames)
        self.closed = False

    def __iter__(self) -> "DummyReader":
        return self

    def __next__(self) -> pd.DataFrame:
        return next(self._frames)

    def close(self) -> None:
        self.closed = True

    def __enter__(self) -> "DummyReader":
        return self

    def __exit__(self, *exc_info: object) -> bool:
        self.close()
        return False


def test_cli_arguments_passed(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    input_csv = tmp_path / "in.csv"
    input_csv.write_text("a|b\n1|2\n", encoding="latin1")
    output_csv = tmp_path / "out.csv"
    called: dict[str, object] = {}

    def fake_read_csv(
        path: Path | str,
        sep: str,
        encoding: str,
        chunksize: int,
        *args: Any,
        **kwargs: Any,
    ) -> DummyReader:
        called["path"] = path
        called["sep"] = sep
        called["encoding"] = encoding
        called["chunksize"] = chunksize

        return DummyReader([pd.DataFrame({"a": [1], "b": [2]})])

    def fake_write(
        chunks: Iterable[pd.DataFrame],
        output: Path | str,
        col_order: list[str] | None = None,
        key_cols: list[str] | None = None,
        chunksize: int | None = None,
        merge_chunksize: int | None = None,
        drop_unexpected_cols: bool = False,
        *args: Any,
        **kwargs: Any,
    ) -> None:
        called["output"] = output
        called["write_chunksize"] = chunksize
        called["merge_chunksize"] = merge_chunksize
        called["write_sep"] = kwargs.get("sep")
        called["write_encoding"] = kwargs.get("encoding")
        list(chunks)  # exhaust generator

    monkeypatch.setattr(pd, "read_csv", fake_read_csv)
    monkeypatch.setattr(cli, "write_csv_chunks_deterministic", fake_write)

    rc = cli.main(
        [
            "--input",
            str(input_csv),
            "--output",
            str(output_csv),
            "--sep",
            "|",
            "--encoding",
            "latin1",
            "--key-cols",
            "a",
            "--chunk-size",
            "2",
            "--merge-chunk-size",
            "3",
            "--log-level",
            "INFO",
        ]
    )
    assert rc == 0
    assert called["path"] == input_csv
    assert called["output"] == output_csv
    assert called["sep"] == "|"
    assert called["encoding"] == "latin1"

    assert called["chunksize"] == 2
    assert called["write_chunksize"] == 2
    assert called["merge_chunksize"] == 3
    assert called["write_sep"] == "|"
    assert called["write_encoding"] == "latin1"


def test_cli_generates_output_path(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    input_csv = tmp_path / "input.csv"
    input_csv.write_text("a,b\n1,2\n", encoding="utf8")
    called: dict[str, Path | str] = {}

    def fake_read_csv(
        path: Path | str,
        sep: str,
        encoding: str,
        chunksize: int,
        *args: Any,
        **kwargs: Any,
    ) -> DummyReader:
        called["path"] = path
        return DummyReader([pd.DataFrame({"a": [1], "b": [2]})])

    def fake_write(
        chunks: Iterable[pd.DataFrame],
        output: Path | str,
        col_order: list[str] | None = None,
        key_cols: list[str] | None = None,
        chunksize: int | None = None,
        merge_chunksize: int | None = None,
        drop_unexpected_cols: bool = True,
        *args: Any,
        **kwargs: Any,
    ) -> None:
        called["output"] = output
        list(chunks)

    monkeypatch.setattr(pd, "read_csv", fake_read_csv)
    monkeypatch.setattr(cli, "write_csv_chunks_deterministic", fake_write)
    monkeypatch.setattr(
        cli,
        "date",
        type("D", (), {"today": staticmethod(lambda: date(2024, 1, 2))}),
    )

    rc = cli.main(["--input", str(input_csv), "--key-cols", "a"])
    assert rc == 0
    expected = input_csv.with_name("output.input_20240102.csv")
    assert called["output"] == expected


def test_cli_uses_configured_delimiter(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    input_csv = tmp_path / "in.csv"
    input_csv.write_text("a,b\n1,2\n", encoding="utf8")
    config_path = tmp_path / "cfg.yaml"
    config_data = yaml.safe_load(Path("config/config.yaml").read_text(encoding="utf8"))
    config_data["local"]["io"]["csv_sep"] = ";"
    config_data["local"]["io"]["csv_encoding"] = "latin1"
    config_data["local"]["io"]["csv_chunksize"] = 256
    config_path.write_text(yaml.safe_dump(config_data), encoding="utf8")

    called: dict[str, object] = {}

    def fake_read_csv(
        path: Path | str,
        sep: str,
        encoding: str,
        chunksize: int,
        *args: Any,
        **kwargs: Any,
    ) -> DummyReader:
        called["path"] = path
        called["sep"] = sep
        called["encoding"] = encoding
        called["chunksize"] = chunksize

        return DummyReader([pd.DataFrame({"a": [1], "b": [2]})])

    def fake_write(
        chunks: Iterable[pd.DataFrame],
        output: Path | str,
        col_order: list[str] | None = None,
        key_cols: list[str] | None = None,
        chunksize: int | None = None,
        merge_chunksize: int | None = None,
        drop_unexpected_cols: bool = False,
        *args: Any,
        **kwargs: Any,
    ) -> None:
        called["output"] = output
        called["write_chunksize"] = chunksize
        called["merge_chunksize"] = merge_chunksize
        called["write_sep"] = kwargs.get("sep")
        called["write_encoding"] = kwargs.get("encoding")
        list(chunks)

    monkeypatch.setattr(pd, "read_csv", fake_read_csv)
    monkeypatch.setattr(cli, "write_csv_chunks_deterministic", fake_write)

    rc = cli.main(
        [
            "--input",
            str(input_csv),
            "--key-cols",
            "a",
            "--config",
            str(config_path),
        ]
    )
    assert rc == 0
    assert called["path"] == input_csv
    assert called["sep"] == ";"
    assert called["encoding"] == "latin1"
    assert called["chunksize"] == 256
    assert called["write_sep"] == ";"
    assert called["write_encoding"] == "latin1"
    assert called["write_chunksize"] == 256


def test_cli_does_not_leak_file_descriptors(tmp_path: Path) -> None:
    input_csv = tmp_path / "input.csv"
    input_csv.write_text("id,value\n1,a\n2,b\n", encoding="utf8")

    baseline = _fd_count()
    base_args = [
        "--input",
        str(input_csv),
        "--key-cols",
        "id",
        "--chunk-size",
        "1",
    ]

    for run in range(3):
        output = tmp_path / f"out_{run}.csv"
        rc = cli.main([*base_args, "--output", str(output)])
        assert rc == 0
        assert output.exists()

    assert _fd_count() == baseline
