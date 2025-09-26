from collections.abc import Iterable, Iterator
from datetime import date
from pathlib import Path
from typing import Any

import pandas as pd
import pytest

from library.utils.cli_tools import csv_utils_main as cli


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
    ) -> Iterator[pd.DataFrame]:
        called["path"] = path
        called["sep"] = sep
        called["encoding"] = encoding
        called["chunksize"] = chunksize

        def gen() -> Iterator[pd.DataFrame]:
            yield pd.DataFrame({"a": [1], "b": [2]})

        return gen()

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
    ) -> Iterator[pd.DataFrame]:
        called["path"] = path

        def gen() -> Iterator[pd.DataFrame]:
            yield pd.DataFrame({"a": [1], "b": [2]})

        return gen()

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
    expected = input_csv.with_name("output_input_20240102.csv")
    assert called["output"] == expected
