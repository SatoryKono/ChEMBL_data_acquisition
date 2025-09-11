from datetime import date
from pathlib import Path

import pandas as pd

import csv_utils_main as cli


def test_cli_arguments_passed(monkeypatch, tmp_path: Path) -> None:
    input_csv = tmp_path / "in.csv"
    input_csv.write_text("a|b\n1|2\n", encoding="latin1")
    output_csv = tmp_path / "out.csv"
    called: dict[str, object] = {}

    def fake_read_csv(path, sep, encoding, chunksize):  # type: ignore[override]
        called["path"] = path
        called["sep"] = sep
        called["encoding"] = encoding
        called["chunksize"] = chunksize


        def gen():
            yield pd.DataFrame({"a": [1], "b": [2]})

        return gen()

    def fake_write(
        chunks,
        output,
        col_order=None,
        key_cols=None,
        chunksize=None,
        drop_unexpected_cols=False,

    ):  # type: ignore[override]
        called["output"] = output
        called["write_chunksize"] = chunksize
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
            "--chunk-size",
            "2",
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



def test_cli_generates_output_path(monkeypatch, tmp_path: Path) -> None:
    input_csv = tmp_path / "input.csv"
    input_csv.write_text("a,b\n1,2\n", encoding="utf8")
    called: dict[str, Path] = {}

    def fake_read_csv(path, sep, encoding):  # type: ignore[override]
        called["path"] = path
        return pd.DataFrame({"a": [1], "b": [2]})

    def fake_write(
        df,
        output,
        col_order=None,
        key_cols=None,
        drop_unexpected_cols=True,
    ):  # type: ignore[override]
        called["output"] = output

    monkeypatch.setattr(pd, "read_csv", fake_read_csv)
    monkeypatch.setattr(cli, "write_csv_deterministic", fake_write)
    monkeypatch.setattr(
        cli,
        "date",
        type("D", (), {"today": staticmethod(lambda: date(2024, 1, 2))}),
    )

    rc = cli.main(["--input", str(input_csv)])
    assert rc == 0
    expected = input_csv.with_name("output_input_20240102.csv")
    assert called["output"] == expected

