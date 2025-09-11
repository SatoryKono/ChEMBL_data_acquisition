from pathlib import Path

import pandas as pd

from scripts import csv_utils_main as cli


def test_cli_arguments_passed(monkeypatch, tmp_path: Path) -> None:
    input_csv = tmp_path / "in.csv"
    input_csv.write_text("a|b\n1|2\n", encoding="latin1")
    output_csv = tmp_path / "out.csv"
    called: dict[str, object] = {}

    def fake_read_csv(path, sep, encoding):  # type: ignore[override]
        called["path"] = path
        called["sep"] = sep
        called["encoding"] = encoding
        return pd.DataFrame({"a": [1], "b": [2]})

    def fake_write(df, output, col_order=None, key_cols=None):  # type: ignore[override]
        called["output"] = output

    monkeypatch.setattr(pd, "read_csv", fake_read_csv)
    monkeypatch.setattr(cli, "write_csv_deterministic", fake_write)

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
            "--log-level",
            "INFO",
        ]
    )
    assert rc == 0
    assert called["path"] == input_csv
    assert called["output"] == output_csv
    assert called["sep"] == "|"
    assert called["encoding"] == "latin1"
