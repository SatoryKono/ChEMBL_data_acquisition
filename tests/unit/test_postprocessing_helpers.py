from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from library.postprocessing import helpers


@pytest.mark.unit
def test_read_csv_strict__applies_encoding_sep_and_na(tmp_path: Path) -> None:
    path = tmp_path / "sample.csv"
    path.write_text(
        "name;flag;count\n"
        "Alpha;TRUE;1\n"
        "Brav\xf3;[#N/A];\n"
        "Gamma;false;[#N/A]\n",
        encoding="cp1252",
    )

    frame = helpers.read_csv_strict(
        path,
        encoding=("utf-8", "cp1252"),
        dtype_map={"name": "Text", "flag": "Logical", "count": "Int64"},
        na_values=("[#N/A]",),
        separators=(";",),
    )

    assert frame["name"].dtype == "string"
    assert frame["name"].tolist() == ["Alpha", "Brav\xf3", "Gamma"]

    flags = frame["flag"]
    assert flags.dtype == "boolean"
    assert bool(flags.iloc[0]) is True
    assert flags.iloc[1] is pd.NA
    assert bool(flags.iloc[2]) is False

    counts = frame["count"]
    assert counts.dtype == "Int64"
    assert counts.iloc[0] == 1
    assert counts.iloc[1] is pd.NA
    assert counts.iloc[2] is pd.NA


@pytest.mark.unit
def test_coerce_types__invalid_logical_raises() -> None:
    frame = pd.DataFrame({"flag": ["maybe"]})

    with pytest.raises(ValueError):
        helpers.coerce_types(frame, {"flag": "Logical"})


@pytest.mark.unit
def test_stable_sort__preserves_duplicate_order() -> None:
    frame = pd.DataFrame(
        {
            "key": [2, 1, 1, 2],
            "value": ["delta", "alpha", "beta", "gamma"],
            "seq": [0, 1, 2, 3],
        }
    )

    sorted_frame = helpers.stable_sort(frame, ["key"])

    subset = sorted_frame[sorted_frame["key"] == 1]
    assert subset["seq"].tolist() == [1, 2]

    untouched = helpers.stable_sort(frame, [])
    assert untouched.equals(frame)


@pytest.mark.unit
def test_read_csv_with_fallbacks__uses_iso_8859_alias(
    monkeypatch, tmp_path: Path
) -> None:
    path = tmp_path / "document.csv"
    # 0x81 cannot be decoded using cp1252 but is valid in ISO-8859-1.
    raw = "document_chembl_id\tcompleted\nCHEMBL1\t\u0081\n"
    path.write_bytes(raw.encode("iso-8859-1"))

    original_read_csv = pd.read_csv

    def fake_read_csv(*args, **kwargs):
        encoding = kwargs.get("encoding")
        if encoding in {"utf-8", "utf-8-sig", "cp1252", "latin-1"}:
            raise UnicodeDecodeError(str(encoding), b"", 0, 1, "mock failure")
        return original_read_csv(*args, **kwargs)

    monkeypatch.setattr(pd, "read_csv", fake_read_csv)

    frame = helpers.read_csv_with_fallbacks(path, sep="\t")

    assert frame.columns.tolist() == ["document_chembl_id", "completed"]
    assert frame.iloc[0]["completed"] == "\u0081"
