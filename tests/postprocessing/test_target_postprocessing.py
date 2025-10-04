from __future__ import annotations

import io
from pathlib import Path

import pandas as pd
import pytest

from library.postprocessing import target


@pytest.mark.unit
def test_normalise_series__strips_and_normalises():
    series = pd.Series(["  Alpha  ", "n/a", None, " NA "])

    result = target._normalise_series(series, lowercase=False)

    assert result.tolist() == ["Alpha", "", "", ""]


@pytest.mark.unit
def test_normalise_series__applies_lowercase():
    series = pd.Series(["Beta", " Mixed Case ", "NONE"])

    result = target._normalise_series(series, lowercase=True)

    assert result.tolist() == ["beta", "mixed case", ""]


@pytest.mark.unit
def test_normalise_series__empty_series_returns_copy():
    series = pd.Series(dtype=str)

    result = target._normalise_series(series, lowercase=True)

    assert result.empty
    assert result.dtype == series.dtype


@pytest.mark.unit
def test_split_pipe__filters_whitespace_and_blanks():
    assert target._split_pipe(" value | |second|") == ["value", "second"]
    assert target._split_pipe("") == []
    assert target._split_pipe(None) == []  # type: ignore[arg-type]


@pytest.mark.unit
def test_split_synonyms__deduplicates_and_skips_placeholders():
    value = "Alpha : beta : NA : Beta : none"

    assert target._split_synonyms(value) == ["Alpha", "beta", "Beta"]
    assert target._split_synonyms(None) == []  # type: ignore[arg-type]


@pytest.mark.unit
def test_syn_expand__tokenises_and_deduplicates():
    term = "Alpha-beta gamma Alpha"

    assert target._syn_expand(term) == ["alpha", "beta", "gamma"]


@pytest.mark.unit
def test_read_csv__tries_multiple_encodings(monkeypatch: pytest.MonkeyPatch, tmp_path: Path):
    attempted: list[str] = []

    def fake_read_csv(path: Path, *, dtype, encoding: str, sep: str):  # type: ignore[override]
        attempted.append(encoding)
        if encoding != "utf-8":
            raise UnicodeDecodeError("utf-8", b"", 0, 1, "boom")
        return pd.DataFrame({"column": ["value"]})

    monkeypatch.setattr(pd, "read_csv", fake_read_csv)

    csv_path = tmp_path / "input.csv"
    csv_path.write_text("column\nvalue\n", encoding="utf-8")

    frame = target._read_csv(csv_path, encodings=("utf-8-sig", "utf-8"), sep=",")

    assert attempted == ["utf-8-sig", "utf-8"]
    assert frame.to_dict(orient="list") == {"column": ["value"]}


@pytest.mark.unit
def test_read_csv__raises_last_unicode_error(monkeypatch: pytest.MonkeyPatch, tmp_path: Path):
    def failing_read_csv(path: Path, *, dtype, encoding: str, sep: str):  # type: ignore[override]
        raise UnicodeDecodeError(encoding, b"", 0, 1, "boom")

    monkeypatch.setattr(pd, "read_csv", failing_read_csv)

    csv_path = tmp_path / "input.csv"
    csv_path.write_text("column\nvalue\n", encoding="utf-8")

    with pytest.raises(UnicodeDecodeError):
        target._read_csv(csv_path, encodings=("utf-8",), sep=",")


@pytest.mark.unit
def test_process_targets__produces_expected_csv_structure(tmp_path: Path):
    input_path = tmp_path / "targets.csv"
    frame = pd.DataFrame(
        [
            {
                "target_chembl_id": "CHEMBL_T1",
                "isoform_ids": "ENSP0001",
                "isoform_names": "Alpha",
                "isoform_synonyms": "Alpha:Alpha Variant",
            }
        ]
    )
    frame.to_csv(input_path, index=False, encoding="utf-8")

    output_path = target.process_targets(input_path)
    result = pd.read_csv(output_path)

    assert list(result.columns) == [
        "target_chembl_id",
        "isoform_id",
        "isoform_name",
        "term",
        "token",
    ]
    assert result.to_dict(orient="records") == [
        {
            "target_chembl_id": "CHEMBL_T1",
            "isoform_id": "ENSP0001",
            "isoform_name": "Alpha",
            "term": "Alpha",
            "token": "alpha",
        },
        {
            "target_chembl_id": "CHEMBL_T1",
            "isoform_id": "ENSP0001",
            "isoform_name": "Alpha",
            "term": "alpha",
            "token": "alpha",
        },
        {
            "target_chembl_id": "CHEMBL_T1",
            "isoform_id": "ENSP0001",
            "isoform_name": "Alpha",
            "term": "alpha variant",
            "token": "alpha",
        },
        {
            "target_chembl_id": "CHEMBL_T1",
            "isoform_id": "ENSP0001",
            "isoform_name": "Alpha",
            "term": "alpha variant",
            "token": "variant",
        },
    ]


@pytest.mark.unit
def test_process_targets__filters_rows_without_terms(tmp_path: Path):
    input_path = tmp_path / "targets.csv"
    frame = pd.DataFrame(
        [
            {
                "target_chembl_id": "CHEMBL_T2",
                "isoform_ids": "",
                "isoform_names": "",
                "isoform_synonyms": "",
            }
        ]
    )
    frame.to_csv(input_path, index=False, encoding="utf-8")

    output_path = target.process_targets(input_path)
    result = pd.read_csv(output_path)

    assert result.empty
    assert list(result.columns) == [
        "target_chembl_id",
        "isoform_id",
        "isoform_name",
        "term",
        "token",
    ]


@pytest.mark.unit
def test_process_targets__drops_duplicate_rows(tmp_path: Path):
    input_path = tmp_path / "targets.csv"
    frame = pd.DataFrame(
        [
            {
                "target_chembl_id": "CHEMBL_T3",
                "isoform_ids": "ENSP0003",
                "isoform_names": "Gamma",
                "isoform_synonyms": "Gamma",
            },
            {
                "target_chembl_id": "CHEMBL_T3",
                "isoform_ids": "ENSP0003",
                "isoform_names": "Gamma",
                "isoform_synonyms": "Gamma",
            },
        ]
    )
    frame.to_csv(input_path, index=False, encoding="utf-8")

    output_path = target.process_targets(input_path)
    result = pd.read_csv(output_path)

    assert result.to_dict(orient="records") == [
        {
            "target_chembl_id": "CHEMBL_T3",
            "isoform_id": "ENSP0003",
            "isoform_name": "Gamma",
            "term": "Gamma",
            "token": "gamma",
        },
        {
            "target_chembl_id": "CHEMBL_T3",
            "isoform_id": "ENSP0003",
            "isoform_name": "Gamma",
            "term": "gamma",
            "token": "gamma",
        },
    ]


@pytest.mark.unit
def test_process_targets__sorts_output_rows(tmp_path: Path):
    input_path = tmp_path / "targets.csv"
    frame = pd.DataFrame(
        [
            {
                "target_chembl_id": "CHEMBL_T4",
                "isoform_ids": "ENSP0005|ENSP0004",
                "isoform_names": "Beta|Alpha",
                "isoform_synonyms": "Beta Alt:Beta|Alpha Alt:Alpha",
            }
        ]
    )
    frame.to_csv(input_path, index=False, encoding="utf-8")

    output_path = target.process_targets(input_path)
    result = pd.read_csv(output_path)

    assert result["isoform_name"].tolist() == [
        "Alpha",
        "Alpha",
        "Alpha",
        "Alpha",
        "Beta",
        "Beta",
        "Beta",
        "Beta",
    ]


@pytest.mark.unit
def test_process_targets__writes_to_custom_output_dir(tmp_path: Path):
    input_path = tmp_path / "targets.csv"
    frame = pd.DataFrame(
        [
            {
                "target_chembl_id": "CHEMBL_T5",
                "isoform_ids": "ENSP0006",
                "isoform_names": "Delta",
                "isoform_synonyms": "Delta",
            }
        ]
    )
    frame.to_csv(input_path, index=False, encoding="utf-8")

    destination = tmp_path / "out"
    output_path = target.process_targets(input_path, output_dir=destination)

    assert output_path.parent == destination
    assert output_path.exists()


@pytest.mark.unit
def test_process_targets__supports_custom_separator(tmp_path: Path):
    input_path = tmp_path / "targets.csv"
    buffer = io.StringIO()
    frame = pd.DataFrame(
        [
            {
                "target_chembl_id": "CHEMBL_T6",
                "isoform_ids": "ENSP0007",
                "isoform_names": "Epsilon",
                "isoform_synonyms": "Epsilon",
            }
        ]
    )
    frame.to_csv(buffer, index=False, sep=";")
    input_path.write_text(buffer.getvalue(), encoding="utf-8")

    output_path = target.process_targets(input_path, sep=";")
    result = pd.read_csv(output_path)

    assert result["isoform_name"].tolist() == ["Epsilon", "Epsilon"]


@pytest.mark.unit
def test_process_targets__handles_isoform_mismatch_lists(tmp_path: Path):
    input_path = tmp_path / "targets.csv"
    frame = pd.DataFrame(
        [
            {
                "target_chembl_id": "CHEMBL_T7",
                "isoform_ids": "ENSP0008|ENSP0009",
                "isoform_names": "Zeta|",
                "isoform_synonyms": "Zeta Alt:Zeta|",
            }
        ]
    )
    frame.to_csv(input_path, index=False, encoding="utf-8")

    output_path = target.process_targets(input_path)
    result = pd.read_csv(output_path)

    assert result.to_dict(orient="records") == [
        {
            "target_chembl_id": "CHEMBL_T7",
            "isoform_id": "ENSP0008",
            "isoform_name": "Zeta",
            "term": "Zeta",
            "token": "zeta",
        },
        {
            "target_chembl_id": "CHEMBL_T7",
            "isoform_id": "ENSP0008",
            "isoform_name": "Zeta",
            "term": "zeta",
            "token": "zeta",
        },
        {
            "target_chembl_id": "CHEMBL_T7",
            "isoform_id": "ENSP0008",
            "isoform_name": "Zeta",
            "term": "zeta alt",
            "token": "alt",
        },
        {
            "target_chembl_id": "CHEMBL_T7",
            "isoform_id": "ENSP0008",
            "isoform_name": "Zeta",
            "term": "zeta alt",
            "token": "zeta",
        },
    ]
