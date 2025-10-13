from __future__ import annotations

from pathlib import Path

import pytest

from library.cli.entrypoints import activity


@pytest.mark.unit
def test_derive_standard_output_labels__handles_hidden_temp_file(
    tmp_path: Path,
) -> None:
    dataset_path = tmp_path / ".output.activities_20240101.csv.tmp"

    table_name, date_tag = activity._derive_standard_output_labels(dataset_path)

    assert table_name == "activity"
    assert date_tag == "20240101"


@pytest.mark.unit
def test_derive_standard_output_labels__deduplicates_output_prefix(
    tmp_path: Path,
) -> None:
    dataset_path = tmp_path / "output..output.targets_20240101.csv"

    table_name, date_tag = activity._derive_standard_output_labels(dataset_path)

    assert table_name == "targets"
    assert date_tag == "20240101"


@pytest.mark.unit
def test_derive_standard_output_labels__falls_back_to_current_date(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    dataset_path = tmp_path / "summary.csv"
    monkeypatch.setattr(activity, "_current_date_token", lambda: "20991231")

    table_name, date_tag = activity._derive_standard_output_labels(dataset_path)

    assert table_name == "summary"
    assert date_tag == "20991231"


@pytest.mark.unit
def test_derive_standard_output_labels__uses_default_stem_when_missing(
    tmp_path: Path,
) -> None:
    dataset_path = tmp_path / ".output.20240101.csv.tmp"

    table_name, date_tag = activity._derive_standard_output_labels(dataset_path)

    assert table_name == activity.DEFAULT_OUTPUT_STEM
    assert date_tag == "20240101"


@pytest.mark.unit
def test_derive_standard_output_labels__normalises_dotted_stem(tmp_path: Path) -> None:
    dataset_path = tmp_path / "output.activity.20240101.csv"

    table_name, date_tag = activity._derive_standard_output_labels(dataset_path)

    assert table_name == activity.DEFAULT_OUTPUT_STEM
    assert date_tag == "20240101"
