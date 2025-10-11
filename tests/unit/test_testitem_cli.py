from __future__ import annotations

from pathlib import Path

import pytest

from library.pipelines.testitem import cli as testitem_cli


@pytest.mark.unit
@pytest.mark.parametrize(
    "output, fallback_date, expected",
    [
        (Path(".output.testitem_20240101.csv.tmp"), "19991231", ("testitem", "20240101")),
        ("output.testitems_20240101.csv", "19991231", ("testitem", "20240101")),
        ("output.testitem.csv", "19991231", ("testitem", "19991231")),
    ],
)
def test_normalise_output_labels__derives_table_and_date(
    output: Path | str, fallback_date: str, expected: tuple[str, str]
) -> None:
    """Intermediate filenames should map to canonical artefact labels."""

    assert (
        testitem_cli._normalise_output_labels(output, fallback_date=fallback_date)  # noqa: SLF001
        == expected
    )
