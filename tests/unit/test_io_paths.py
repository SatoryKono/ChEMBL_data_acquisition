from __future__ import annotations

import pytest

from library.io.paths import derive_output_labels


@pytest.mark.unit
@pytest.mark.parametrize(
    "source, default_table, fallback_date, expected",
    [
        (
            ".output.activities_20240101.csv.tmp",
            "activities",
            "19700101",
            ("activities", "20240101"),
        ),
        (
            "output.activities.csv",
            "activities",
            "19991231",
            ("activities", "19991231"),
        ),
        (
            "OUTPUT.ASSAYS_20240101.csv",
            "assays",
            "19700101",
            ("ASSAYS", "20240101"),
        ),
        (
            "..output.output.targets_20200315.csv.tmp",
            "targets",
            "20000101",
            ("targets", "20200315"),
        ),
    ],
)
def test_derive_output_labels__normalises_source_names(
    source: str,
    default_table: str,
    fallback_date: str,
    expected: tuple[str, str],
) -> None:
    """The helper must normalise hidden stems and redundant prefixes."""

    assert derive_output_labels(
        source,
        default_table=default_table,
        fallback_date=fallback_date,
    ) == expected
