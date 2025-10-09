"""Unit tests for the document pipeline utilities."""

from __future__ import annotations

import pytest

from library.pipelines.document.pipeline import normalise_doi


@pytest.mark.unit
@pytest.mark.parametrize(
    ("value", "expected"),
    [
        (None, ""),
        ("", ""),
        ("  ", ""),
        ("10.1002/cmdc.202100547", "10.1002/cmdc.202100547"),
        ("DOI:10.1002/cmdc.202100547", "10.1002/cmdc.202100547"),
        ("doi:10.1002/cmdc.202100547", "10.1002/cmdc.202100547"),
        (
            "  https://doi.org/10.1002/cmdc.202100547  ",
            "10.1002/cmdc.202100547",
        ),
        (
            "http://doi.org/10.1002/cmdc.202100547",
            "10.1002/cmdc.202100547",
        ),
        (
            "doi.org/10.1002/cmdc.202100547",
            "10.1002/cmdc.202100547",
        ),
        (
            "10.1021/acs.jmedchem.5b00015",
            "10.1021/acs.jmedchem.5b00015",
        ),
        # Check that casing is preserved.
        ("10.1002/CMDc.202100547", "10.1002/CMDc.202100547"),
    ],
)
def test_normalise_doi(value: str | None, expected: str) -> None:
    """The ``normalise_doi`` helper should return canonical DOIs."""

    assert normalise_doi(value) == expected
