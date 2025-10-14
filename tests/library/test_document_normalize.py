"""Unit tests covering document normalization helpers."""

from __future__ import annotations

import pandas as pd
import pytest

from library.pipelines.document.pipeline import normalise_doi
from library.postprocessing.documents.steps import enrich_document_publication_year


@pytest.mark.unit
@pytest.mark.parametrize(
    ("raw", "expected"),
    [
        (" 10.1234/example ", "10.1234/example"),
        ("https://doi.org/10.1000/abc", "10.1000/abc"),
        ("http://doi.org/10.2000/xyz", "10.2000/xyz"),
        ("doi.org/10.3000/test", "10.3000/test"),
        ("DOI:10.4000/demo", "10.4000/demo"),
        ("", ""),
    ],
)
def test_normalise_doi_strips_known_prefixes(raw: str, expected: str) -> None:
    """Normalization should remove DOI prefixes and whitespace consistently."""

    assert normalise_doi(raw) == expected


@pytest.mark.unit
def test_enrich_document_publication_year_prioritises_external_sources() -> None:
    """CrossRef and OpenAlex years should override the native export when valid."""

    frame = pd.DataFrame(
        {
            "year": ["", "2018", "invalid"],
            "crossref.published": ["2015", "", ""],
            "openalex.publication_year": [pd.NA, "2021", "2022"],
        }
    )

    result = enrich_document_publication_year(
        frame,
        fallback_year=1990,
        prefer_doi_year=True,
    )

    expected = pd.Series([2015, 2021, 2022], dtype="Int64")
    pd.testing.assert_series_equal(result["publication_year"], expected, check_names=False)
