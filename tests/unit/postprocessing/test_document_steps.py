from __future__ import annotations

import pandas as pd
import pytest

from library.postprocess.documents.steps import normalize_document_fields


@pytest.mark.unit
def test_normalize_document_fields__applies_whitespace_and_unicode() -> None:
    frame = pd.DataFrame(
        {
            "Title": ["  Cafe\u0301  "],
            "Journal": ["  Journal   of   Testing  "],
            "Doc_Type": ["  PREPRINT  "],
        }
    )

    result = normalize_document_fields(
        frame,
        trim_whitespace=True,
        normalise_unicode=True,
    )

    assert list(result.columns) == ["title", "journal", "doc_type"]
    assert result.loc[0, "title"] == "Café"
    assert result.loc[0, "journal"] == "Journal of Testing"
    assert result.loc[0, "doc_type"] == "PREPRINT"
    assert frame.loc[0, "Title"] == "  Cafe\u0301  "


@pytest.mark.unit
def test_normalize_document_fields__skips_whitespace_when_disabled() -> None:
    frame = pd.DataFrame({"Title": ["  Example  "]})

    result = normalize_document_fields(frame, trim_whitespace=False)

    assert result.loc[0, "title"] == "  Example  "
