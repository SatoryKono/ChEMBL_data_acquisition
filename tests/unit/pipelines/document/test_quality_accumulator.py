from __future__ import annotations

import pandas as pd
import pytest

from library.pipelines.document.pipeline import DocumentQualityAccumulator


@pytest.mark.unit
def test_document_quality_accumulator_preserves_metrics() -> None:
    frame = pd.DataFrame(
        {
            "doi": ["10.1/abc", "", ""],
            "publication_class": ["review", None, " "],
            "PubMed.Error": ["", "timeout", None],
            "scholar.Error": [None, "", ""],
            "OpenAlex.Error": ["err", "", ""],
            "crossref.Error": ["", "", "x"],
            "fetch_status": ["error", "success", "error"],
            "error_source": ["PubMed", "Crossref", ""],
        }
    )

    accumulator = DocumentQualityAccumulator()
    accumulator.consume(frame)
    summary = accumulator.build()

    assert summary["rows_total"] == 3
    assert summary["doi_coverage"] == pytest.approx(1 / 3)
    assert summary["publication_class_counts"] == {"review": 1, "unknown": 2}
    assert summary["error_counts"] == {
        "pubmed": 2,
        "semantic_scholar": 1,
        "openalex": 1,
        "crossref": 1,
    }
    assert summary["error_placeholder_counts"] == {
        "pubmed": 1,
        "semantic_scholar": 0,
        "openalex": 0,
        "crossref": 0,
        "unknown": 1,
    }

