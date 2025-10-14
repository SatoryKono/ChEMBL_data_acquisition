import numpy as np
import pandas as pd
import pytest

from library.postprocessing.target.cellularity import Cellularity


@pytest.mark.unit
def test_classify_by_fetch__fungal_phylum_without_literal_token():
    classifier = Cellularity(
        fetcher=lambda tax_id, email: ("Eukaryota", "Ascomycota"),
    )

    result = classifier.classify_by_fetch(tax_id=9606)

    assert result == "ambiguous"


@pytest.mark.unit
def test_get_lineage_names__series_missing_values_skip_fetch() -> None:
    calls: list[tuple[object, object | None]] = []

    def _fetcher(tax_id: object, email: str | None) -> list[str]:
        calls.append((tax_id, email))
        return ["unexpected"]

    classifier = Cellularity(fetcher=_fetcher)
    tax_id = pd.Series([pd.NA, np.nan], dtype="object")

    assert classifier.get_lineage_names(tax_id) == []
    assert calls == []


@pytest.mark.unit
def test_get_lineage_names__series_with_non_missing_triggers_fetch() -> None:
    calls: list[tuple[object, object | None]] = []

    def _fetcher(tax_id: object, email: str | None) -> list[str]:
        calls.append((tax_id, email))
        return ["value"]

    classifier = Cellularity(fetcher=_fetcher)
    tax_id = pd.Series([pd.NA, "9606"], dtype="object")

    assert classifier.get_lineage_names(tax_id) == ["value"]
    assert len(calls) == 1
    called_tax_id, called_email = calls[0]
    assert isinstance(called_tax_id, pd.Series)
    assert called_email is None
