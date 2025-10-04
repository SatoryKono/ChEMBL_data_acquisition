import pytest

from library.postprocessing.target.cellularity import Cellularity


@pytest.mark.unit
def test_classify_by_fetch__fungal_phylum_without_literal_token():
    classifier = Cellularity(
        fetcher=lambda tax_id, email: ("Eukaryota", "Ascomycota"),
    )

    result = classifier.classify_by_fetch(tax_id=9606)

    assert result == "ambiguous"
