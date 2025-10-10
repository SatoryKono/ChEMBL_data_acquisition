from __future__ import annotations

import pandas as pd
import pytest

from library.postprocessing.target.cellularity import (
    Cellularity,
    add_cellularity_smart,
)


@pytest.mark.unit
def test_classify_by_lineage__resolves_known_phyla() -> None:
    classifier = Cellularity(fetcher=lambda tax_id, email: [])
    assert classifier.classify_by_lineage("Eukaryota", "Chordata") == "multicellular"
    assert classifier.classify_by_lineage("Bacteria", "Proteobacteria") == "unicellular"
    assert classifier.classify_by_lineage("Viruses", "") == "acellular (virus)"


@pytest.mark.unit
def test_add_cellularity_smart__prefers_lineage_over_fetch() -> None:
    frame = pd.DataFrame(
        [
            {
                "taxon_id": "9606",
                "lineage_superkingdom": "Eukaryota",
                "lineage_class": "Chordata",
            }
        ]
    )
    result = add_cellularity_smart(
        frame,
        "taxon_id",
        "lineage_superkingdom",
        "lineage_class",
        fetcher=lambda tax_id, email: ["bacteria"],
    )
    assert result.loc[0, "cellularity"] == "multicellular"


@pytest.mark.unit
def test_add_cellularity_smart__falls_back_to_fetch(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    captured: dict[str, str | None] = {}

    def _fetch(tax_id: str, email: str | None) -> list[str]:
        captured["tax_id"] = tax_id
        captured["email"] = email
        return ["Eukaryota", "Metazoa"]

    frame = pd.DataFrame(
        [
            {
                "taxon_id": "12345",
                "lineage_superkingdom": "",
                "lineage_class": "",
            }
        ]
    )
    result = add_cellularity_smart(
        frame,
        "taxon_id",
        "lineage_superkingdom",
        "lineage_class",
        email="contact@example.org",
        fetcher=_fetch,
    )
    assert captured == {"tax_id": "12345", "email": "contact@example.org"}
    assert result.loc[0, "cellularity"] == "multicellular"
