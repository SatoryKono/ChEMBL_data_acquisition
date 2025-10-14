from __future__ import annotations

from collections.abc import Iterable

import pandas as pd

from library.postprocessing.target.steps import (
    fetch_normalize_target,
    generate_target_reports,
)


class DummyChemblClient:
    def __init__(self, records: Iterable[dict[str, object]]) -> None:
        self._records = list(records)

    def list_targets(self, *, limit: int, offset: int = 0) -> dict[str, object]:
        slice_end = offset + limit
        chunk = self._records[offset:slice_end]
        next_url: str | None = None
        if slice_end < len(self._records):
            next_url = f"https://chembl/api/target.json?offset={slice_end}"
        return {"targets": chunk, "page_meta": {"next": next_url}}


class DummyUniProtClient:
    def __init__(self, families: dict[str, str], synonyms: dict[str, list[str]]) -> None:
        self._families = families
        self._synonyms = synonyms

    def get_entry(self, accession: str) -> dict[str, object]:
        return {
            "protein_families": [self._families.get(accession, "")],
            "synonyms": self._synonyms.get(accession, []),
        }


class DummyGtoPdbClient:
    def __init__(self, classes: dict[str, str]) -> None:
        self._classes = classes

    def get_target_by_uniprot(self, accession: str) -> dict[str, object]:
        return {"targetClass": self._classes.get(accession, "")}


def _build_target(index: int) -> dict[str, object]:
    chembl_id = f"CHEMBL{index:05d}"
    accession = f"P{index:05d}"
    return {
        "target_chembl_id": chembl_id,
        "pref_name": f"Target {index}",
        "target_type": "PROTEIN",
        "organism": "Homo sapiens",
        "target_components": [
            {
                "accession": accession,
                "gene_symbol": f"GENE{index}",
                "component_synonyms": [
                    {"syn_type": "GENE_SYMBOL", "synonym": f"GENE{index}"},
                    f"Synonym {index}",
                ],
            }
        ],
    }


def test_target_pipeline_enrichment() -> None:
    records = [_build_target(i) for i in range(1, 11)]
    chembl_client = DummyChemblClient(records)
    families = {f"P{i:05d}": f"Family {i}" for i in range(1, 11)}
    synonyms = {f"P{i:05d}": [f"Alias {i}"] for i in range(1, 11)}
    uniprot_client = DummyUniProtClient(families, synonyms)
    gtopdb_client = DummyGtoPdbClient({f"P{i:05d}": f"Class {i}" for i in range(1, 11)})

    dataset = fetch_normalize_target(
        10,
        chembl_client=chembl_client,
        uniprot_client=uniprot_client,
        gtopdb_client=gtopdb_client,
    )

    assert len(dataset) == 10
    assert dataset["target_class"].notna().mean() >= 0.9
    assert dataset["protein_family"].notna().mean() >= 0.9
    assert dataset["synonyms"].notna().mean() >= 0.9

    target_data = generate_target_reports(dataset)
    assert target_data.qc_summary["row_count"] == 10
    assert isinstance(target_data.quality_report, pd.DataFrame)
    assert isinstance(target_data.correlation_report, pd.DataFrame)
