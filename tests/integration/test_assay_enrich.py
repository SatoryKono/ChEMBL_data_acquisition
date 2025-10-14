"""Integration tests covering dictionary enrichment for assays."""

from __future__ import annotations

from library.postprocessing.assay.steps import AssayFetcher, fetch_normalize_assay


def test_dictionary_enrichment_overrides_missing_values(tmp_path, monkeypatch):
    dictionary_dir = tmp_path / "_assay"
    dictionary_dir.mkdir()

    (dictionary_dir / "assay_group.csv").write_text(
        "assay_chembl_id,assay_group\nCHEMBL1,group-a\n",
        encoding="utf-8",
    )
    (dictionary_dir / "assay_strain.csv").write_text(
        "assay_chembl_id,assay_strain\nCHEMBL1,strain-a\n",
        encoding="utf-8",
    )
    (dictionary_dir / "accession.csv").write_text(
        "assay_chembl_id,accession\nCHEMBL1,P12345\n",
        encoding="utf-8",
    )

    records = [
        {
            "assay_chembl_id": "CHEMBL1",
            "assay_type": "BINDING",
            "assay_test_type": "IN VITRO",
            "target_chembl_id": "CHEMBL10",
            "created_on": "2020-01-01",
            "updated_on": None,
            "assay_group": None,
            "assay_strain": None,
            "accession": None,
        }
    ]

    def fake_fetch(self, limit=None):  # type: ignore[unused-argument]
        return list(records)

    monkeypatch.setattr(AssayFetcher, "fetch", fake_fetch, raising=False)

    frame = fetch_normalize_assay(dictionary_dir=dictionary_dir)

    row = frame.loc[0]
    assert row["assay_group"] == "group-a"
    assert row["assay_strain"] == "strain-a"
    assert row["accession"] == "P12345"
