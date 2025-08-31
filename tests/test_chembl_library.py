import json
from pathlib import Path

import requests_mock

import sys
sys.path.append(str(Path(__file__).resolve().parents[1]))
from library import chembl_library as cl


def test_get_assays_returns_dataframe(requests_mock: requests_mock.Mocker) -> None:
    data_path = Path(__file__).parent / "data/assays.json"
    payload = json.loads(data_path.read_text())
    url = (
        "https://www.ebi.ac.uk/chembl/api/data/assay.json?format=json&assay_chembl_id__in=CHEMBL123"
    )
    requests_mock.get(url, json=payload)

    df = cl.get_assays(["CHEMBL123"])

    assert not df.empty
    assert df.loc[0, "assay_chembl_id"] == "CHEMBL123"


def test_get_assays_with_variants_adds_filter(requests_mock: requests_mock.Mocker) -> None:
    url = (
        "https://www.ebi.ac.uk/chembl/api/data/assay.json?format=json&variant_sequence__isnull=false&assay_chembl_id__in=CHEMBL123"
    )
    requests_mock.get(url, json={"assays": []})

    df = cl.get_assays_with_variants(["CHEMBL123"])

    assert df.empty
    assert requests_mock.called
    assert requests_mock.last_request.url == url


def test_get_assays_empty_ids_returns_empty_dataframe() -> None:
    df = cl.get_assays(["", "#N/A"])
    assert df.empty
