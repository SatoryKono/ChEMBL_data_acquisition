import pandas as pd
import pytest
import responses

from library import iuphar_library as ii
from library.config import IupharCfg


@responses.activate
def test_websearch_gene_to_id_uses_cfg(monkeypatch):
    data = ii.IUPHARData(target_df=pd.DataFrame(), family_df=pd.DataFrame())
    called: dict[str, tuple[int, int]] = {}
    orig_get = ii.requests.get

    def capture(url: str, timeout: tuple[int, int]):
        called["timeout"] = timeout
        return orig_get(url, timeout=timeout)

    sleeps: list[float] = []
    monkeypatch.setattr(ii.requests, "get", capture)
    monkeypatch.setattr(ii.time, "sleep", lambda s: sleeps.append(s))

    cfg = IupharCfg(
        base="https://example.org/services",
        timeout_connect=1,
        timeout_read=2,
        rps=5,
        burst=5,
    )
    url = "https://example.org/services/targets/?geneSymbol=GENE"
    responses.add(responses.GET, url, json=[{"id": 1}], status=200)
    result = data.websearch_gene_to_id("GENE", cfg)
    assert result == {"id": 1}
    assert responses.calls[0].request.url == url
    assert called["timeout"] == (1, 2)
    assert sleeps and sleeps[0] == pytest.approx(0.2)


@responses.activate
def test_iuphar_upload_uses_cfg(monkeypatch):
    target_df = pd.DataFrame({"target_id": [1], "family_id": ["F1"]})
    family_df = pd.DataFrame(
        {"family_id": ["F1"], "family_name": ["Fam"], "parent_family_id": [pd.NA]}
    )
    data = ii.IUPHARData(target_df=target_df, family_df=family_df)
    calls: list[tuple[str, tuple[int, int]]] = []
    orig_get = ii.requests.get

    uni_csv = "GtoPdb IUPHAR ID,IUPHAR ID,UniProtKB ID\n1,1,P12345\n"
    hgnc_csv = "GtoPdb IUPHAR ID,HGNC ID,IUPHAR Name\n1,HG1,Name\n"

    def capture(url: str, timeout: tuple[int, int]):
        calls.append((url, timeout))
        return orig_get(url, timeout=timeout)

    sleeps: list[float] = []
    monkeypatch.setattr(ii.requests, "get", capture)
    monkeypatch.setattr(ii.time, "sleep", lambda s: sleeps.append(s))

    cfg = IupharCfg(
        base="https://example.org/services",
        timeout_connect=1,
        timeout_read=2,
        rps=10,
        burst=5,
    )
    responses.add(
        responses.GET,
        "https://example.org/DATA/GtP_to_UniProt_mapping.csv",
        body=uni_csv,
        status=200,
    )
    responses.add(
        responses.GET,
        "https://example.org/DATA/GtP_to_HGNC_mapping.csv",
        body=hgnc_csv,
        status=200,
    )
    df = data.iuphar_upload(cfg)
    assert not df.empty
    assert calls[0][0] == "https://example.org/DATA/GtP_to_UniProt_mapping.csv"
    assert calls[0][1] == (1, 2)
    assert calls[1][0] == "https://example.org/DATA/GtP_to_HGNC_mapping.csv"
    assert sleeps and sleeps[0] == pytest.approx(0.1)
