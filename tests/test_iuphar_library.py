import pandas as pd
import pytest

from library import iuphar_library as ii
from library.config import IupharCfg


def test_websearch_gene_to_id_uses_cfg(monkeypatch) -> None:
    """``websearch_gene_to_id`` should honour the supplied configuration."""
    data = ii.IUPHARData(target_df=pd.DataFrame(), family_df=pd.DataFrame())
    called: dict[str, object] = {}

    def fake_get(url: str, timeout: tuple[int, int]) -> object:
        called["url"] = url
        called["timeout"] = timeout

        class Resp:
            def raise_for_status(self) -> None:
                return None

            def json(self) -> list[dict[str, int]]:
                return [{"id": 1}]

        return Resp()

    sleeps: list[float] = []
    monkeypatch.setattr(ii.requests, "get", fake_get)
    monkeypatch.setattr(ii.time, "sleep", lambda s: sleeps.append(s))

    cfg = IupharCfg(
        base="https://example.org/services",
        timeout_connect=1,
        timeout_read=2,
        rps=5,
        burst=5,
    )
    result = data.websearch_gene_to_id("GENE", cfg)
    assert result == {"id": 1}
    assert called["url"] == "https://example.org/services/targets/?geneSymbol=GENE"
    assert called["timeout"] == (1, 2)
    assert sleeps and sleeps[0] == pytest.approx(0.2)


def test_iuphar_upload_uses_cfg(monkeypatch) -> None:
    """Uploading auxiliary data should use configured HTTP timeouts."""
    target_df = pd.DataFrame({"target_id": [1], "family_id": ["F1"]})
    family_df = pd.DataFrame(
        {"family_id": ["F1"], "family_name": ["Fam"], "parent_family_id": [pd.NA]}
    )
    data = ii.IUPHARData(target_df=target_df, family_df=family_df)
    calls: list[tuple[str, tuple[int, int]]] = []

    uni_csv = "GtoPdb IUPHAR ID,IUPHAR ID,UniProtKB ID\n1,1,P12345\n"
    hgnc_csv = "GtoPdb IUPHAR ID,HGNC ID,IUPHAR Name\n1,HG1,Name\n"

    def fake_get(url: str, timeout: tuple[int, int]) -> object:
        calls.append((url, timeout))
        text = uni_csv if "UniProt" in url else hgnc_csv

        class Resp:
            def __init__(self, text: str) -> None:
                self.text = text

            def raise_for_status(self) -> None:
                return None

        return Resp(text)

    sleeps: list[float] = []
    monkeypatch.setattr(ii.requests, "get", fake_get)
    monkeypatch.setattr(ii.time, "sleep", lambda s: sleeps.append(s))

    cfg = IupharCfg(
        base="https://example.org/services",
        timeout_connect=1,
        timeout_read=2,
        rps=10,
        burst=5,
    )
    df = data.iuphar_upload(cfg)
    assert not df.empty
    assert calls[0][0] == "https://example.org/DATA/GtP_to_UniProt_mapping.csv"
    assert calls[0][1] == (1, 2)
    assert calls[1][0] == "https://example.org/DATA/GtP_to_HGNC_mapping.csv"
    assert sleeps and sleeps[0] == pytest.approx(0.1)
