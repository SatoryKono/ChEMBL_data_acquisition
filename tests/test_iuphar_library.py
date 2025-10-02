import pytest

from pathlib import Path

pd = pytest.importorskip("pandas")
requests = pytest.importorskip("requests")
responses = pytest.importorskip("responses")

from library import iuphar_library as ii  # noqa: E402
from library.clients import iuphar as ci  # noqa: E402
from library.config import IupharCfg, RetryCfg  # noqa: E402


def test_websearch_gene_to_id_uses_cfg2(monkeypatch) -> None:
    """``websearch_gene_to_id`` should honour the supplied configuration."""
    data = ii.IUPHARData(target_df=pd.DataFrame(), family_df=pd.DataFrame())
    called: dict[str, object] = {}

    def fake_get(url: str, timeout: tuple[int, int]) -> object:
        called["url"] = url
        called["timeout"] = timeout

        class Resp:
            def __enter__(self):
                return self

            def __exit__(self, *exc):
                return False

            def raise_for_status(self):
                return None

            def json(self) -> list[dict[str, int]]:
                return [{"id": 1}]

        return Resp()

    sleeps: list[float] = []
    monkeypatch.setattr(ci._session, "get", fake_get)

    monkeypatch.setattr(ci, "sleep", lambda s: sleeps.append(s))

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
    assert called["timeout"] == (1, 2)
    assert sleeps == []


@responses.activate
def test_iuphar_upload_uses_cfg(monkeypatch):
    target_df = pd.DataFrame({"target_id": ["001"], "family_id": ["F01"]})
    family_df = pd.DataFrame(
        {"family_id": ["F01"], "family_name": ["Fam"], "parent_family_id": [pd.NA]}
    )
    data = ii.IUPHARData(target_df=target_df, family_df=family_df)
    calls: list[tuple[str, tuple[int, int]]] = []
    uni_csv = "GtoPdb IUPHAR ID,IUPHAR ID,UniProtKB ID\n001,001,P12345\n"
    hgnc_csv = "GtoPdb IUPHAR ID,HGNC ID,IUPHAR Name\n001,HG01,Name\n"

    def fake_get(url: str, timeout: tuple[int, int]) -> object:
        calls.append((url, timeout))
        text = uni_csv if "UniProt" in url else hgnc_csv

        class Resp:
            def __init__(self, text: str) -> None:
                self.text = text

            def __enter__(self):
                return self

            def __exit__(self, *exc):
                return False

            def raise_for_status(self):
                return None

        return Resp(text)

    sleeps: list[float] = []
    monkeypatch.setattr(ci._session, "get", fake_get)

    monkeypatch.setattr(ci, "sleep", lambda s: sleeps.append(s))

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
    assert df["GtoPdb IUPHAR ID"].dropna().tolist() == ["001"]
    assert df["Target id"].dropna().tolist() == ["001"]
    assert df["chembl_hgnc_id"].dropna().tolist() == ["HG01"]
    assert sleeps == []


def test_query_gene_symbol_backoff(monkeypatch):
    cfg = IupharCfg(
        base="https://example.org/services",
        timeout_connect=1,
        timeout_read=2,
        rps=5,
        burst=5,
    )
    retry = RetryCfg(max_attempts=2, backoff_factor=1)
    calls: list[tuple[str, tuple[int, int]]] = []

    def fake_get(url: str, timeout: tuple[int, int]):
        calls.append((url, timeout))
        if len(calls) == 1:
            raise requests.RequestException("boom")

        class Resp:
            def __enter__(self):
                return self

            def __exit__(self, *exc):
                return False

            def raise_for_status(self):
                return None

            def json(self):
                return [{"id": 1}]

        return Resp()

    sleeps: list[float] = []
    monkeypatch.setattr(ci._session, "get", fake_get)
    monkeypatch.setattr(ci, "sleep", lambda s: sleeps.append(s))
    monkeypatch.setattr(ci.random, "uniform", lambda a, b: 0)

    result = ci.query_gene_symbol("GENE", cfg, retry)
    assert result == {"id": 1}
    assert calls[0][0] == "https://example.org/services/targets/?geneSymbol=GENE"
    assert sleeps == [pytest.approx(1.0)]


def test_iuphar_upload_retries(monkeypatch) -> None:
    """``iuphar_upload`` retries failed downloads with backoff."""

    target_df = pd.DataFrame({"target_id": ["001"], "family_id": ["F01"]})
    family_df = pd.DataFrame(
        {"family_id": ["F01"], "family_name": ["Fam"], "parent_family_id": [pd.NA]}
    )
    data = ii.IUPHARData(target_df=target_df, family_df=family_df)

    cfg = IupharCfg(
        base="https://example.org/services",
        timeout_connect=1,
        timeout_read=2,
        rps=10,
        burst=5,
    )
    retry = RetryCfg(max_attempts=2, backoff_factor=1)

    sleeps: list[float] = []
    monkeypatch.setattr(ci, "sleep", lambda s: sleeps.append(s))
    monkeypatch.setattr(ci.random, "uniform", lambda a, b: 0)

    calls = {"n": 0}

    def fake_get(url: str, timeout: tuple[int, int]):
        calls["n"] += 1
        if calls["n"] == 1:
            raise requests.RequestException("boom")

        text = (
            "GtoPdb IUPHAR ID,IUPHAR ID,UniProtKB ID\n001,001,P12345\n"
            if "UniProt" in url
            else "GtoPdb IUPHAR ID,HGNC ID,IUPHAR Name\n001,HG01,Name\n"
        )

        class Resp:
            def __init__(self, text: str) -> None:
                self.text = text

            def __enter__(self):
                return self

            def __exit__(self, *exc):
                return False

            def raise_for_status(self) -> None:
                return None

        return Resp(text)

    monkeypatch.setattr(ci._session, "get", fake_get)
    monkeypatch.setattr(
        ci,
        "get_limiter",
        lambda *a, **k: type("L", (), {"acquire": lambda self: None})(),
    )

    df = data.iuphar_upload(cfg, retry)
    assert df["Target id"].dropna().tolist() == ["001"]

    assert not df.empty
    assert sleeps == [pytest.approx(1.0)]
    assert calls["n"] == 3


def test_map_uniprot_file_uses_mapping_uniprot(tmp_path: Path) -> None:
    target_df = pd.DataFrame(
        {
            "target_id": ["T1"],
            "uniprot_id": ["Q99999"],
            "family_id": ["F1"],
            "type": ["Enzyme.Transferase"],
            "target_name": ["Example target"],
        }
    )
    family_df = pd.DataFrame(
        {
            "family_id": ["F1"],
            "parent_family_id": [pd.NA],
            "family_name": ["Example family"],
            "type": ["Enzyme.Transferase"],
        }
    )
    data = ii.IUPHARData(target_df=target_df, family_df=family_df)

    input_df = pd.DataFrame(
        {
            "uniprot_id": [""],
            "mapping_uniprot_id": ["P12345|Q99999"],
        }
    )
    input_csv = tmp_path / "input.csv"
    output_csv = tmp_path / "output.csv"
    input_df.to_csv(input_csv, index=False)

    result = data.map_uniprot_file(input_csv, output_csv)

    assert result.loc[0, "target_id"] == "T1"
    assert result.loc[0, "IUPHAR_class"] == "Enzyme"
    assert result.loc[0, "IUPHAR_subclass"] == "Transferase"
    assert result.loc[0, "mapping_uniprot_id"] == "P12345|Q99999"


def test_target_id_from_row_uses_mapping_accessions() -> None:
    target_df = pd.DataFrame(
        {
            "target_id": ["T1"],
            "uniprot_id": ["Q99999"],
            "family_id": ["F1"],
            "hgnc_id": ["HGNC:1"],
            "hgnc_name": ["GENE1"],
            "gene_name": ["GENE1"],
        }
    )
    family_df = pd.DataFrame(
        {
            "family_id": ["F1"],
            "parent_family_id": [pd.NA],
            "family_name": ["Example family"],
            "type": ["Enzyme.Transferase"],
        }
    )
    data = ii.IUPHARData(target_df=target_df, family_df=family_df)

    row = pd.Series({"uniprot_id": "", "mapping_uniprot_id": "ALT|Q99999"})

    assert data.target_id_from_row(row) == "T1"


def test_target_id_by_uniprot_accepts_iterable_values() -> None:
    target_df = pd.DataFrame(
        {
            "target_id": ["T1"],
            "uniprot_id": ["Q99999"],
            "family_id": ["F1"],
            "hgnc_id": ["HGNC:1"],
            "hgnc_name": ["GENE1"],
            "gene_name": ["GENE1"],
        }
    )
    family_df = pd.DataFrame(
        {
            "family_id": ["F1"],
            "parent_family_id": [pd.NA],
            "family_name": ["Example family"],
            "type": ["Enzyme.Transferase"],
        }
    )
    data = ii.IUPHARData(target_df=target_df, family_df=family_df)

    assert data.target_id_by_uniprot(["ALT", "Q99999"]) == "T1"
