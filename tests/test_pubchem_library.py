import responses

from library import pubchem_library as pl


@responses.activate
def test_get_cid_from_smiles_uses_base() -> None:
    cfg = pl.PubChemCfg(base="https://example.org/api", delay=0)
    url = "https://example.org/api/compound/smiles/C/cids/JSON"
    responses.add(responses.GET, url, json={"IdentifierList": {"CID": [1]}}, status=200)
    cid = pl.get_cid_from_smiles("C", cfg)
    assert responses.calls[0].request.url == url
    assert cid == "1"


@responses.activate
def test_make_request_uses_timeout(monkeypatch) -> None:
    called: dict[str, tuple[int, int]] = {}
    orig_get = pl._session.get

    def capture(url: str, timeout: tuple[int, int]):
        called["timeout"] = timeout
        return orig_get(url, timeout=timeout)

    monkeypatch.setattr(pl._session, "get", capture)
    monkeypatch.setattr(pl.time, "sleep", lambda s: None)
    cfg = pl.PubChemCfg(timeout_connect=1, timeout_read=2, delay=0, retries=1)
    responses.add(responses.GET, "https://example.org", json={}, status=200)
    pl.make_request("https://example.org", cfg)
    assert called["timeout"] == (1, 2)
