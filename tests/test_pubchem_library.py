from library import pubchem_library as pl


def test_get_cid_from_smiles_uses_base(monkeypatch):
    called = {}

    def fake_request(url: str, cfg: pl.PubChemCfg, delay: float = 3.0):
        called["url"] = url
        return {"IdentifierList": {"CID": [1]}}

    monkeypatch.setattr(pl, "make_request", fake_request)
    cfg = pl.PubChemCfg(base="https://example.org/api")
    cid = pl.get_cid_from_smiles("C", cfg)
    assert called["url"] == "https://example.org/api/compound/smiles/C/cids/JSON"
    assert cid == "1"


def test_make_request_uses_timeout(monkeypatch):
    called = {}

    def fake_get(url: str, timeout: tuple[int, int]):
        called["timeout"] = timeout

        class Resp:
            status_code = 200

            def raise_for_status(self):
                return None

            def json(self):
                return {}

        return Resp()

    monkeypatch.setattr(pl._session, "get", fake_get)
    monkeypatch.setattr(pl.time, "sleep", lambda s: None)
    cfg = pl.PubChemCfg(timeout_connect=1, timeout_read=2)
    pl.make_request("https://example.org", cfg, delay=0)
    assert called["timeout"] == (1, 2)
