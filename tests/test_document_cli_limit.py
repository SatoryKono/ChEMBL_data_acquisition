import pytest

from scripts import get_document_data as gdd


def test_zero_limit_allowed(monkeypatch: pytest.MonkeyPatch) -> None:
    """``--limit 0`` should be accepted by the document CLI."""

    called: dict[str, object] = {}

    def fake_run(cfg, args):  # type: ignore[no-untyped-def]
        called["cfg_limit"] = cfg.document.chembl.limit
        called["args_limit"] = args.limit
        return 0

    monkeypatch.setattr(gdd, "run_chembl", fake_run)
    monkeypatch.setattr(gdd, "ensure_dirs", lambda cfg: None)

    exit_code = gdd.main(["chembl", "--limit", "0"])

    assert exit_code == 0
    assert called["cfg_limit"] == 0
    assert called["args_limit"] == 0
