import pytest

from scripts import get_assay_data as gas


def test_zero_limit_allowed(monkeypatch: pytest.MonkeyPatch) -> None:
    """``--limit 0`` should be accepted by the assay CLI."""

    called: dict[str, object] = {}

    def fake_run(cfg, args):  # type: ignore[no-untyped-def]
        called["cfg_limit"] = cfg.assay.limit
        called["args_limit"] = args.limit
        return 0

    monkeypatch.setattr(gas, "run_chembl", fake_run)
    monkeypatch.setattr(gas, "ensure_dirs", lambda cfg: None)

    exit_code = gas.main(["--limit", "0"])

    assert exit_code == 0
    assert called["cfg_limit"] == 0
    assert called["args_limit"] == 0
