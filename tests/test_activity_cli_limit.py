import pytest

from scripts import get_activity_data as gad


def test_negative_limit_rejected(capsys: pytest.CaptureFixture[str]) -> None:
    """Ensure ``--limit`` rejects negative integers."""
    with pytest.raises(SystemExit) as excinfo:
        gad.main(["--limit", "-1"])
    assert excinfo.value.code == 2
    err = capsys.readouterr().err
    assert "--limit must be zero or a positive integer" in err


def test_zero_limit_allowed(monkeypatch: pytest.MonkeyPatch) -> None:
    """``--limit 0`` should succeed and propagate to configuration overrides."""

    called: dict[str, object] = {}

    def fake_run(cfg, args):  # type: ignore[no-untyped-def]
        called["cfg_limit"] = cfg.activity.limit
        called["args_limit"] = args.limit
        return 0

    monkeypatch.setattr(gad, "run_chembl", fake_run)
    monkeypatch.setattr(gad, "ensure_dirs", lambda cfg: None)

    exit_code = gad.main(["--limit", "0", "--dry-run"])

    assert exit_code == 0
    assert called["cfg_limit"] == 0
    assert called["args_limit"] == 0
