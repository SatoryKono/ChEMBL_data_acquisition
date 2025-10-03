from __future__ import annotations

from pathlib import Path

import pytest

from scripts import get_activity_data as gad


def test_negative_limit_rejected(capsys: pytest.CaptureFixture[str]) -> None:
    """Ensure ``--limit`` rejects negative integers."""
    with pytest.raises(SystemExit) as excinfo:
        gad.main(["--limit", "-1"])
    assert excinfo.value.code == 2
    err = capsys.readouterr().err
    assert "--limit must be zero or a positive integer" in err


def test_zero_limit_skips_pipeline(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """``--limit 0`` should short-circuit execution without touching outputs."""

    recorded: list[tuple[str, dict[str, object]]] = []

    def capture_info(event: str, **kwargs: object) -> None:
        recorded.append((event, kwargs))

    monkeypatch.setattr(gad.logger, "info", capture_info)

    def fail_run_cli_command(*_: object, **__: object) -> int:
        pytest.fail("run_cli_command must not execute when limit is zero")

    monkeypatch.setattr(gad, "run_cli_command", fail_run_cli_command)

    output_path = tmp_path / "activities.csv"

    exit_code = gad.main(["--limit", "0", "--final-out", str(output_path)])

    assert exit_code == 0
    assert recorded == [("pipeline_skip_limit", {"limit": 0})]
    assert not output_path.exists()
    assert not Path(f"{output_path}.meta.yaml").exists()
