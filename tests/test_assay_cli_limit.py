from pathlib import Path

import pytest

from scripts import get_assay_data as gas


def test_zero_limit_skips_pipeline(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """``--limit 0`` should short-circuit the assay pipeline."""

    recorded: list[tuple[str, dict[str, object]]] = []

    def capture_info(event: str, **kwargs: object) -> None:
        recorded.append((event, kwargs))

    monkeypatch.setattr(gas.logger, "info", capture_info)

    def fail_run_cli_command(*_: object, **__: object) -> int:
        pytest.fail("run_cli_command must not execute when limit is zero")

    monkeypatch.setattr(gas, "run_cli_command", fail_run_cli_command)

    output_path = tmp_path / "assays.csv"

    exit_code = gas.main(["--limit", "0", "--output", str(output_path)])

    assert exit_code == 0
    assert recorded == [("pipeline_skip_limit", {"limit": 0})]
    assert not output_path.exists()
    assert not Path(f"{output_path}.meta.yaml").exists()
