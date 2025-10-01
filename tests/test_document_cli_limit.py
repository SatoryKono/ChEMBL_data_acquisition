from pathlib import Path

import pytest

from scripts import get_document_data as gdd


def test_zero_limit_skips_pipeline(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """``--limit 0`` should short-circuit the document pipeline."""

    recorded: list[tuple[str, dict[str, object]]] = []

    def capture_info(event: str, **kwargs: object) -> None:
        recorded.append((event, kwargs))

    monkeypatch.setattr(gdd.logger, "info", capture_info)

    def fail_run_cli_command(*_: object, **__: object) -> int:
        pytest.fail("run_cli_command must not execute when limit is zero")

    monkeypatch.setattr(gdd, "run_cli_command", fail_run_cli_command)

    output_path = tmp_path / "documents.csv"

    exit_code = gdd.main(
        ["chembl", "--limit", "0", "--output", str(output_path)]
    )

    assert exit_code == 0
    assert recorded == [("pipeline_skip_limit", {"limit": 0})]
    assert not output_path.exists()
    assert not Path(f"{output_path}.meta.yaml").exists()
