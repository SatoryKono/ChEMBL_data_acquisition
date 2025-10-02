from pathlib import Path

import pytest

from scripts import get_testitem_data as gtd


def test_zero_limit_skips_pipeline(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """``--limit 0`` should short-circuit the test item pipeline."""

    recorded: list[tuple[str, dict[str, object]]] = []

    def capture_info(event: str, **kwargs: object) -> None:
        recorded.append((event, kwargs))

    monkeypatch.setattr(gtd.logger, "info", capture_info)

    def fail_run_cli_command(*_: object, **__: object) -> int:
        pytest.fail("run_cli_command must not execute when limit is zero")

    monkeypatch.setattr(gtd, "run_cli_command", fail_run_cli_command)

    input_csv = tmp_path / "input.csv"
    input_csv.write_text("molecule_chembl_id\nCHEMBL1\n")
    output_csv = tmp_path / "out.csv"

    exit_code = gtd.main(
        [
            "--input",
            str(input_csv),
            "--final-out",
            str(output_csv),
            "--limit",
            "0",
        ]
    )

    assert exit_code == 0
    assert recorded == [("pipeline_skip_limit", {"limit": 0})]
    assert not output_csv.exists()
    assert not Path(f"{output_csv}.meta.yaml").exists()
