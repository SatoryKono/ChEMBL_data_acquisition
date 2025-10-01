from pathlib import Path

import pytest

from scripts import get_testitem_data as gtd


def test_zero_limit_skips_execution(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """``--limit 0`` should skip configuration and avoid creating artefacts."""

    def fail_run(*_args: object, **_kwargs: object) -> int:
        raise AssertionError("run_chembl should not execute when limit is zero")

    monkeypatch.setattr(gtd, "run_chembl", fail_run)

    def fail_ensure(_: object) -> None:
        raise AssertionError("ensure_dirs should not run when limit is zero")

    monkeypatch.setattr(gtd, "ensure_dirs", fail_ensure)

    def fail_apply(*_args: object, **_kwargs: object) -> None:
        pytest.fail("apply_config_overrides should not run when limit is zero")

    monkeypatch.setattr(gtd.cli, "apply_config_overrides", fail_apply)

    input_csv = tmp_path / "input.csv"
    input_csv.write_text("molecule_chembl_id\nCHEMBL1\n", encoding="utf8")
    output_csv = tmp_path / "out.csv"

    exit_code = gtd.main(
        [
            "--input",
            str(input_csv),
            "--output",
            str(output_csv),
            "--limit",
            "0",
        ]
    )

    assert exit_code == 0
    assert not output_csv.exists()
    assert not Path(f"{output_csv}.meta.yaml").exists()
    assert not output_csv.with_name(f"{output_csv.stem}_failure_cases.csv").exists()
    assert not output_csv.with_suffix(".quality.json").exists()
    stemless = output_csv.with_suffix("")
    assert not Path(f"{stemless}_quality_report_table.csv").exists()
    assert not Path(f"{stemless}_data_correlation_report_table.csv").exists()
