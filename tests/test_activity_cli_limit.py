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


def test_zero_limit_skips_execution(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """``--limit 0`` should log a skip and avoid touching storage or configuration."""

    input_csv = tmp_path / "activity.csv"
    input_csv.write_text("activity_chembl_id\nCHEMBL1\n", encoding="utf8")
    output_csv = tmp_path / "out.csv"

    def fail_run(*_args: object, **_kwargs: object) -> int:
        raise AssertionError("run_chembl should not be invoked when limit is zero")

    monkeypatch.setattr(gad, "run_chembl", fail_run)

    def fail_ensure(_: object) -> None:
        raise AssertionError("ensure_dirs should not run when limit is zero")

    monkeypatch.setattr(gad, "ensure_dirs", fail_ensure)

    def fail_apply(*_args: object, **_kwargs: object) -> None:
        pytest.fail("apply_config_overrides should not run when limit is zero")

    monkeypatch.setattr(gad.cli, "apply_config_overrides", fail_apply)

    exit_code = gad.main(
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
