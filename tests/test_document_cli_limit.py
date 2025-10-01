from pathlib import Path

import pytest

from scripts import get_document_data as gdd


def test_zero_limit_skips_execution(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """``--limit 0`` should short-circuit before configuring the pipeline."""

    input_csv = tmp_path / "documents.csv"
    input_csv.write_text("document_chembl_id\nCHEMBL1\n", encoding="utf8")
    output_csv = tmp_path / "out.csv"

    def fail_run(*_args: object, **_kwargs: object) -> int:
        raise AssertionError("run_chembl should not execute when limit is zero")

    monkeypatch.setattr(gdd, "run_chembl", fail_run)

    def fail_ensure(_: object) -> None:
        raise AssertionError("ensure_dirs should not run when limit is zero")

    monkeypatch.setattr(gdd, "ensure_dirs", fail_ensure)

    def fail_apply(*_args: object, **_kwargs: object) -> None:
        pytest.fail("apply_config_overrides should not run when limit is zero")

    monkeypatch.setattr(gdd.cli, "apply_config_overrides", fail_apply)

    exit_code = gdd.main(
        [
            "chembl",
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
