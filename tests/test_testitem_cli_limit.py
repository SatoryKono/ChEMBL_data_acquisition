from pathlib import Path

import pytest

from scripts import get_testitem_data as gtd


def test_zero_limit_allowed(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """``--limit 0`` should be accepted by the test item CLI."""

    called: dict[str, object] = {}

    def fake_run(cfg, args):  # type: ignore[no-untyped-def]
        called["cfg_limit"] = cfg.testitem.limit
        called["args_limit"] = args.limit
        return 0

    monkeypatch.setattr(gtd, "run_chembl", fake_run)
    monkeypatch.setattr(gtd, "ensure_dirs", lambda cfg: None)

    input_csv = tmp_path / "input.csv"
    input_csv.write_text("molecule_chembl_id\nCHEMBL1\n")
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
    assert called["cfg_limit"] == 0
    assert called["args_limit"] == 0
