"""Unit tests for :mod:`scripts.check_determinism`."""

from __future__ import annotations

import csv
import subprocess
from pathlib import Path

import yaml

from scripts import check_determinism


def test_default_input_csv__matches_activity_column(tmp_path: Path) -> None:
    """Ensure the fallback CSV mirrors the configured activity column."""

    original_path_cls = check_determinism.Path
    script_repo_root = original_path_cls(check_determinism.__file__).resolve().parents[1]
    candidate_path = script_repo_root / "data" / "input" / "activity.csv"
    candidate_str = str(candidate_path)

    path_base = type(original_path_cls())

    class _TestPath(path_base):
        def __new__(cls, *args, **kwargs):  # pragma: no cover - delegating constructor
            return super().__new__(cls, *args, **kwargs)

        def exists(self) -> bool:  # pragma: no cover - exercised indirectly
            if str(self) == candidate_str:
                return False
            return super().exists()

    check_determinism.Path = _TestPath
    try:
        fallback = check_determinism._default_input_csv(tmp_path)
    finally:
        check_determinism.Path = original_path_cls

    assert fallback.parent == tmp_path

    with fallback.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.reader(handle))

    assert rows, "Fallback CSV must contain at least a header row"

    header, *data_rows = rows
    repo_root = Path(__file__).resolve().parents[2]
    config_path = repo_root / "config" / "config.yaml"
    with config_path.open(encoding="utf-8") as config_file:
        config_data = yaml.safe_load(config_file) or {}

    expected_column = (
        config_data
        .get("sources", {})
        .get("chembl", {})
        .get("pipelines", {})
        .get("activity", {})
        .get("column")
    )

    assert expected_column, "Activity pipeline column must be configured"

    assert header == [expected_column]
    assert data_rows == [["ACT1"], ["ACT2"]]


def test_run_activity__passes_dry_run_flag(monkeypatch, tmp_path: Path) -> None:
    """Ensure the activity runner forwards the --dry-run flag."""

    captured: dict[str, object] = {}

    def _fake_run(cmd, *, text, capture_output, env):
        captured["cmd"] = cmd
        captured["env"] = env
        return subprocess.CompletedProcess(cmd, 0, stdout="", stderr="")

    monkeypatch.setattr(check_determinism.subprocess, "run", _fake_run)

    destination = tmp_path / "out.csv"
    input_csv = tmp_path / "input.csv"

    check_determinism._run_activity(
        limit=5,
        destination=destination,
        input_csv=input_csv,
        dry_run=True,
    )

    assert "cmd" in captured, "subprocess.run must be invoked"
    assert captured["cmd"].count("--dry-run") == 1
