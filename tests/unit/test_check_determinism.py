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


def test_main__metadata_sidecars_verified(
    tmp_path: Path, monkeypatch, capsys
) -> None:
    """Report success when metadata sidecars are identical."""

    work_dir = tmp_path / "runs"
    work_dir.mkdir()
    monkeypatch.setattr(
        check_determinism.tempfile, "mkdtemp", lambda prefix: str(work_dir)
    )

    input_csv = tmp_path / "input.csv"
    input_csv.write_text("activity_chembl_id\nACT0\n", encoding="utf-8")

    meta_contents = iter(["key: value\n", "key: value\n"])

    def _fake_run(limit: int, destination: Path, input_file: Path):
        destination.write_text(
            "activity_chembl_id\nACT0\n", encoding="utf-8"
        )
        meta_path = destination.with_suffix(destination.suffix + ".meta.yaml")
        meta_path.write_text(next(meta_contents), encoding="utf-8")
        return subprocess.CompletedProcess(
            args=["test"], returncode=0, stdout="", stderr=""
        )

    monkeypatch.setattr(check_determinism, "_run_activity", _fake_run)

    exit_code = check_determinism.main(["--limit", "1", "--input", str(input_csv)])
    captured = capsys.readouterr()

    assert exit_code == 0
    assert "Metadata deterministic output confirmed" in captured.out
    assert "Metadata SHA256:" in captured.out


def test_main__metadata_sidecars_mismatch(
    tmp_path: Path, monkeypatch, capsys
) -> None:
    """Fail with a warning when metadata sidecars diverge."""

    work_dir = tmp_path / "runs"
    work_dir.mkdir()
    monkeypatch.setattr(
        check_determinism.tempfile, "mkdtemp", lambda prefix: str(work_dir)
    )

    input_csv = tmp_path / "input.csv"
    input_csv.write_text("activity_chembl_id\nACT0\n", encoding="utf-8")

    meta_contents = iter(["key: first\n", "key: second\n"])

    def _fake_run(limit: int, destination: Path, input_file: Path):
        destination.write_text(
            "activity_chembl_id\nACT0\n", encoding="utf-8"
        )
        meta_path = destination.with_suffix(destination.suffix + ".meta.yaml")
        meta_path.write_text(next(meta_contents), encoding="utf-8")
        return subprocess.CompletedProcess(
            args=["test"], returncode=0, stdout="", stderr=""
        )

    monkeypatch.setattr(check_determinism, "_run_activity", _fake_run)

    exit_code = check_determinism.main(["--limit", "1", "--input", str(input_csv)])
    captured = capsys.readouterr()

    assert exit_code == 1
    assert "WARNING: Metadata mismatch detected:" in captured.out
