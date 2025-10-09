"""Unit tests for :mod:`scripts.check_determinism`."""

from __future__ import annotations

import csv
import subprocess
from pathlib import Path
from subprocess import CompletedProcess

import pytest
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


def _write_run_payload(destination: Path, payload: str, metadata: dict[str, object]) -> None:
    destination.write_text(payload, encoding="utf-8")
    meta_path = destination.with_name(destination.name + ".meta.yaml")
    with meta_path.open("w", encoding="utf-8") as handle:
        yaml.safe_dump(metadata, handle, sort_keys=False)


def _patch_mkdtemp(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> list[Path]:
    created_dirs: list[Path] = []

    def _mkdtemp(prefix: str, dir: str | None = None) -> str:  # pragma: no cover - helper
        del dir
        temp_dir = tmp_path / f"{prefix}{len(created_dirs)}"
        temp_dir.mkdir()
        created_dirs.append(temp_dir)
        return str(temp_dir)

    monkeypatch.setattr(check_determinism.tempfile, "mkdtemp", _mkdtemp)
    return created_dirs


def test_main__metadata_equivalent_runs_pass(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    """Metadata sidecars differing only in timestamps should be accepted."""

    input_csv = tmp_path / "activity.csv"
    input_csv.write_text("activity_chembl_id\nCHEMBL1\n", encoding="utf-8")

    created_dirs = _patch_mkdtemp(tmp_path, monkeypatch)

    payload = "activity_id\n1\n"
    base_metadata = {
        "generated_at": "2025-01-01T00:00:00+00:00",
        "stats": {"rows_total": 1, "rows_kept": 1},
        "schema": "activity",
    }

    runs: list[tuple[int, Path, Path]] = []

    def _fake_run_activity(limit: int, destination: Path, observed_input: Path) -> CompletedProcess[str]:
        metadata = dict(base_metadata)
        metadata["generated_at"] = f"2025-01-01T00:00:0{len(runs)}+00:00"
        _write_run_payload(destination, payload, metadata)
        runs.append((limit, destination, observed_input))
        return CompletedProcess(args=["python"], returncode=0, stdout="ok\n", stderr="")

    monkeypatch.setattr(check_determinism, "_run_activity", _fake_run_activity)

    exit_code = check_determinism.main(["--limit", "2", "--input", str(input_csv)])

    assert exit_code == 0

    captured = capsys.readouterr()
    assert "Deterministic output confirmed" in captured.out

    assert len(runs) == 2
    for limit, _destination, observed_input in runs:
        assert limit == 2
        assert observed_input == input_csv

    for directory in created_dirs:
        assert not directory.exists()


def test_main__metadata_mismatch_returns_error(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    """Non-identical metadata must cause a determinism failure."""

    input_csv = tmp_path / "activity.csv"
    input_csv.write_text("activity_chembl_id\nCHEMBL1\n", encoding="utf-8")

    created_dirs = _patch_mkdtemp(tmp_path, monkeypatch)

    payload = "activity_id\n1\n"
    metadata_runs = [
        {
            "generated_at": "2025-01-01T00:00:00+00:00",
            "stats": {"rows_total": 1, "rows_kept": 1},
            "schema": "activity",
        },
        {
            "generated_at": "2025-01-01T00:00:01+00:00",
            "stats": {"rows_total": 2, "rows_kept": 2},
            "schema": "activity",
        },
    ]

    runs: list[tuple[int, Path, Path]] = []

    def _fake_run_activity(limit: int, destination: Path, observed_input: Path) -> CompletedProcess[str]:
        metadata = metadata_runs[len(runs)]
        _write_run_payload(destination, payload, metadata)
        runs.append((limit, destination, observed_input))
        return CompletedProcess(args=["python"], returncode=0, stdout="ok\n", stderr="")

    monkeypatch.setattr(check_determinism, "_run_activity", _fake_run_activity)

    exit_code = check_determinism.main(["--limit", "2", "--input", str(input_csv)])

    assert exit_code == 1

    captured = capsys.readouterr()
    assert "Metadata mismatch detected" in captured.out

    assert len(runs) == 2
    for limit, _destination, observed_input in runs:
        assert limit == 2
        assert observed_input == input_csv

    for directory in created_dirs:
        assert not directory.exists()
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
