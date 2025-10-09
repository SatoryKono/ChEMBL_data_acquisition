"""Integration tests for :mod:`scripts.check_determinism`."""

from __future__ import annotations

import hashlib
from pathlib import Path
from subprocess import CompletedProcess

import pytest

from scripts import check_determinism


@pytest.mark.integration
def test_main__reports_deterministic_hash(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    """The determinism CLI should detect matching hashes across repeated runs."""

    input_csv = tmp_path / "activity.csv"
    input_csv.write_text("activity_chembl_id\nCHEMBL1\n", encoding="utf-8")

    created_dirs: list[Path] = []

    def _mkdtemp(
        prefix: str, dir: str | None = None
    ) -> str:  # pragma: no cover - helper
        del dir
        temp_dir = tmp_path / f"{prefix}{len(created_dirs)}"
        temp_dir.mkdir()
        created_dirs.append(temp_dir)
        return str(temp_dir)

    monkeypatch.setattr(check_determinism.tempfile, "mkdtemp", _mkdtemp)

    payload = "activity_id\n1\n2\n"
    expected_hash = hashlib.sha256(payload.encode("utf-8")).hexdigest()

    runs: list[tuple[int, Path, Path]] = []

    def _fake_run_activity(
        limit: int, destination: Path, input_path: Path, *, dry_run: bool
    ) -> CompletedProcess[str]:
        assert dry_run is True
        destination.write_text(payload, encoding="utf-8")
        runs.append((limit, destination, input_path))
        return CompletedProcess(args=["python"], returncode=0, stdout="ok\n", stderr="")

    monkeypatch.setattr(check_determinism, "_run_activity", _fake_run_activity)

    exit_code = check_determinism.main(["--limit", "2", "--input", str(input_csv)])

    assert exit_code == 0

    captured = capsys.readouterr()
    assert "Deterministic output confirmed" in captured.out
    assert "Metadata hash check: skipped" in captured.out
    assert f"SHA256: {expected_hash}" in captured.out

    assert len(runs) == 2
    for limit, _destination, observed_input in runs:
        assert limit == 2
        assert observed_input == input_csv

    for directory in created_dirs:
        assert not directory.exists()
