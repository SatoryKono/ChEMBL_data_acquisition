from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest


def test_determinism_script_returns_zero() -> None:
    """Run determinism check script and assert a zero exit code.

    The script now writes a small DataFrame through both deterministic CSV
    writers and compares their SHA-256 hashes. When the output is deterministic
    the script exits with code ``0``.
    """

    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "library.utils.cli_tools.check_determinism",
            "--log-level",
            "DEBUG",
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stderr


def test_determinism_script_reports_mismatch(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The script should exit with ``1`` when hashes differ."""

    from library.utils.cli_tools import check_determinism as impl
    script = impl

    def fake_write_csv_deterministic(
        df: pd.DataFrame,
        path: Path,
        **_: object,
    ) -> Path:
        target = Path(path)
        content = "first" if "first" in target.name else "second"
        target.write_text(content, encoding="utf-8")
        return target

    monkeypatch.setattr(impl, "write_csv_deterministic", fake_write_csv_deterministic)
    monkeypatch.setattr(
        impl,
        "sha256_file",
        lambda p: "hash1" if "first" in str(p) else "hash2",
    )
    monkeypatch.setattr(sys, "argv", ["check_determinism", "--log-level", "DEBUG"])

    rc = script.main()
    assert rc == 1
