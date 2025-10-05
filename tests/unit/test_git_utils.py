from __future__ import annotations

import importlib
import logging
from pathlib import Path

import pytest


@pytest.mark.unit
@pytest.mark.parametrize("module_name", ["library.git_utils", "library.common.git"])
def test_git_sha__missing_git_directory_logged_as_info(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    caplog: pytest.LogCaptureFixture,
    module_name: str,
) -> None:
    """Missing ``.git`` directories should not emit warning-level logs."""

    module = importlib.import_module(module_name)

    cache_clear = getattr(module._git_sha, "cache_clear", None)
    if cache_clear is None:  # pragma: no cover - defensive guard
        pytest.skip("_git_sha missing cache_clear helper")
    cache_clear()

    monkeypatch.setattr(module, "_repo_root", lambda: tmp_path)

    caplog.set_level(logging.DEBUG, logger="chembl")

    sha = module._git_sha()

    assert sha == "UNKNOWN"

    directory_logs = [
        (record.levelno, record.message)
        for record in caplog.records
        if "git_directory_missing" in record.message
    ]

    assert directory_logs, "expected git_directory_missing log entry"
    assert all(level < logging.WARNING for level, _ in directory_logs)
