"""Unit tests covering logging setup for ``scripts.get_data``."""

from __future__ import annotations

from collections.abc import Callable
from pathlib import Path
from types import SimpleNamespace

import pytest
import scripts.get_data as get_data


@pytest.mark.unit
def test_configure_logging__creates_file_and_registers_shutdown(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """``configure_logging`` should reuse CLI helpers for file creation."""

    monkeypatch.setattr(get_data, "LOGS_DIR", tmp_path)
    monkeypatch.setattr(get_data, "_LOGGING_CONTEXT_MANAGER", None)
    monkeypatch.setattr(get_data, "_LOGGING_CONTEXT", None)
    monkeypatch.setattr(get_data, "_LOGGING_SHUTDOWN_REGISTERED", False)

    registered_callbacks: list[Callable[[], None]] = []

    def _register(callback: Callable[[], None]) -> Callable[[], None]:
        registered_callbacks.append(callback)
        return callback

    monkeypatch.setattr(
        get_data,
        "atexit",
        SimpleNamespace(register=_register),
    )

    created_configs: list[object] = []
    original_create = get_data.create_logger_config

    def _capture(level: str, *, run_id: str | None = None):
        cfg = original_create(level, run_id=run_id)
        created_configs.append(cfg)
        return cfg

    monkeypatch.setattr(get_data, "create_logger_config", _capture)

    log_path = get_data.configure_logging("warning")

    assert log_path.parent == tmp_path
    assert log_path.exists()
    assert created_configs, "expected create_logger_config to be invoked"
    cfg = created_configs[0]
    assert cfg.level == "WARNING"
    assert cfg.run_id, "run identifier should be populated"
    assert registered_callbacks, "shutdown hook should be registered"

    for callback in registered_callbacks:
        callback()

    assert get_data._LOGGING_CONTEXT_MANAGER is None
    assert get_data._LOGGING_SHUTDOWN_REGISTERED is False
