from __future__ import annotations

from typing import Any

import pytest

from library.config.env import (
    _apply_env_overrides,
    _expand_config_placeholders,
    _resolve_placeholder_base_path,
)


@pytest.mark.unit
def test_apply_env_overrides__alias(monkeypatch):
    monkeypatch.setenv("CHEMBL_DA_LOG_LEVEL", "DEBUG")
    data = {"system": {"log": {"level": "INFO"}}}

    overrides, warnings = _apply_env_overrides(data)

    assert data["system"]["log"]["level"] == "DEBUG"
    assert overrides[("system", "log", "level")] == "CHEMBL_DA_LOG_LEVEL"
    assert warnings == []


@pytest.mark.unit
def test_apply_env_overrides__explicit_path(monkeypatch):
    monkeypatch.setenv("CHEMBL_DA__SYSTEM__LOG__LEVEL", "WARNING")
    data = {"system": {"log": {"level": "INFO"}}}

    overrides, warnings = _apply_env_overrides(data)

    assert data["system"]["log"]["level"] == "WARNING"
    assert overrides[("system", "log", "level")] == "CHEMBL_DA__SYSTEM__LOG__LEVEL"
    assert warnings == []


@pytest.mark.unit
def test_apply_env_overrides__unknown_path_records_warning(monkeypatch, caplog):
    monkeypatch.setenv("CHEMBL_DA__SYSTEM__UNKNOWN", "42")
    data: dict[str, Any] = {}

    with caplog.at_level("WARNING"):
        overrides, warnings = _apply_env_overrides(data)

    assert overrides == {}
    assert warnings == [
        "CHEMBL_DA__SYSTEM__UNKNOWN: unknown configuration path 'system.unknown'"
    ]
    assert any("config_env_unknown_path" in message for message in caplog.messages)


@pytest.mark.unit
def test_resolve_placeholder_base_path__env_override(monkeypatch, tmp_path):
    base = tmp_path / "chembl"
    monkeypatch.setenv("CHEMBL_DA_BASE_PATH", str(base))

    resolved = _resolve_placeholder_base_path(None)

    assert resolved == base.resolve()


@pytest.mark.unit
def test_expand_config_placeholders__replaces_marker(tmp_path):
    base = tmp_path / "chembl"
    base.mkdir()
    data = {"paths": {"output": "$CHEMBL_DA_BASE_PATH/results"}}

    expanded = _expand_config_placeholders(data, base_path=base)

    assert expanded["paths"]["output"] == str(base / "results")
