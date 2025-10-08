"""Unit tests for :mod:`library.cli.parser`."""

from __future__ import annotations

import argparse
from types import SimpleNamespace
from pathlib import Path

import pytest

from library.cli import parser as parser_module
from library.cli.parser import apply_config_overrides
from library.config import ConfigMetadata


class _SpyLogger:
    """Test double capturing warning calls."""

    def __init__(self) -> None:
        self.events: list[tuple[str, dict[str, object]]] = []

    def warning(self, event: str, **kwargs: object) -> None:  # pragma: no cover - thin
        self.events.append((event, kwargs))

    # Compatibility methods used elsewhere in the module.  They are no-ops for
    # the scope of these tests.
    def info(self, *args: object, **kwargs: object) -> None:  # pragma: no cover
        return None

    def error(self, *args: object, **kwargs: object) -> None:  # pragma: no cover
        return None

    def bind(self, **kwargs: object) -> "_SpyLogger":  # pragma: no cover
        return self


@pytest.mark.unit
def test_apply_config_overrides__missing_config_attribute(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """Gracefully ignore CLI mappings that point to absent config attributes."""

    parser = argparse.ArgumentParser()
    cfg_path = tmp_path / "config.yaml"
    cfg_path.write_text(
        "sources:\n  chembl:\n    pipelines:\n      target:\n        all:\n          column: target_chembl_id\n",
        encoding="utf-8",
    )
    parser.add_argument("--config", default=str(cfg_path))
    parser.add_argument("--column", default="target_chembl_id")
    args = parser.parse_args([])

    spy = _SpyLogger()
    monkeypatch.setattr(parser_module, "logger", spy)
    def _load_stub(*_args, **_kwargs):
        cfg = SimpleNamespace(
            sources=SimpleNamespace(
                chembl=SimpleNamespace(
                    pipelines=SimpleNamespace(
                        target=SimpleNamespace(
                            all=SimpleNamespace(column="target_chembl_id")
                        )
                    )
                )
            ),
            local=SimpleNamespace(io=SimpleNamespace(output_stamp_mode="omit")),
        )
        return cfg, ConfigMetadata(snapshot={}, sources={})

    monkeypatch.setattr(parser_module, "load_config", _load_stub)

    cfg = apply_config_overrides(
        args,
        parser,
        cfg_path,
        mapping={"column": "target.all.legacy_column"},
    )

    assert cfg.sources.chembl.pipelines.target.all.column == "target_chembl_id"
    assert args.column == "target_chembl_id"
    event_name = "config_attribute_missing"
    matching = [
        payload
        for name, payload in (spy.events)
        if name == event_name
        and payload.get("argument") == "column"
        and payload.get("path")
        == "sources.chembl.pipelines.target.all.legacy_column"
    ]
    assert matching, f"expected {event_name!r} log entry not found in {spy.events!r}"


@pytest.mark.unit
def test_apply_config_overrides__require_date_enforced(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    parser = argparse.ArgumentParser()
    cfg_path = tmp_path / "config.yaml"
    cfg_path.write_text(
        "local:\n  io:\n    output_stamp_mode: require\n",
        encoding="utf-8",
    )
    parser.add_argument("--config", default=str(cfg_path))
    parser.add_argument("--date", default=None)
    args = parser.parse_args([])

    def _fake_load_config(*_args, **_kwargs):
        cfg = SimpleNamespace(
            sources=SimpleNamespace(),
            local=SimpleNamespace(io=SimpleNamespace(output_stamp_mode="require")),
        )
        return cfg, ConfigMetadata(snapshot={}, sources={})

    monkeypatch.setattr(parser_module, "load_config", _fake_load_config)

    with pytest.raises(ValueError, match="--date is required"):
        apply_config_overrides(args, parser, cfg_path)
