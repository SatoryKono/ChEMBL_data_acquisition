"""Unit tests for :mod:`library.cli.parser`."""

from __future__ import annotations

import argparse
from pathlib import Path

import pytest

from library.cli import parser as parser_module
from library.cli.parser import apply_config_overrides


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

    def bind(self, **kwargs: object) -> _SpyLogger:  # pragma: no cover
        return self


@pytest.mark.unit
def test_apply_config_overrides__missing_config_attribute(monkeypatch: pytest.MonkeyPatch) -> None:
    """Gracefully ignore CLI mappings that point to absent config attributes."""

    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config/config.yaml")
    parser.add_argument("--column", default="target_chembl_id")
    args = parser.parse_args([])

    spy = _SpyLogger()
    monkeypatch.setattr(parser_module, "logger", spy)

    cfg = apply_config_overrides(
        args,
        parser,
        Path("config/config.yaml"),
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
