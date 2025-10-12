"""Unit tests for :mod:`library.cli.parser`."""

from __future__ import annotations

import argparse
from pathlib import Path
from types import SimpleNamespace
from uuid import UUID

import pytest

from library.cli import parser as parser_module
from library.cli.parser import apply_config_overrides
from library.config import ConfigMetadata
from library.io.paths import default_output_path


class _SpyLogger:
    """Test double capturing warning calls."""

    def __init__(self) -> None:
        self.events: list[tuple[str, dict[str, object]]] = []

    def warning(self, event: str, **kwargs: object) -> None:  # pragma: no cover - thin
        self.events.append((event, kwargs))

    # Compatibility methods used elsewhere in the module.  They record events so
    # assertions can inspect the emitted structured logs.
    def info(self, event: str, **kwargs: object) -> None:  # pragma: no cover
        self.events.append((event, kwargs))

    def error(self, *args: object, **kwargs: object) -> None:  # pragma: no cover
        return None

    def bind(self, **kwargs: object) -> _SpyLogger:  # pragma: no cover
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
        and payload.get("path") == "sources.chembl.pipelines.target.all.legacy_column"
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

    cfg = apply_config_overrides(args, parser, cfg_path)

    assert cfg.local.io.output_stamp_mode == "require"
    assert args.output_stamp_mode == "require"
    assert getattr(args, "date", None) is None


@pytest.mark.unit
def test_apply_config_overrides__uses_default_config_when_none(monkeypatch):
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default=None)
    parser.add_argument("--log-level", default="INFO")
    args = parser.parse_args([])

    spy = _SpyLogger()
    monkeypatch.setattr(parser_module, "logger", spy)

    calls: dict[str, object] = {}

    def _load_stub(path, **_kwargs):
        calls["path"] = path
        cfg = SimpleNamespace(
            sources=SimpleNamespace(),
            local=SimpleNamespace(io=SimpleNamespace(output_stamp_mode="omit")),
        )
        return cfg, ConfigMetadata(snapshot={}, sources={})

    monkeypatch.setattr(parser_module, "load_config", _load_stub)

    cfg = apply_config_overrides(args, parser, None)

    assert cfg.local.io.output_stamp_mode == "omit"
    assert calls["path"] == parser_module.DEFAULT_CONFIG_PATH
    assert args.config == parser_module.DEFAULT_CONFIG_PATH
    assert any(
        event == "config_default_path_used"
        and payload.get("config") == str(parser_module.DEFAULT_CONFIG_PATH)
        for event, payload in spy.events
    )


@pytest.mark.unit
def test_create_logger_config__default_run_id_unique(monkeypatch):
    """Ensure default run IDs remain unique across invocations."""

    uuid_values = iter(
        [
            UUID("00000000-0000-0000-0000-000000000001"),
            UUID("00000000-0000-0000-0000-000000000002"),
        ]
    )

    monkeypatch.setattr(parser_module.uuid, "uuid4", lambda: next(uuid_values))

    first = parser_module.create_logger_config("info")
    second = parser_module.create_logger_config("info")

    assert first.run_id != second.run_id
    # Both identifiers should remain parseable as UUIDs for compatibility.
    UUID(first.run_id)
    UUID(second.run_id)


@pytest.mark.unit
def test_create_logger_config__explicit_run_id_retained():
    """User-specified identifiers must bypass the default generation."""

    cfg = parser_module.create_logger_config("info", run_id="custom-id")

    assert cfg.run_id == "custom-id"


@pytest.mark.unit
def test_default_output_path__uses_config_prefix(tmp_path: Path) -> None:
    cfg = SimpleNamespace(
        output_dir=tmp_path,
        default_date_prefix="19991231",
    )

    output_path = default_output_path(tmp_path / "input.csv", cfg)

    assert output_path == tmp_path / "output.input_19991231.csv"


@pytest.mark.unit
def test_default_output_path__prefers_cli_date_over_config_prefix(tmp_path: Path) -> None:
    cfg = SimpleNamespace(
        output_dir=tmp_path,
        default_date_prefix="19991231",
    )

    output_path = default_output_path(
        tmp_path / "input.csv",
        cfg,
        date="20240615",
    )

    assert output_path == tmp_path / "output.input_20240615.csv"


@pytest.mark.unit
def test_default_output_path__falls_back_to_explicit_date(tmp_path: Path) -> None:
    cfg = SimpleNamespace(output_dir=tmp_path, default_date_prefix=None)

    output_path = default_output_path(
        tmp_path / "input.csv",
        cfg,
        date="20240615",
    )

    assert output_path == tmp_path / "output.input_20240615.csv"


@pytest.mark.unit
def test_default_output_path__omit_mode_without_suffix(tmp_path: Path) -> None:
    cfg = SimpleNamespace(
        output_dir=tmp_path,
        default_date_prefix="19991231",
        output_stamp_mode="omit",
    )

    output_path = default_output_path(tmp_path / "input.csv", cfg)

    assert output_path == tmp_path / "output.input.csv"


@pytest.mark.unit
def test_default_output_path__require_mode_needs_date(tmp_path: Path) -> None:
    cfg = SimpleNamespace(
        output_dir=tmp_path,
        default_date_prefix=None,
        output_stamp_mode="require",
    )

    with pytest.raises(ValueError, match="date must be provided"):
        default_output_path(tmp_path / "input.csv", cfg)


@pytest.mark.unit
def test_default_output_path__strips_redundant_output_prefix(tmp_path: Path) -> None:
    cfg = SimpleNamespace(
        output_dir=tmp_path,
        default_date_prefix=None,
        output_stamp_mode="require",
    )

    output_path = default_output_path(
        tmp_path / "output.documents_20240101.csv",
        cfg,
        date="20251010",
    )

    assert output_path == tmp_path / "output.documents_20240101_20251010.csv"


@pytest.mark.unit
def test_default_output_path__handles_hidden_temp_files(tmp_path: Path) -> None:
    cfg = SimpleNamespace(
        output_dir=tmp_path,
        default_date_prefix=None,
        output_stamp_mode="omit",
    )

    output_path = default_output_path(
        tmp_path / ".output.assay_20240101.csv.tmp",
        cfg,
    )

    assert output_path == tmp_path / "output.assay_20240101.csv"


@pytest.mark.unit
def test_default_output_path__deduplicates_intermediate_suffix(tmp_path: Path) -> None:
    cfg = SimpleNamespace(
        output_dir=tmp_path,
        default_date_prefix=None,
        output_stamp_mode="omit",
    )

    output_path = default_output_path(
        tmp_path / "output.assay_20240101.csv.tmp",
        cfg,
    )

    assert output_path == tmp_path / "output.assay_20240101.csv"


@pytest.mark.unit
def test_default_output_path__normalises_chained_hidden_prefix(tmp_path: Path) -> None:
    cfg = SimpleNamespace(
        output_dir=tmp_path,
        default_date_prefix=None,
        output_stamp_mode="omit",
    )

    output_path = default_output_path(
        tmp_path / "output..output.activities_20240101.csv.tmp",
        cfg,
    )

    assert output_path == tmp_path / "output.activities_20240101.csv"


@pytest.mark.unit
def test_prepare_io_paths__strips_csv_suffix_from_output_stem(tmp_path: Path) -> None:
    args = SimpleNamespace(
        base_path=tmp_path,
        input_dir=None,
        output_dir=tmp_path,
        cache_dir=None,
        input_csv=tmp_path / "activities.csv",
        final_out=None,
        output_csv=None,
        raw_format=None,
        output_stamp_mode="date",
        raw_out=None,
        checkpoint=None,
        date="20251011",
        force=False,
        skip_existing=False,
    )

    parser_module.prepare_io_paths(
        args,
        output_stem="output.testitems_20240101.csv",
        suffix=".csv",
    )

    expected = (tmp_path / "output.testitems_20240101_20251011.csv").resolve()
    assert args.final_out == expected
    assert args.output_csv == expected
    assert expected.name.count(".csv") == 1
    assert "..output" not in expected.name


@pytest.mark.unit
def test_prepare_io_paths__normalises_prefixed_stem_without_date(tmp_path: Path) -> None:
    args = SimpleNamespace(
        base_path=tmp_path,
        input_dir=None,
        output_dir=tmp_path,
        cache_dir=None,
        input_csv=tmp_path / "activities.csv",
        final_out=None,
        output_csv=None,
        raw_format=None,
        output_stamp_mode=None,
        raw_out=None,
        checkpoint=None,
        date=None,
        force=False,
        skip_existing=False,
    )

    parser_module.prepare_io_paths(
        args,
        output_stem="output..output.testitems_20240101.csv",
        suffix=".csv",
    )

    expected = (tmp_path / "output.testitems_20240101.csv").resolve()
    assert args.final_out == expected
    assert args.output_csv == expected
    assert expected.name.count(".csv") == 1
    assert "..output" not in expected.name
