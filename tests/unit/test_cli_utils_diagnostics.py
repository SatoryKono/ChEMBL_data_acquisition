"""Unit tests for diagnostic toggles in :mod:`library.cli.utils`."""

from __future__ import annotations

import argparse
from types import SimpleNamespace

import pytest

from library.cli import LoggerConfig
from library.cli.utils import run_cli_command


class _LoggerStub:
    def __init__(self) -> None:
        self.calls: list[tuple[str, tuple[object, ...], dict[str, object]]] = []

    def info(self, *args: object, **kwargs: object) -> None:  # pragma: no cover - helper
        self.calls.append(("info", args, kwargs))

    def error(self, *args: object, **kwargs: object) -> None:  # pragma: no cover - helper
        self.calls.append(("error", args, kwargs))

    def warning(self, *args: object, **kwargs: object) -> None:  # pragma: no cover - helper
        self.calls.append(("warning", args, kwargs))

    def debug(self, *args: object, **kwargs: object) -> None:  # pragma: no cover - helper
        self.calls.append(("debug", args, kwargs))


@pytest.mark.unit
def test_run_cli_command__disables_doc_quality_without_diagnostics(
    tmp_path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Default execution should disable QC profiling when diagnostics stay off."""

    parser = argparse.ArgumentParser(prog="activity")
    config_path = tmp_path / "config.yaml"
    config_path.write_text("{}", encoding="utf-8")
    parser.add_argument("--config", default=str(config_path))

    logger_stub = _LoggerStub()

    def _fake_configure_logger(cfg: LoggerConfig) -> _LoggerStub:
        return logger_stub

    cfg = SimpleNamespace(
        system=SimpleNamespace(doc_quality=SimpleNamespace(enable=True)),
    )

    def _fake_apply_config_overrides(*_args, **_kwargs):
        return cfg

    monkeypatch.setattr("library.cli.utils.configure_logger", _fake_configure_logger)
    monkeypatch.setattr("library.cli.utils.prepare_io_paths", lambda _args: None)
    monkeypatch.setattr("library.cli.utils.apply_config_overrides", _fake_apply_config_overrides)
    monkeypatch.setattr("library.cli.utils.ensure_dirs", lambda _cfg: None)

    captured: dict[str, object] = {}

    def _fake_run(cfg_arg, args_arg):
        captured["doc_quality"] = cfg_arg.system.doc_quality.enable
        return 0

    args = argparse.Namespace(
        config=str(config_path),
        log_level="INFO",
        verbose=False,
        run_id=None,
        emit_legacy_artifacts=False,
        debug=False,
        keep_intermediate=False,
        postprocess=False,
        print_config=False,
    )

    log_cfg = LoggerConfig(level="INFO", run_id="test-run")

    exit_code = run_cli_command(
        args=args,
        parser=parser,
        log_cfg=log_cfg,
        mapping={},
        run=_fake_run,
        logger=logger_stub,
    )

    assert exit_code == 0
    assert captured["doc_quality"] is False
