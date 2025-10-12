from __future__ import annotations

import argparse
from collections.abc import Callable
from pathlib import Path
from types import SimpleNamespace

import pytest

import library.cli as cli_module
import library.cli_utils as cli_utils_module
from library.utils.cli_tools import check_determinism, chunk_io_main, mapper_main


@pytest.mark.unit
@pytest.mark.parametrize(
    ("module", "argv_factory"),
    [
        (
            chunk_io_main,
            lambda tmp_path: [
                "--input",
                str(_ensure_file(tmp_path, "chunk_io_input.csv", "value\n")),
                "--final-out",
                str(tmp_path / "chunk_io_output.csv"),
                "--checkpoint",
                str(tmp_path / "chunk_io_state.json"),
                "--chunk-size",
                "2",
            ],
        ),
        (
            mapper_main,
            lambda tmp_path: [
                "--input",
                str(_ensure_file(tmp_path, "mapper_input.csv", "chembl_id\nCHEMBL1\n")),
                "--final-out",
                str(tmp_path / "mapper_output.csv"),
            ],
        ),
        (check_determinism, lambda _: []),
    ],
    ids=["chunk_io", "mapper", "check_determinism"],
)
def test_cli_wrapper__run_id_and_generated_at_stable(
    module: object,
    argv_factory: Callable[[Path], list[str]],
    tmp_path_factory: pytest.TempPathFactory,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    base_path = tmp_path_factory.mktemp("cli-wrapper")

    captured_run_ids: list[str | None] = []
    captured_configs: list[SimpleNamespace] = []

    def _dummy_logger() -> SimpleNamespace:
        return SimpleNamespace(
            info=lambda *_, **__: None,
            error=lambda *_, **__: None,
            warning=lambda *_, **__: None,
            debug=lambda *_, **__: None,
        )

    def _configure_logger(cfg: cli_module.LoggerConfig) -> SimpleNamespace:
        snapshot = SimpleNamespace(run_id=cfg.run_id, generated_at=cfg.generated_at)
        captured_configs.append(snapshot)
        return _dummy_logger()

    monkeypatch.setattr(cli_module, "configure_logger", _configure_logger)
    monkeypatch.setattr(cli_utils_module.cli, "configure_logger", _configure_logger)

    cache_root = base_path / "cache"
    output_root = base_path / "output"
    cache_root.mkdir(parents=True, exist_ok=True)
    output_root.mkdir(parents=True, exist_ok=True)

    def _fake_apply_config(
        args: argparse.Namespace,
        parser: argparse.ArgumentParser,
        config_path: object,
        mapping: dict[str, str] | None = None,
        *,
        base_parser: argparse.ArgumentParser | None = None,
    ) -> SimpleNamespace:
        del parser, config_path, mapping, base_parser
        args._config_metadata = SimpleNamespace(snapshot={})
        cfg = SimpleNamespace(
            io=SimpleNamespace(
                exist_ok=True,
                output_dir=output_root,
                cache_dir=cache_root,
            )
        )
        return cfg

    monkeypatch.setattr(cli_utils_module, "apply_config_overrides", _fake_apply_config)

    def _fake_run(cfg: object, args: argparse.Namespace) -> int:
        del cfg
        captured_run_ids.append(getattr(args, "run_id", None))
        return 0

    monkeypatch.setattr(module, "run", _fake_run)

    def _invoke_once() -> tuple[str | None, SimpleNamespace]:
        start_config_len = len(captured_configs)
        start_run_len = len(captured_run_ids)
        argv = list(argv_factory(base_path))
        exit_code = module.main(argv)
        assert exit_code == 0
        assert len(captured_run_ids) == start_run_len + 1
        run_id_value = captured_run_ids[-1]
        snapshots = captured_configs[start_config_len:]
        assert snapshots, "configure_logger must be called"
        return run_id_value, snapshots[-1]

    first_run_id, first_cfg = _invoke_once()
    second_run_id, second_cfg = _invoke_once()

    assert first_run_id == second_run_id
    assert first_cfg.generated_at == second_cfg.generated_at


def _ensure_file(base: Path, name: str, content: str) -> Path:
    target = base / name
    target.parent.mkdir(parents=True, exist_ok=True)
    if not target.exists():
        target.write_text(content, encoding="utf-8")
    return target
