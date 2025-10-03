"""Tests for keyboard layout aliases in :mod:`scripts.get_target_data`."""

from __future__ import annotations

import argparse
from pathlib import Path

import pytest

import library.cli.utils as cli_utils
from library.config import Config
from scripts import get_target_data as gtd


@pytest.mark.parametrize(
    "alias, expected",
    [
        ("фдд", "all"),
        ("СРУЬИД", "chembl"),
        ("гтшзкще", "uniprot"),
        ("Шгзрфк", "iuphar"),
    ],
)
def test_target_subcommand_keyboard_alias(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    cfg: Config,
    alias: str,
    expected: str,
) -> None:
    """Running the CLI with a Russian keyboard alias resolves to the command."""

    input_csv = tmp_path / "targets.csv"
    config_yaml = tmp_path / "config.yaml"
    config_yaml.write_text("dummy: true\n")

    resolved_command: dict[str, str] = {}

    def fake_apply(
        args: argparse.Namespace,
        parser: argparse.ArgumentParser,
        cfg_path: str | Path,
        mapping: dict[str, str] | None = None,
        *,
        base_parser: argparse.ArgumentParser | None = None,
    ) -> Config:
        return cfg

    def fake_run(cfg: Config, args: argparse.Namespace) -> int:
        resolved_command["command"] = args.command
        return 0

    monkeypatch.setattr(cli_utils, "apply_config_overrides", fake_apply)
    monkeypatch.setattr(cli_utils, "ensure_dirs", lambda cfg: None)
    monkeypatch.setattr(gtd, "run", fake_run)

    exit_code = gtd.main([
        "--config",
        str(config_yaml),
        "--input",
        str(input_csv),
        alias,
    ])

    assert exit_code == 0
    assert resolved_command["command"] == expected
