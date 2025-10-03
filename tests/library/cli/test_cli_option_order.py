"""Tests for global options specified before sub-commands."""

from __future__ import annotations

import argparse
from pathlib import Path

import pytest

import library.cli.utils as cli_utils
from library.config import Config
from scripts import get_document_data as gdd
from scripts import get_target_data as gtd


def test_target_input_before_subcommand(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, cfg: Config
) -> None:
    """Options provided before the sub-command should be honoured."""
    input_csv = tmp_path / "targets.csv"
    config_yaml = tmp_path / "config.yaml"
    config_yaml.write_text("dummy: true\n")

    captured: dict[str, Path] = {}

    def fake_apply(
        args: argparse.Namespace,
        parser: argparse.ArgumentParser,
        cfg_path: str | Path,
        mapping: dict[str, str] | None = None,
        *,
        base_parser: argparse.ArgumentParser | None = None,
    ) -> Config:
        return cfg

    def fake_run_all(cfg: Config, args: argparse.Namespace) -> int:
        captured["input_csv"] = Path(args.input_csv)
        return 0

    monkeypatch.setattr(cli_utils, "apply_config_overrides", fake_apply)
    monkeypatch.setattr(gtd, "run_all", fake_run_all)
    monkeypatch.setattr(cli_utils, "ensure_dirs", lambda cfg: None)

    rc = gtd.main(["--config", str(config_yaml), "--input", str(input_csv), "all"])
    assert rc == 0
    assert captured["input_csv"] == input_csv


def test_document_input_before_subcommand(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, cfg: Config
) -> None:
    """Global options before the sub-command work for document CLI."""
    input_csv = tmp_path / "docs.csv"
    config_yaml = tmp_path / "config.yaml"
    config_yaml.write_text("dummy: true\n")

    captured: dict[str, Path] = {}

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
        captured["input_csv"] = Path(args.input_csv)
        return 0

    monkeypatch.setattr(cli_utils, "apply_config_overrides", fake_apply)
    monkeypatch.setattr(gdd, "run", fake_run)
    monkeypatch.setattr(cli_utils, "ensure_dirs", lambda cfg: None)

    rc = gdd.main(["--config", str(config_yaml), "--input", str(input_csv), "chembl"])
    assert rc == 0
    assert captured["input_csv"] == input_csv
