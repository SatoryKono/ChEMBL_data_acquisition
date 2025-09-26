"""Smoke tests for :mod:`scripts.get_activity_data`."""

from importlib import import_module


def test_get_activity_data_cli_exposes_main() -> None:
    module = import_module("scripts.get_activity_data")
    assert hasattr(module, "main")


def test_get_activity_data_cli_print_config() -> None:
    module = import_module("scripts.get_activity_data")
    # ``--print-config`` ensures the pipeline stops before network or IO work.
    exit_code = module.main(["--print-config"])
    assert exit_code == 0
