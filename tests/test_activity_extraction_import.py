"""Smoke tests for :mod:`activity_extraction_main`."""

from importlib import import_module


def test_activity_extraction_main_exposes_main() -> None:
    module = import_module("activity_extraction_main")
    assert hasattr(module, "main")


def test_activity_extraction_main_print_config() -> None:
    module = import_module("activity_extraction_main")
    # ``--print-config`` ensures the pipeline stops before network or IO work.
    exit_code = module.main(["--print-config"])
    assert exit_code == 0
