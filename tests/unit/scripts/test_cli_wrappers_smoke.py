"""Smoke tests ensuring CLI wrapper modules expose the expected API."""

from __future__ import annotations

import argparse
import inspect
from importlib import import_module

import pytest

_MAIN_MODULES = (
    "scripts.get_data",
    "scripts.get_activity_data",
    "scripts.get_assay_data",
    "scripts.get_document_data",
    "scripts.get_target_data",
    "scripts.get_testitem_data",
)

_BUILD_PARSER_MODULES = (
    "scripts.get_activity_data",
    "scripts.get_assay_data",
    "scripts.get_document_data",
    "scripts.get_target_data",
    "scripts.get_testitem_data",
)


@pytest.mark.unit
@pytest.mark.parametrize("module_name", _MAIN_MODULES)
def test_cli_wrapper_main__signature(module_name: str) -> None:
    module = import_module(module_name)
    main = getattr(module, "main")
    signature = inspect.signature(main)

    argv_parameter = signature.parameters.get("argv")
    assert argv_parameter is not None
    assert argv_parameter.default is None

    return_annotation = signature.return_annotation
    assert str(return_annotation) in {"int", "<class 'int'>"}


@pytest.mark.unit
@pytest.mark.parametrize("module_name", _BUILD_PARSER_MODULES)
def test_cli_wrapper_build_parser__returns_parser(module_name: str) -> None:
    module = import_module(module_name)
    build_parser = getattr(module, "build_parser")

    parser, *_ = build_parser()

    assert isinstance(parser, argparse.ArgumentParser)
