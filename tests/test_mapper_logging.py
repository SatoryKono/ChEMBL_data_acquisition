import importlib
import logging
import sys


def test_mapper_library_has_no_logging_side_effect() -> None:
    root = logging.getLogger()
    for handler in list(root.handlers):
        root.removeHandler(handler)
    assert not root.handlers
    module_name = "library.mapper_library"
    if module_name in sys.modules:
        del sys.modules[module_name]
    importlib.import_module(module_name)
    assert not root.handlers
