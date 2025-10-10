import logging
import warnings
from unittest import mock

from library.common.logging_setup import LoggerConfig, configure_logger


def test_configure_logger__closes_existing_handlers() -> None:
    root_logger = logging.getLogger()
    warnings_logger = logging.getLogger("py.warnings")

    original_root_handlers = list(root_logger.handlers)
    original_warnings_handlers = list(warnings_logger.handlers)
    original_root_level = root_logger.level
    original_warnings_level = warnings_logger.level
    original_warnings_propagate = warnings_logger.propagate
    was_capturing_warnings = warnings.showwarning is getattr(
        logging, "_showwarning", object()
    )

    root_handler = mock.create_autospec(logging.FileHandler, instance=True)
    warnings_handler = mock.create_autospec(logging.FileHandler, instance=True)

    root_logger.handlers = [root_handler]
    warnings_logger.handlers = [warnings_handler]

    try:
        configure_logger(LoggerConfig())
    finally:
        root_logger.handlers = original_root_handlers
        root_logger.setLevel(original_root_level)
        warnings_logger.handlers = original_warnings_handlers
        warnings_logger.setLevel(original_warnings_level)
        warnings_logger.propagate = original_warnings_propagate
        logging.captureWarnings(was_capturing_warnings)

    root_handler.close.assert_called_once_with()
    warnings_handler.close.assert_called_once_with()
