"""Stable console script entry points for data acquisition commands.

The goal of this module is to provide a single import location that can be
referenced from :mod:`pyproject.toml` so that all console scripts resolve to
functions shipped within the :mod:`library` package. Historically the project
exposed a mixture of ``scripts`` modules and thin ``library.cli.commands``
wrappers which made it difficult to discover the supported entry points. By
normalising everything here we keep compatibility with the existing command
implementations while surfacing them under a predictable namespace.

Each public function matches a console command and delegates to the
corresponding implementation from :mod:`library.cli.commands` or the legacy
``scripts``/``tools`` packages. The helpers lean on
:func:`library.cli.commands._run` to preserve the historical command-line
behaviour (including optional ``argv`` injection for tests).
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import Final

from library.cli import logging as _cli_logging
from library.cli.commands import _run as _run_command
from library.cli.commands import get_target_data as _target_command
from library.cli.entrypoints import (
    activity as _activity_entrypoint,
    assay as _assay_entrypoint,
    document as _document_entrypoint,
    get_data as _orchestrator_entrypoint,
    testitem as _testitem_entrypoint,
)

# ---------------------------------------------------------------------------
# Generic dispatch helper
# ---------------------------------------------------------------------------


def _dispatch(module: str, argv: Sequence[str] | None = None) -> int:
    """Import ``module`` and execute its ``main`` function.

    ``module`` may reference either a fully-qualified module path or one of the
    historical short names understood by :func:`library.cli.commands._run`.
    ``argv`` is forwarded to the resolved ``main`` callable when it accepts
    parameters so that tests can provide deterministic argument vectors.
    """

    return _run_command(module, argv)


# ---------------------------------------------------------------------------
# Public entry points (sorted alphabetically by command name)
# ---------------------------------------------------------------------------


def check_determinism_main(argv: Sequence[str] | None = None) -> int:
    """Entry point for the ``check-determinism`` CLI."""

    return _dispatch("library.utils.cli_tools.check_determinism", argv)


def chunk_io_main(argv: Sequence[str] | None = None) -> int:
    """Entry point for the ``chunk-io`` CLI."""

    return _dispatch("library.utils.cli_tools.chunk_io_main", argv)


def csv_utils_main(argv: Sequence[str] | None = None) -> int:
    """Entry point for the ``csv-utils`` CLI."""

    return _dispatch("library.utils.cli_tools.csv_utils_main", argv)


def get_activities_main(argv: Sequence[str] | None = None) -> int:
    """Entry point for the ``get-activities`` CLI."""

    return _dispatch("library.utils.cli_tools.get_activities", argv)


def get_activity_data_main(argv: Sequence[str] | None = None) -> int:
    """Entry point for the ``get-activity-data`` pipeline CLI."""

    raw_args: list[str]
    if argv is None:
        raw_args = []
    else:
        raw_args = [str(arg) for arg in argv]

    def _has_date_override(items: list[str]) -> bool:
        for item in items:
            if item == "--date" or item.startswith("--date="):
                return True
        return False

    if not _has_date_override(raw_args):
        raw_args.extend(["--date", _cli_logging._current_date_str()])

    return _activity_entrypoint.main(raw_args)


def get_assay_data_main(argv: Sequence[str] | None = None) -> int:
    """Entry point for the ``get-assay-data`` pipeline CLI."""

    return _assay_entrypoint.main(argv)


def get_data_main(argv: Sequence[str] | None = None) -> int:
    """Entry point for the orchestrating ``get-data`` CLI."""

    return _orchestrator_entrypoint.main(argv)


def get_document_data_main(argv: Sequence[str] | None = None) -> int:
    """Entry point for the ``get-document-data`` pipeline CLI."""

    return _document_entrypoint.main(argv)


def get_document_type_main(argv: Sequence[str] | None = None) -> int:
    """Entry point for the ``get-document-type`` helper CLI."""

    return _dispatch("library.utils.cli_tools.get_document_type", argv)


def get_input_initialisation_main(argv: Sequence[str] | None = None) -> int:
    """Entry point for the ``get-input-initialisation`` helper CLI."""

    return _dispatch("library.utils.cli_tools.get_input_initialisation", argv)


def get_target_data_main(argv: Sequence[str] | None = None) -> int:
    """Entry point for the ``get-target-data`` pipeline CLI."""

    return _target_command.main(argv)


def get_testitem_data_main(argv: Sequence[str] | None = None) -> int:
    """Entry point for the ``get-testitem-data`` pipeline CLI."""

    return _testitem_entrypoint.main(argv)


def make_md_summary_main(argv: Sequence[str] | None = None) -> int:
    """Entry point for the ``make-md-summary`` reporting CLI."""

    return _dispatch("tools.make_md_summary", argv)


def mapper_main(argv: Sequence[str] | None = None) -> int:
    """Entry point for the ``mapper`` helper CLI."""

    return _dispatch("library.utils.cli_tools.mapper_main", argv)


def table_quality_main(argv: Sequence[str] | None = None) -> int:
    """Entry point for the ``table-quality`` QA CLI."""

    return _dispatch("library.utils.cli_tools.table_quality_main", argv)


__all__: Final[list[str]] = [
    "check_determinism_main",
    "chunk_io_main",
    "csv_utils_main",
    "get_activities_main",
    "get_activity_data_main",
    "get_assay_data_main",
    "get_data_main",
    "get_document_data_main",
    "get_document_type_main",
    "get_input_initialisation_main",
    "get_target_data_main",
    "get_testitem_data_main",
    "make_md_summary_main",
    "mapper_main",
    "table_quality_main",
]
