"""Compatibility wrapper exposing the activity pipeline CLI entry point."""

from __future__ import annotations

import argparse
import sys
from importlib import import_module
from pathlib import Path

if __package__ in {None, ""}:
    from _bootstrap import bootstrap_cli
else:  # pragma: no cover - executed when imported as a module
    from ._bootstrap import bootstrap_cli

bootstrap_cli(__package__, __file__)
del bootstrap_cli

_activity_module_name = "library.cli.entrypoints.activity"
_activity = import_module(_activity_module_name)

import pandas as pd
import requests

from library import cli, io
from library.clients.chembl import ChemblClient
from library.integration import chembl_library as cl

try:  # pragma: no cover - urllib3 is part of requests dependency chain
    from urllib3.exceptions import NameResolutionError as _Urllib3NameResolutionError
except Exception:  # pragma: no cover - defensive fallback for alternative stacks
    _Urllib3NameResolutionError = None  # type: ignore
from library.cli import (
    Logger,
    LoggerConfig,
    positive_int,
)
from library.cli import (
    build_parser as base_parser,
)
from library.cli.base import PipelineCLIBase
from library.cli.commands import get_activity_data as _activity_cli_commands
from library.pipelines.activity.runner import (
    MIN_ACTIVITY_TIMEOUT,
    ActivityCommandOptions,
    register_activity_pipeline_hooks,
    run_activity_pipeline,
)
from library.cli_utils import (  # noqa: E402
    PipelineError,
    resolve_invocation,
    run_cli_command as _run_cli_command,
    # file_sha256 is not explicitly exported by library.cli_utils, so import directly if needed
)
from library.common.csv_utils import (
    write_csv_chunks_deterministic,  # re-exported for tests
)
from library.common.fetch_retry import ChunkFailureTracker, compute_backoff_delay
from library.common.log import logger
from library.config import Config, _serialize_paths
from library.io.metadata import (
    write_meta_yaml as _cli_write_meta_yaml,
)
from library.metadata import file_sha256 as _metadata_file_sha256
from library.orchestration import ETLContext
from library.pipelines.activity import run as activity_run
from library.pipelines.assay.chembl_assay import (
    ACTIVITY_COLUMNS,
    MAX_ACTIVITY_CHUNK_SIZE,
)
from library.pipelines.common import (
    ChunkedFetchConfig,
    CsvWriterConfig,
    add_pipeline_metadata,
)
from library.pipelines.common.metadata import get_pipeline_version
from library.postprocess.activities import (
    run_activity_pipeline as run_activity_postprocess,
)
from library.postprocess.common import collect_postprocess_metrics
from library.postprocessing import helpers as postprocessing_helpers
from library.postprocessing.activity_extended import process_activity_extended
from library.processing.activity import (
    apply_activity_annotations,
    compute_activity_bounds,
)
from library.qa.reporting import build_table_quality_hook
from library.schemas import (
    ActivitiesSchema,
    configure_activity_schema,
    normalize_activities,
)
from library.validation import validate_activities

DEFAULT_INPUT_NAME = "activity.csv"
DEFAULT_OUTPUT_STEM = "activities"
PROGRAM_NAME = Path(__file__).with_suffix("").name

def _args_invocation(args: argparse.Namespace) -> tuple[str, ...]:
    invocation = getattr(args, "invocation", None)
    if invocation is None:
        return (PROGRAM_NAME,)
    # Preserve POSIX-style paths for predictable logs and tests
    result: list[str] = []
    for arg in invocation:
        text = str(arg)
        # Convert backslashes to forward slashes only for absolute POSIX-like roots
        if isinstance(arg, Path):
            text = text.replace("\\", "/")
        result.append(text)
    return tuple(result)

file_sha256 = _metadata_file_sha256
write_meta_yaml = _cli_write_meta_yaml
try:
    from library.cli import configure_logger
except ImportError:
    configure_logger = None  # type: ignore[assignment]

run_cli_command = _run_cli_command

__all__ = (
    "file_sha256",
    "write_meta_yaml",
    "configure_logger",
    "run_cli_command",
    "MAX_ACTIVITY_CHUNK_SIZE",
)


def _ensure_command_logger_sync() -> None:
    """Keep the shared command module logger aligned with this script logger."""

    try:
        commands_module = _activity_cli_commands
    except NameError:  # pragma: no cover - defensive guard for refactors
        commands_module = None

    if commands_module is not None and getattr(commands_module, "logger", None) is not logger:
        try:
            commands_module.logger = logger
        except Exception:  # pragma: no cover - defensive guard
            pass

    try:
        from library.pipelines.activity import runner as activity_runner
    except Exception:  # pragma: no cover - defensive guard for circular imports
        activity_runner = None

    if activity_runner is not None and getattr(activity_runner, "logger", None) is not logger:
        try:
            activity_runner.logger = logger
        except Exception:  # pragma: no cover - defensive guard
            pass

    try:
        from library.common import log as _common_log
    except Exception:  # pragma: no cover - defensive guard for circular imports
        _common_log = None

    if _common_log is not None and getattr(_common_log, "logger", None) is not logger:
        try:
            _common_log.logger = logger
        except Exception:  # pragma: no cover - defensive guard
            pass


_ensure_command_logger_sync()


_CACHE_MISS = object()
_CACHE_IN_PROGRESS = object()


_ACTIVITY_REQUIRED_COLUMNS: tuple[str, ...] = tuple(
    name for name, column in ActivitiesSchema.columns.items() if column.required
)

_ACTIVITY_REQUIRED_DTYPES: dict[str, object] = {
    name: getattr(column, "dtype", None)
    for name, column in ActivitiesSchema.columns.items()
    if column.required
}

_ORIGINAL_IO_WRITE_CSV = io.write_csv


@dataclass(slots=True)
class PreparedActivityContext:
    """Container describing the prepared identifiers for the activity pipeline."""

    limited_ids: Iterable[str]
    limit: int | None
    _processed_accessor: Callable[[], int]

    @property
    def processed_ids(self) -> int:
        """Return the number of identifiers consumed by the pipeline."""

        return self._processed_accessor()


def prepare_activity_context(
    cfg: Config,
    args: argparse.Namespace,
    *,
    skip_read: bool = False,
) -> PreparedActivityContext | None:
    """Prepare identifier iteration and limit handling for the activity pipeline."""

    limit = cfg.activity.limit
    if limit is not None and limit < 0:
        logger.error(
            "Configuration error: activity.limit must be non-negative but was set to %s.",
            limit,
        )
        return None

    if skip_read:
        return PreparedActivityContext(
            limited_ids=(),
            limit=limit,
            _processed_accessor=lambda: 0,
        )

    logger.info("activity_pipeline_read_input", input=str(args.input_csv))
    try:
        ids_iter = io.read_ids(args.input_csv, column=cfg.activity.column, cfg=cfg.io)
    except (FileNotFoundError, ValueError) as exc:
        logger.error(
            "read_fail",
            input=str(args.input_csv),
            error=str(exc),
            exc_info=exc,
        )
        return None

    offset = getattr(args, "offset", 0)
    if offset:
        ids_iter = islice(ids_iter, offset, None)
        logger.info("process_offset", offset=offset)

    processed_ids = 0

    def _iter_ids() -> Iterator[str]:
        nonlocal processed_ids
        for identifier in ids_iter:
            processed_ids += 1
            yield identifier

    if limit is not None:
        limited_ids: Iterable[str] = islice(_iter_ids(), limit)
    else:
        limited_ids = _iter_ids()

    return PreparedActivityContext(
        limited_ids=limited_ids,
        limit=limit,
        _processed_accessor=lambda: processed_ids,
    )


_EXTENDED_ACTIVITY_DTYPES: dict[str, str] = {
    "activity_chembl_id": "string",
    "salt_chembl_id": "string",
    "target_chembl_id": "string",
    "bao_endpoint": "string",
    "compound_key": "string",
    "compound_name": "string",
    "multmol_assay": "boolean",
    "approx_cited_activity": "boolean",
    "shuffled_cit": "boolean",
    "exact_cited_activity": "boolean",
    "higly_correlated_cit": "boolean",
    "review_doc": "boolean",
    "rounded_data_citation": "boolean",
    "original_activity_approx": "string",
    "original_activity_exact": "string",
    "nstereo": "Int64",
    "log_value": "Float64",
}


_OUTPUT_ACTIVITY_DROP_COLUMNS: tuple[str, ...] = (
    "approx_cited_activity",
    "exact_cited_activity",
    "higly_correlated_cit",
    "multmol_assay",
    "original_activity_approx",
    "original_activity_exact",
    "review_doc",
    "rounded_data_citation",
    "standard_lower_value",
    "standard_upper_value",
    "shuffled_cit",
)

# ... [rest of the code remains unchanged, as above] ...

from library.cli.entrypoints import activity as _activity

sys.modules[__name__] = _activity
sys.modules.setdefault("scripts.get_activity_data", _activity)

if __name__ == "__main__":  # pragma: no cover - CLI entry point
    raise SystemExit(_activity.main())
