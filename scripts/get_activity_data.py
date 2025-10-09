"""Thin wrapper exposing the activity pipeline CLI entry point."""

from __future__ import annotations

from functools import wraps

if __package__ in {None, ""}:
    from _bootstrap import bootstrap_cli
else:  # pragma: no cover - executed when imported as a module
    from ._bootstrap import bootstrap_cli

bootstrap_cli(__package__, __file__)
del bootstrap_cli

from library.cli.entrypoints import activity as _activity
from library.pipelines.activity import runner as _activity_runner
from library.pipelines.activity.runner import register_activity_pipeline_hooks
from library.pipelines.assay.chembl_assay import MAX_ACTIVITY_CHUNK_SIZE

ActivityCommandOptions = _activity.ActivityCommandOptions
ActivityPipelineCLI = _activity.ActivityPipelineCLI
PreparedActivityContext = _activity.PreparedActivityContext

DEFAULT_INPUT_NAME = _activity.DEFAULT_INPUT_NAME
DEFAULT_OUTPUT_STEM = _activity.DEFAULT_OUTPUT_STEM
PROGRAM_NAME = _activity.PROGRAM_NAME
MIN_ACTIVITY_TIMEOUT = _activity.MIN_ACTIVITY_TIMEOUT

cl = _activity.cl
configure_logger = _activity.configure_logger
file_sha256 = _activity.file_sha256
io = _activity.io
logger = _activity.logger
run_cli_command = _activity.run_cli_command
sleep = _activity.sleep
write_csv_chunks_deterministic = _activity.write_csv_chunks_deterministic
write_meta_yaml = _activity.write_meta_yaml

_args_invocation = _activity._args_invocation


def _sync_runtime_dependencies() -> None:
    module_globals = globals()
    logger_obj = module_globals["logger"]
    _activity.logger = logger_obj
    _activity_runner.logger = logger_obj
    _activity.sleep = module_globals["sleep"]
    _activity.write_csv_chunks_deterministic = module_globals[
        "write_csv_chunks_deterministic"
    ]
    _activity.run_chembl = module_globals["run_chembl"]
    _activity._emit_completion_message = module_globals["_emit_completion_message"]


def _with_runtime_sync(func):
    @wraps(func)
    def wrapper(*args, **kwargs):
        _sync_runtime_dependencies()
        return func(*args, **kwargs)

    return wrapper


_emit_completion_message = _with_runtime_sync(_activity._emit_completion_message)
_ensure_extended_activity_columns = _with_runtime_sync(
    _activity._ensure_extended_activity_columns
)
_ensure_molecule_pref_name = _with_runtime_sync(_activity._ensure_molecule_pref_name)
_ensure_src_assay_id = _with_runtime_sync(_activity._ensure_src_assay_id)
_filter_activity_output_columns = _with_runtime_sync(
    _activity._filter_activity_output_columns
)
_load_assay_src_lookup = _with_runtime_sync(_activity._load_assay_src_lookup)
main = _with_runtime_sync(_activity.main)
prepare_activity_context = _with_runtime_sync(_activity.prepare_activity_context)
run = _with_runtime_sync(_activity.run)
run_chembl = _with_runtime_sync(_activity.run_chembl)

register_activity_pipeline_hooks(
    runner=run_chembl,
    emit_completion_message=_emit_completion_message,
)

__all__ = (
    "ActivityCommandOptions",
    "ActivityPipelineCLI",
    "DEFAULT_INPUT_NAME",
    "DEFAULT_OUTPUT_STEM",
    "MAX_ACTIVITY_CHUNK_SIZE",
    "MIN_ACTIVITY_TIMEOUT",
    "PROGRAM_NAME",
    "PreparedActivityContext",
    "cl",
    "configure_logger",
    "file_sha256",
    "io",
    "logger",
    "main",
    "prepare_activity_context",
    "run",
    "run_chembl",
    "run_cli_command",
    "sleep",
    "write_csv_chunks_deterministic",
    "write_meta_yaml",
    "_args_invocation",
    "_emit_completion_message",
    "_ensure_extended_activity_columns",
    "_ensure_molecule_pref_name",
    "_ensure_src_assay_id",
    "_filter_activity_output_columns",
    "_load_assay_src_lookup",
)


def __getattr__(name: str) -> object:
    """Proxy attribute access to :mod:`library.cli.entrypoints.activity`."""

    try:
        return getattr(_activity, name)
    except AttributeError as exc:  # pragma: no cover - passthrough for missing attrs
        raise AttributeError(name) from exc


if __name__ == "__main__":  # pragma: no cover - integration execution path
    raise SystemExit(main())
