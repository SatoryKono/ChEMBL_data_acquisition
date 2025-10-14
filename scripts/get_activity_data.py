"""Compatibility wrapper exposing the activity pipeline CLI entry point."""

from __future__ import annotations

import sys
from collections.abc import Iterable

import argparse
from dataclasses import dataclass
from datetime import datetime
from importlib import import_module
from pathlib import Path
from types import ModuleType
from typing import TYPE_CHECKING, Sequence

# ruff: noqa: E402  # bootstrap alters import order for script compatibility
if TYPE_CHECKING:
    from . import _bootstrap as _bootstrap_module
elif __package__ in {None, ""}:
    import _bootstrap as _bootstrap_module  # pragma: no cover - CLI fallback
else:  # pragma: no cover - executed when imported as a module
    from . import _bootstrap as _bootstrap_module

bootstrap_cli = _bootstrap_module.bootstrap_cli

bootstrap_cli(__package__, __file__)
del bootstrap_cli
del _bootstrap_module


def _export_module_api(
    module: ModuleType, *, extra: Iterable[str] = ()
) -> tuple[str, ...]:
    """Expose ``module`` attributes in the wrapper namespace."""

    exported: dict[str, object] = {}
    for name, value in module.__dict__.items():
        if name.startswith("__"):
            continue
        exported[name] = value

    for name in extra:
        if hasattr(module, name):
            exported[name] = getattr(module, name)

    globals().update(exported)

    module_all = tuple(getattr(module, "__all__", ()))
    if module_all:
        ordered = list(module_all)
        for name in exported:
            if name not in ordered:
                ordered.append(name)
    else:
        ordered = list(exported)
    return tuple(ordered)


def _augment_activity_module(module: ModuleType) -> tuple[str, ...]:
    """Install backwards compatible exports on the activity CLI module."""

    extra_exports: dict[str, object] = {}

    try:  # pragma: no cover - defensive guard for optional dependencies
        from library.pipelines.activity import runner as activity_runner
    except Exception:  # pragma: no cover - ensure import failures do not break CLI
        activity_runner = None
    else:
        extra_exports.update(
            {
                "MAX_ACTIVITY_CHUNK_SIZE": activity_runner.MAX_ACTIVITY_CHUNK_SIZE,
                "register_activity_pipeline_hooks": activity_runner.register_activity_pipeline_hooks,
            }
        )

    try:  # pragma: no cover - defensive guard for optional dependencies
        from library.cli.entrypoints import activity as activity_entrypoint
    except (
        Exception
    ):  # pragma: no cover - the entrypoint may not be importable during tests
        activity_entrypoint = None
    else:
        emit_completion = getattr(activity_entrypoint, "_emit_completion_message", None)
        if emit_completion is not None:
            extra_exports.setdefault("_emit_completion_message", emit_completion)

    for name, value in extra_exports.items():
        setattr(module, name, value)

    existing_all = tuple(getattr(module, "__all__", ()))
    missing = tuple(name for name in extra_exports if name not in existing_all)
    if missing:
        module.__all__ = existing_all + missing  # type: ignore[attr-defined]
    return tuple(extra_exports)


def _load_activity_module() -> ModuleType:
    module = import_module("library.cli.commands.get_activity_data")
    return module


_MODULE = _load_activity_module()
_EXTRA_EXPORTS = _augment_activity_module(_MODULE)
_MODULE._synchronise_wrapper_module()
__all__ = _export_module_api(_MODULE, extra=_EXTRA_EXPORTS)


from library.io.config_loader import load_config
from library.io.metadata_writer import save_metadata
from library.io.output_writer import save_standard_outputs
from library.postprocessing.activity.steps import (
    ACTIVITY_COLUMN_ORDER,
    fetch_normalize_activity,
    generate_activity_reports,
)
from library.utils.logging import get_logger

TABLE_NAME = "activity"


def _positive_int(value: str) -> int:
    parsed = int(value)
    if parsed <= 0:
        raise argparse.ArgumentTypeError("limit must be positive")
    return parsed


def _date_tag(value: str) -> str:
    try:
        datetime.strptime(value, "%Y%m%d")
    except ValueError as exc:  # noqa: TRY003 - convert into argparse error
        raise argparse.ArgumentTypeError("date tag must be formatted as YYYYMMDD") from exc
    return value


@dataclass(frozen=True)
class PipelineArgs:
    """Normalised CLI arguments for the simplified activity pipeline."""

    limit: int
    date_tag: str
    output_dir: Path
    config: Path | None


def parse_activity_args(
    argv: Sequence[str] | None = None,
) -> tuple[argparse.Namespace, PipelineArgs]:
    """Parse simplified CLI arguments returning namespace and dataclass view."""

    parser = argparse.ArgumentParser(
        description="Fetch and normalise ChEMBL activity data",
        prog="activity-normalize",
    )
    parser.add_argument(
        "--limit",
        type=_positive_int,
        default=1000,
        help="Number of activity rows to fetch",
    )
    parser.add_argument(
        "--date-tag",
        type=_date_tag,
        default=datetime.utcnow().strftime("%Y%m%d"),
        help="Execution date token in YYYYMMDD format",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("data/output"),
        help="Directory for generated artefacts",
    )
    parser.add_argument(
        "--config",
        type=Path,
        default=None,
        help="Optional configuration path (defaults to config/config.yaml)",
    )

    namespace = parser.parse_args(argv)
    args = PipelineArgs(
        limit=namespace.limit,
        date_tag=namespace.date_tag,
        output_dir=Path(namespace.output_dir),
        config=Path(namespace.config) if namespace.config is not None else None,
    )
    return namespace, args


def run_activity_cli(argv: Sequence[str] | None = None) -> int:
    """Execute the simplified activity pipeline returning an exit status."""

    namespace, args = parse_activity_args(argv)
    logger = get_logger(__name__).bind(cli_variant="activity_simplified")

    module_load_config = getattr(_MODULE, "load_config", load_config)
    try:
        config = (
            module_load_config(args.config)
            if args.config is not None
            else module_load_config(None)
        )
    except Exception as exc:  # noqa: BLE001 - surface configuration failures
        logger.error("activity_cli_config_error", error=str(exc))
        return 1

    args.output_dir.mkdir(parents=True, exist_ok=True)

    try:
        module_fetch = getattr(_MODULE, "fetch_normalize_activity", fetch_normalize_activity)
        module_generate = getattr(
            _MODULE, "generate_activity_reports", generate_activity_reports
        )
        module_save_outputs = getattr(
            _MODULE, "save_standard_outputs", save_standard_outputs
        )
        module_save_metadata = getattr(_MODULE, "save_metadata", save_metadata)
        column_order = getattr(_MODULE, "ACTIVITY_COLUMN_ORDER", ACTIVITY_COLUMN_ORDER)

        dataset = module_fetch(args.limit)
        activity_data = module_generate(dataset)
        artifacts = module_save_outputs(
            activity_data.dataset,
            activity_data.correlation_report,
            activity_data.quality_report,
            TABLE_NAME,
            args.date_tag,
            output_dir=args.output_dir,
            column_order=column_order,
            key_columns=["activity_id"],
        )
        artifact_paths = [
            artifacts.dataset,
            artifacts.correlation_report,
            artifacts.quality_report,
        ]
        module_save_metadata(
            TABLE_NAME,
            args.date_tag,
            namespace,
            qc_summary=activity_data.qc_summary,
            output_dir=args.output_dir,
            artifacts=artifact_paths,
            sources=[config.path.name],
        )
    except Exception as exc:  # noqa: BLE001 - top-level failure boundary
        logger.exception("activity_cli_pipeline_failed", error=str(exc))
        return 1

    logger.info("activity_cli_pipeline_success", rows=int(activity_data.dataset.shape[0]))
    return 0


for _export_name, _export_value in {
    "PipelineArgs": PipelineArgs,
    "parse_activity_args": parse_activity_args,
    "run_activity_cli": run_activity_cli,
    "fetch_normalize_activity": fetch_normalize_activity,
    "generate_activity_reports": generate_activity_reports,
    "save_standard_outputs": save_standard_outputs,
    "save_metadata": save_metadata,
    "load_config": load_config,
    "ACTIVITY_COLUMN_ORDER": ACTIVITY_COLUMN_ORDER,
}.items():
    setattr(_MODULE, _export_name, _export_value)


__all__ = tuple(
    dict.fromkeys(
        (*__all__, "PipelineArgs", "parse_activity_args", "run_activity_cli")
    )
)


def __getattr__(name: str) -> object:
    """Proxy attribute access to :mod:`library.cli.commands.get_activity_data`."""

    try:
        return getattr(_MODULE, name)
    except AttributeError as exc:  # pragma: no cover - passthrough for missing attrs
        raise AttributeError(name) from exc


def __dir__() -> list[str]:  # pragma: no cover - passthrough helper
    return sorted({*globals().keys(), *dir(_MODULE)})


class _Adapter(ModuleType):
    """Proxy module syncing attribute writes back to ``_MODULE``."""

    def __getattr__(self, name: str) -> object:  # pragma: no cover - passthrough helper
        return __getattr__(name)

    def __setattr__(
        self, name: str, value: object
    ) -> None:  # pragma: no cover - passthrough helper
        setattr(_MODULE, name, value)
        super().__setattr__(name, value)

    def __dir__(self) -> list[str]:  # pragma: no cover - passthrough helper
        return __dir__()


module_proxy = sys.modules[__name__]
if module_proxy is not _MODULE:
    module_proxy.__class__ = _Adapter


if __name__ == "__main__":  # pragma: no cover - CLI entry point
    raise SystemExit(_MODULE.main())
