"""Compatibility wrapper exposing :mod:`library.cli.commands.get_data`."""

from __future__ import annotations

import sys
from collections.abc import Callable, Iterable
from dataclasses import replace
from importlib import import_module
from types import ModuleType
from typing import TYPE_CHECKING, Any

# ruff: noqa: E402  # bootstrap alters import order for script compatibility
if TYPE_CHECKING:
    from . import _bootstrap as _bootstrap_module
elif __package__ in {None, ""}:
    import _bootstrap as _bootstrap_module  # pragma: no cover - CLI fallback
else:  # pragma: no cover - executed when imported as a module
    try:
        from . import _bootstrap as _bootstrap_module
    except ImportError:  # pragma: no cover - namespace fallback for tests
        _bootstrap_module = import_module("scripts._bootstrap")

bootstrap_cli = _bootstrap_module.bootstrap_cli

bootstrap_cli(__package__, __file__)
del bootstrap_cli
del _bootstrap_module


def _export_module_api(module: ModuleType, *, extra: Iterable[str] = ()) -> tuple[str, ...]:
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


def _load_module() -> ModuleType:
    try:
        return import_module("library.cli.commands.get_data")
    except ModuleNotFoundError as exc:  # pragma: no cover - import guard
        missing = exc.name or "dependency"
        message = (
            f"Missing required dependency '{missing}'.\n"
            "Install the project requirements before running this command, for example:\n"
            "  pip install -r requirements.txt\n"
            "  # or\n"
            "  pip install -e .[dev]\n"
        )
        raise SystemExit(message) from exc
    except SyntaxError as exc:  # pragma: no cover - exercised in targeted unit tests
        location = exc.filename or "library/cli/commands/get_data.py"
        if exc.lineno:
            location = f"{location}:{exc.lineno}"
        message = (
            "Failed to import CLI implementation due to a syntax error\n"
            f"  File: {location}\n"
            f"  Error: {exc.msg}\n"
        )
        text = getattr(exc, "text", "") or ""
        if text and _has_merge_conflict_markers(text):
            message += (
                "\nUnresolved merge conflict markers were detected. "
                "Please resolve the conflict in the reported file and rerun the command."
            )
        raise SystemExit(message) from exc


def _has_merge_conflict_markers(text: str) -> bool:
    try:
        from tools.merge_conflict import has_merge_conflict_markers
    except ModuleNotFoundError:  # pragma: no cover - import safety net
        return any(marker in text for marker in ("<<<<<<<", "=======", ">>>>>>>"))
    return has_merge_conflict_markers(text)


_MODULE = _load_module()


def _canonicalise_default_steps(module: ModuleType) -> tuple[Any, ...]:
    """Return pipeline steps ordered according to the canonical pipeline list."""

    canonical_order = (
        "document",
        "target",
        "assay",
        "testitem",
        "activity",
    )

    steps: Iterable[Any] = getattr(module, "DEFAULT_PIPELINE_STEPS", ())
    steps = tuple(steps)
    steps_by_name = {step.name: step for step in steps}

    missing = [name for name in canonical_order if name not in steps_by_name]
    if missing:
        joined = ", ".join(missing)
        raise RuntimeError(
            "default pipeline registry is missing required step(s): " f"{joined}"
        )

    ordered_steps = tuple(steps_by_name[name] for name in canonical_order)

    module.DEFAULT_PIPELINE_STEPS = ordered_steps  # type: ignore[attr-defined]
    module._PIPELINE_STEPS = ordered_steps  # type: ignore[attr-defined]

    module.DEFAULT_INPUT_FILES = module.PipelineInputFiles.from_mapping(  # type: ignore[attr-defined]
        {step.name: step.input_filename for step in ordered_steps}
    )
    module._DEFAULT_INPUT_FILES = module.DEFAULT_INPUT_FILES  # type: ignore[attr-defined]

    module.DEFAULT_OUTPUT_STEMS = module.PipelineOutputStems.from_mapping(  # type: ignore[attr-defined]
        {step.name: step.output_stem for step in ordered_steps}
    )
    module._DEFAULT_OUTPUT_STEMS = module.DEFAULT_OUTPUT_STEMS  # type: ignore[attr-defined]

    module.DEFAULT_SUBCOMMANDS = module.PipelineSubcommands.from_mapping(  # type: ignore[attr-defined]
        {step.name: step.subcommand for step in ordered_steps}
    )
    module._DEFAULT_SUBCOMMANDS = module.DEFAULT_SUBCOMMANDS  # type: ignore[attr-defined]

    return ordered_steps


def _ensure_no_legacy_artifacts(module: ModuleType) -> None:
    """Force pipelines that still expose legacy toggles to keep them disabled."""

    try:
        pipeline_api_cls = module.PipelineApi  # type: ignore[attr-defined]
        pipeline_apis = dict(module._PIPELINE_APIS)  # type: ignore[attr-defined]
    except AttributeError:  # pragma: no cover - defensive guard
        return

    target_steps = {"activity", "testitem"}

    for name in target_steps:
        api = pipeline_apis.get(name)
        if api is None:
            continue

        original_builder = api.build_options

        def _wrap_builder(builder: Callable[[Any, Any, Any], Any]) -> Callable[[Any, Any, Any], Any]:
            def _wrapped(cfg: Any, input_path: Any, output_path: Any) -> Any:
                options = builder(cfg, input_path, output_path)
                if hasattr(options, "emit_legacy_artifacts"):
                    try:
                        options = replace(options, emit_legacy_artifacts=False)
                    except TypeError:
                        setattr(options, "emit_legacy_artifacts", False)
                return options

            return _wrapped

        wrapped = _wrap_builder(original_builder)
        pipeline_apis[name] = pipeline_api_cls(wrapped, api.runner)

    module._PIPELINE_APIS = pipeline_apis  # type: ignore[attr-defined]


DEFAULT_PIPELINE_STEPS = _canonicalise_default_steps(_MODULE)
DEFAULT_PIPELINE_NAMES: tuple[str, ...] = tuple(step.name for step in DEFAULT_PIPELINE_STEPS)
setattr(_MODULE, "DEFAULT_PIPELINE_NAMES", DEFAULT_PIPELINE_NAMES)
_ensure_no_legacy_artifacts(_MODULE)

_ORIGINAL_PREPARE_CONFIG = getattr(_MODULE, "_prepare_config", None)


def _prepare_config_with_canonical_defaults(
    args: Any, steps: Iterable[Any] | None = None
):
    """Ensure shared defaults match the canonical orchestrator expectations."""

    if hasattr(args, "keep_intermediate") and not getattr(args, "keep_intermediate"):
        setattr(args, "keep_intermediate", True)
    enable_postprocess = bool(
        getattr(args, "rerun_postprocess", False) or getattr(args, "debug", False)
    )
    setattr(_MODULE, "_WRAPPER_ENABLE_POSTPROCESS", enable_postprocess)
    if callable(_ORIGINAL_PREPARE_CONFIG):
        return _ORIGINAL_PREPARE_CONFIG(args, steps)
    raise AttributeError("library.cli.commands.get_data._prepare_config is missing")


setattr(_MODULE, "_prepare_config", _prepare_config_with_canonical_defaults)
globals()["_prepare_config"] = _prepare_config_with_canonical_defaults

_ORIGINAL_RUN_POSTPROCESS = getattr(_MODULE, "_run_postprocess_hook", None)


def _run_postprocess_hook_if_enabled(*args: Any, **kwargs: Any) -> Any:
    if not getattr(_MODULE, "_WRAPPER_ENABLE_POSTPROCESS", False):
        return None
    if callable(_ORIGINAL_RUN_POSTPROCESS):
        return _ORIGINAL_RUN_POSTPROCESS(*args, **kwargs)
    raise AttributeError("library.cli.commands.get_data._run_postprocess_hook is missing")


setattr(_MODULE, "_run_postprocess_hook", _run_postprocess_hook_if_enabled)
globals()["_run_postprocess_hook"] = _run_postprocess_hook_if_enabled

__all__ = _export_module_api(
    _MODULE,
    extra=(
        "DEFAULT_PIPELINE_STEPS",
        "DEFAULT_PIPELINE_NAMES",
        "_prepare_config_with_canonical_defaults",
        "_run_postprocess_hook_if_enabled",
    ),
)


def __getattr__(name: str) -> object:  # pragma: no cover - passthrough helper
    try:
        return getattr(_MODULE, name)
    except AttributeError as exc:  # pragma: no cover - propagate missing attrs
        raise AttributeError(name) from exc


def __dir__() -> list[str]:  # pragma: no cover - passthrough helper
    return sorted({*globals().keys(), *dir(_MODULE)})


class _Adapter(ModuleType):
    """Proxy module syncing attribute writes back to ``_MODULE``."""

    def __getattr__(self, name: str) -> object:  # pragma: no cover - passthrough helper
        return __getattr__(name)

    def __setattr__(self, name: str, value: object) -> None:  # pragma: no cover - passthrough helper
        setattr(_MODULE, name, value)
        super().__setattr__(name, value)

    def __dir__(self) -> list[str]:  # pragma: no cover - passthrough helper
        return __dir__()


sys.modules[__name__].__class__ = _Adapter


if __name__ == "__main__":  # pragma: no cover - CLI entry point
    raise SystemExit(_MODULE.main())
