"""Internal helper to run profiling in a subprocess.

This module is executed via ``python -m tools.profile_worker`` and is not
intended to be invoked directly by end users.  It performs the heavy lifting
of running :mod:`cProfile` (and optionally :mod:`line_profiler`) against one of
our data acquisition entry points and serialises the metrics to JSON so that
:mod:`tools.performance_profiler` can present them in human readable form.
"""

from __future__ import annotations

import argparse
import contextlib
import importlib
import importlib.util
import io
import json
import os
import resource
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Sequence

import cProfile
import pstats


@dataclass
class ProfilingResult:
    """Container for the profiling metrics returned by this worker."""

    module: str
    argv: list[str]
    target_root: str
    return_code: int
    wall_time: float
    cpu_user_time: float
    cpu_system_time: float
    max_rss_kb: int
    io_in_blocks: int
    io_out_blocks: int
    profile_output: str
    line_profile_output: str | None
    stats: list[dict[str, Any]]


def _normalise_module(module: str) -> str:
    module = module.removesuffix(".py").replace("/", ".")
    if module.startswith("scripts."):
        return module
    if module == "scripts":
        raise ValueError("Module name must refer to a concrete script")
    if module.startswith("scripts_"):
        # Support callers passing e.g. ``scripts_get_activity_data``
        module = module.replace("scripts_", "scripts.", 1)
    if module.startswith("get_"):
        return f"scripts.{module}"
    return module


def _load_main(module_name: str) -> Any:
    module = importlib.import_module(module_name)
    main = getattr(module, "main", None)
    if main is None:
        raise RuntimeError(f"Module {module_name!r} does not expose a main() function")
    return main


def _prepare_environment(target_root: Path) -> None:
    sys.path.insert(0, str(target_root))
    os.chdir(target_root)
    importlib.invalidate_caches()


def _run_with_line_profiler(main: Any, argv: Sequence[str]) -> tuple[int, str | None]:
    spec = importlib.util.find_spec("line_profiler")
    if spec is None:  # pragma: no cover - optional dependency
        return main(argv), None

    line_profiler = importlib.import_module("line_profiler")
    profiler = line_profiler.LineProfiler()
    profiler.add_function(main)
    return_code = profiler.runcall(main, argv)
    output_buffer = io.StringIO()
    profiler.print_stats(stream=output_buffer)
    return return_code, output_buffer.getvalue()


def _build_stats(profile_path: Path) -> list[dict[str, Any]]:
    stats = pstats.Stats(str(profile_path))
    stats.sort_stats(pstats.SortKey.CUMULATIVE)
    results: list[dict[str, Any]] = []
    for func, (cc, nc, tt, ct, callers) in stats.stats.items():
        filename, line, name = func
        results.append(
            {
                "function": name,
                "filename": filename,
                "lineno": line,
                "ccalls": cc,
                "ncalls": nc,
                "tottime": tt,
                "cumtime": ct,
            }
        )
    results.sort(key=lambda item: item["cumtime"], reverse=True)
    return results


def _dump_json(result: ProfilingResult, destination: Path) -> None:
    destination.write_text(
        json.dumps(
            {
                "module": result.module,
                "argv": result.argv,
                "target_root": result.target_root,
                "return_code": result.return_code,
                "wall_time": result.wall_time,
                "cpu_user_time": result.cpu_user_time,
                "cpu_system_time": result.cpu_system_time,
                "max_rss_kb": result.max_rss_kb,
                "io_in_blocks": result.io_in_blocks,
                "io_out_blocks": result.io_out_blocks,
                "profile_output": result.profile_output,
                "line_profile_output": result.line_profile_output,
                "stats": result.stats,
            },
            indent=2,
        )
    )


def _parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Execute profiling for a data script")
    parser.add_argument("--module", required=True, help="Target module or script path")
    parser.add_argument(
        "--target-root",
        type=Path,
        default=Path.cwd(),
        help="Root directory containing the target package (defaults to CWD)",
    )
    parser.add_argument(
        "--profile-output",
        type=Path,
        required=True,
        help="File where the binary cProfile dump will be written",
    )
    parser.add_argument(
        "--line-profile-output",
        type=Path,
        help="Optional file capturing line_profiler statistics",
    )
    parser.add_argument(
        "--result-json",
        type=Path,
        required=True,
        help="JSON file where structured metrics will be saved",
    )
    parser.add_argument(
        "--with-line-profiler",
        action="store_true",
        help="Attempt to run line_profiler if available",
    )
    parser.add_argument(
        "script_args",
        nargs=argparse.REMAINDER,
        help="Arguments passed through to the target script",
    )
    return parser.parse_args(argv)


def _normalise_argv(raw: Iterable[str]) -> list[str]:
    argv = list(raw)
    if argv and argv[0] == "--":
        argv = argv[1:]
    return argv


def main(argv: Sequence[str] | None = None) -> int:
    args = _parse_args(argv)
    module_name = _normalise_module(args.module)
    script_argv = _normalise_argv(args.script_args)
    target_root = args.target_root.resolve()

    profile_output = args.profile_output.resolve()
    profile_output.parent.mkdir(parents=True, exist_ok=True)

    line_profile_text: str | None = None
    if args.line_profile_output is not None:
        line_profile_output = args.line_profile_output.resolve()
        line_profile_output.parent.mkdir(parents=True, exist_ok=True)
    else:
        line_profile_output = None

    original_cwd = Path.cwd()
    original_sys_path = list(sys.path)

    with contextlib.ExitStack() as stack:
        stack.callback(lambda: os.chdir(original_cwd))
        stack.callback(lambda: sys.path.clear() or sys.path.extend(original_sys_path))

        _prepare_environment(target_root)

        main_callable = _load_main(module_name)

        profiler = cProfile.Profile()
        start_usage = resource.getrusage(resource.RUSAGE_SELF)
        start_time = time.perf_counter()

        profiler.enable()
        try:
            if args.with_line_profiler and line_profile_output is not None:
                return_code, line_profile_text = _run_with_line_profiler(main_callable, script_argv)
            else:
                return_code = main_callable(script_argv)
        except SystemExit as exc:
            return_code = int(exc.code or 0)
        finally:
            profiler.disable()

        wall_time = time.perf_counter() - start_time
        end_usage = resource.getrusage(resource.RUSAGE_SELF)

    profiler.dump_stats(str(profile_output))

    if line_profile_output is not None and line_profile_text is not None:
        line_profile_output.write_text(line_profile_text)

    stats = _build_stats(profile_output)

    result = ProfilingResult(
        module=module_name,
        argv=script_argv,
        target_root=str(target_root),
        return_code=int(return_code or 0),
        wall_time=wall_time,
        cpu_user_time=end_usage.ru_utime - start_usage.ru_utime,
        cpu_system_time=end_usage.ru_stime - start_usage.ru_stime,
        max_rss_kb=end_usage.ru_maxrss,
        io_in_blocks=end_usage.ru_inblock - start_usage.ru_inblock,
        io_out_blocks=end_usage.ru_oublock - start_usage.ru_oublock,
        profile_output=str(profile_output),
        line_profile_output=str(line_profile_output) if line_profile_output else None,
        stats=stats,
    )

    _dump_json(result, args.result_json.resolve())
    return 0


if __name__ == "__main__":  # pragma: no cover - script entry point
    raise SystemExit(main())
