"""Profiling and benchmarking utilities for ChEMBL data scripts."""

from __future__ import annotations

import argparse
import contextlib
import datetime as dt
import json
import os
import shutil
import subprocess
import sys
import sysconfig
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Sequence

REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_REPORTS_DIR = REPO_ROOT / "reports"
DEFAULT_PROFILES_DIR = DEFAULT_REPORTS_DIR / "profiles"
DEFAULT_OPT_REPORT_TEMPLATE = "performance_optimization_{date}.md"


class ProfilingError(RuntimeError):
    """Raised when the profiling worker fails."""


@dataclass
class ProfilingArtefacts:
    result: dict[str, Any]
    profile_path: Path
    line_profile_path: Path | None
    report_path: Path
    json_path: Path


def _normalise_module(value: str) -> str:
    value = value.removesuffix(".py").replace("/", ".")
    if value.startswith("scripts."):
        return value
    if value.startswith("scripts_"):
        return value.replace("scripts_", "scripts.", 1)
    if value.startswith("get_"):
        return f"scripts.{value}"
    return value


def _ensure_reports_dir() -> None:
    DEFAULT_REPORTS_DIR.mkdir(parents=True, exist_ok=True)
    DEFAULT_PROFILES_DIR.mkdir(parents=True, exist_ok=True)


def _run_worker(
    module: str,
    script_args: Sequence[str],
    *,
    target_root: Path,
    profile_destination: Path,
    enable_line_profiler: bool,
    line_profile_destination: Path | None,
    result_destination: Path,
) -> dict[str, Any]:
    cmd = [
        sys.executable,
        "-m",
        "tools.profile_worker",
        "--module",
        module,
        "--target-root",
        str(target_root),
        "--profile-output",
        str(profile_destination),
        "--result-json",
        str(result_destination),
    ]
    if enable_line_profiler:
        cmd.append("--with-line-profiler")
    if line_profile_destination is not None:
        cmd.extend(["--line-profile-output", str(line_profile_destination)])
    if script_args:
        cmd.append("--")
        cmd.extend(script_args)

    process = subprocess.run(cmd, capture_output=True, text=True, check=False)
    if process.returncode != 0:
        raise ProfilingError(
            "Profiling worker failed with exit code "
            f"{process.returncode}: {process.stderr.strip()}"
        )

    result_payload = json.loads(Path(result_destination).read_text())
    result_payload["stdout"] = process.stdout
    result_payload["stderr"] = process.stderr
    return result_payload


def _format_path(path: str | os.PathLike[str], base: Path) -> str:
    candidate = Path(path)
    try:
        rel = candidate.resolve().relative_to(base.resolve())
        return str(rel)
    except Exception:  # pragma: no cover - defensive
        return str(candidate)


def _operation_type(filename: str, project_root: Path) -> str:
    if filename.startswith("<"):
        return "built-in"

    candidate = Path(filename)
    stdlib = Path(sysconfig.get_paths()["stdlib"]).resolve()
    try:
        resolved = candidate.resolve()
    except Exception:  # pragma: no cover - fallback for transient files
        return "external"

    try:
        if resolved.is_relative_to(project_root.resolve()):
            if "scripts" in resolved.parts:
                return "cli"
            if "library" in resolved.parts:
                return "library"
            return "project"
    except Exception:  # pragma: no cover - safety for non-existent paths
        pass

    try:
        if resolved.is_relative_to(stdlib):
            return "stdlib"
    except Exception:  # pragma: no cover - optional on exotic platforms
        pass

    if "site-packages" in resolved.parts:
        return "third-party"

    return "external"


def _build_summary_table(
    stats: list[dict[str, Any]],
    wall_time: float,
    project_root: Path,
    *,
    top_n: int = 20,
) -> str:
    header = "| Function | Time (s) | % of wall | Operation type |\n|---|---:|---:|---|"
    rows = []
    for entry in stats[:top_n]:
        function = entry["function"]
        filename = entry["filename"]
        lineno = entry["lineno"]
        cumtime = float(entry["cumtime"])
        percent = 0.0 if wall_time == 0 else (cumtime / wall_time) * 100.0
        location = _format_path(filename, project_root)
        op_type = _operation_type(filename, project_root)
        rows.append(
            f"| `{function}` ({location}:{lineno}) | {cumtime:.4f} | {percent:6.2f} | {op_type} |"
        )

    if not rows:
        rows.append("| _no samples recorded_ | 0.0000 | 0.00 | n/a |")

    return "\n".join([header, *rows])


def _render_profile_report(
    *,
    script: str,
    args: Sequence[str],
    result: dict[str, Any],
    report_path: Path,
    top_n: int,
) -> None:
    project_root = Path(result["target_root"]).resolve()
    summary_table = _build_summary_table(result["stats"], float(result["wall_time"]), project_root, top_n=top_n)

    lines = [
        f"# Profiling report for `{script}`",
        "",
        f"- Arguments: `{' '.join(args) if args else '<none>'}`",
        f"- Target root: `{project_root}`",
        f"- Wall time: {result['wall_time']:.4f} s",
        f"- CPU (user): {result['cpu_user_time']:.4f} s",
        f"- CPU (system): {result['cpu_system_time']:.4f} s",
        f"- Max RSS: {result['max_rss_kb'] / 1024:.2f} MiB",
        f"- Block I/O (in/out): {result['io_in_blocks']} / {result['io_out_blocks']}",
        f"- Profile artefact: `{result['profile_output']}`",
    ]
    if result.get("line_profile_output"):
        lines.append(f"- Line-profiler artefact: `{result['line_profile_output']}`")
    lines.extend(
        [
            "",
            "## Top functions by cumulative time",
            "",
            summary_table,
        ]
    )
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text("\n".join(lines))


def _profile_once(
    *,
    script: str,
    args: Sequence[str],
    target_root: Path,
    enable_line_profiler: bool,
    top_n: int,
    artefact_prefix: str,
) -> ProfilingArtefacts:
    timestamp = dt.datetime.utcnow().strftime("%Y%m%dT%H%M%SZ")
    profile_path = DEFAULT_PROFILES_DIR / f"{artefact_prefix}_{timestamp}.prof"
    json_path = DEFAULT_PROFILES_DIR / f"{artefact_prefix}_{timestamp}.json"
    line_profile_path = (
        DEFAULT_PROFILES_DIR / f"{artefact_prefix}_{timestamp}.lprof"
        if enable_line_profiler
        else None
    )

    result = _run_worker(
        _normalise_module(script),
        list(args),
        target_root=target_root,
        profile_destination=profile_path,
        enable_line_profiler=enable_line_profiler,
        line_profile_destination=line_profile_path,
        result_destination=json_path,
    )

    report_path = DEFAULT_PROFILES_DIR / f"{artefact_prefix}_{timestamp}.md"
    _render_profile_report(
        script=script,
        args=args,
        result=result,
        report_path=report_path,
        top_n=top_n,
    )
    return ProfilingArtefacts(result, profile_path, line_profile_path, report_path, json_path)


def _profile_command(args: argparse.Namespace) -> int:
    _ensure_reports_dir()
    artefact_prefix = args.script.replace("/", "_").replace(".", "_")
    artefacts = _profile_once(
        script=args.script,
        args=args.script_args,
        target_root=Path(args.target_root or REPO_ROOT),
        enable_line_profiler=args.line_profiler,
        top_n=args.top,
        artefact_prefix=artefact_prefix,
    )
    if args.output:
        destination = Path(args.output)
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(artefacts.report_path, destination)
        print(f"Profiling report written to {destination}")
    else:
        print(f"Profiling report written to {artefacts.report_path}")
    return int(artefacts.result.get("return_code", 0))


def _format_metrics_row(label: str, result: dict[str, Any]) -> str:
    return (
        f"| {label} | {result['wall_time']:.4f} | {result['cpu_user_time']:.4f} | "
        f"{result['cpu_system_time']:.4f} | {result['max_rss_kb'] / 1024:.2f} | "
        f"{result['io_in_blocks']} | {result['io_out_blocks']} |"
    )


def _resolve_inputs(inputs: Iterable[str], destination: Path) -> list[str]:
    staged: list[str] = []
    destination.mkdir(parents=True, exist_ok=True)
    for item in inputs:
        source = (REPO_ROOT / item).resolve() if not os.path.isabs(item) else Path(item)
        if not source.exists():
            raise FileNotFoundError(f"Input path {source} does not exist")
        target = destination / source.name
        if source.is_dir():
            if target.exists():
                shutil.rmtree(target)
            shutil.copytree(source, target)
        else:
            shutil.copy2(source, target)
        staged.append(str(target))
    return staged


def _render_benchmark_report(
    *,
    script: str,
    args: Sequence[str],
    before: ProfilingArtefacts,
    after: ProfilingArtefacts,
    before_label: str,
    after_label: str,
    optimisations: Sequence[str],
    report_path: Path,
    top_n: int,
    inputs_before: Sequence[str],
    inputs_after: Sequence[str],
) -> None:
    speedup = (
        float(before.result["wall_time"]) / float(after.result["wall_time"])
        if float(after.result["wall_time"]) > 0
        else float("inf")
    )
    metrics_table = "\n".join(
        [
            "| Run | Wall time (s) | CPU user (s) | CPU sys (s) | Max RSS (MiB) | In blocks | Out blocks |",
            "|---|---:|---:|---:|---:|---:|---:|",
            _format_metrics_row(before_label, before.result),
            _format_metrics_row(after_label, after.result),
        ]
    )

    summary_before = _build_summary_table(
        before.result["stats"],
        float(before.result["wall_time"]),
        Path(before.result["target_root"]).resolve(),
        top_n=top_n,
    )
    summary_after = _build_summary_table(
        after.result["stats"],
        float(after.result["wall_time"]),
        Path(after.result["target_root"]).resolve(),
        top_n=top_n,
    )

    lines = [
        f"# Performance optimisation report for `{script}`",
        "",
        f"- Arguments: `{' '.join(args) if args else '<none>'}`",
        f"- Baseline: {before_label}",
        f"- Optimised: {after_label}",
        f"- Speed-up (wall): {speedup:.2f}×",
        "",
        "## Optimisations applied",
        "",
    ]
    if optimisations:
        lines.extend(f"- {item}" for item in optimisations)
    else:
        lines.append("- _Document optimisation steps here_")

    lines.extend(
        [
            "",
            "## Runtime metrics",
            "",
            metrics_table,
            "",
            "## Baseline hot paths",
            "",
            summary_before,
            "",
            "## Optimised hot paths",
            "",
            summary_after,
            "",
            "## Artefacts",
            "",
            f"- Baseline profile: `{before.profile_path}`",
            f"- Baseline JSON: `{before.json_path}`",
            f"- Optimised profile: `{after.profile_path}`",
            f"- Optimised JSON: `{after.json_path}`",
        ]
    )
    if before.line_profile_path:
        lines.append(f"- Baseline line profile: `{before.line_profile_path}`")
    if after.line_profile_path:
        lines.append(f"- Optimised line profile: `{after.line_profile_path}`")

    if inputs_before or inputs_after:
        lines.extend(
            [
                "",
                "## Input staging",
                "",
            ]
        )
        if inputs_before:
            lines.append("**Baseline inputs**")
            lines.extend(f"- `{path}`" for path in inputs_before)
        if inputs_after:
            if inputs_before:
                lines.append("")
            lines.append("**Optimised inputs**")
            lines.extend(f"- `{path}`" for path in inputs_after)

    report_path.parent.mkdir(parents=True, exist_ok=True)
    if report_path.exists():
        report_path.write_text(report_path.read_text() + "\n\n---\n\n" + "\n".join(lines))
    else:
        report_path.write_text("\n".join(lines))


def _benchmark_command(args: argparse.Namespace) -> int:
    _ensure_reports_dir()
    artefact_prefix_base = args.script.replace("/", "_").replace(".", "_")
    top_n = args.top
    before_label: str
    after_label: str

    staged_inputs_before: list[str] = []
    staged_inputs_after: list[str] = []

    with contextlib.ExitStack() as stack:
        if args.before_path:
            before_root = Path(args.before_path).resolve()
            before_label = f"{before_root} ({_git('rev-parse', 'HEAD', cwd=before_root)[:8]})"
        else:
            before_ref = args.before_ref or "HEAD^"
            before_commit = _git("rev-parse", before_ref)
            worktree_cm = stack.enter_context(_temporary_worktree(before_ref))
            before_root = worktree_cm
            before_label = f"{before_ref} ({before_commit[:8]})"

        after_root = Path(args.after_path).resolve() if args.after_path else REPO_ROOT
        after_commit = _git("rev-parse", "HEAD", cwd=after_root)
        after_label = f"HEAD ({after_commit[:8]})"

        if args.input:
            staging_base = Path(tempfile.mkdtemp(prefix="perf-inputs-"))
            stack.callback(lambda: shutil.rmtree(staging_base, ignore_errors=True))
            staged_inputs_before = _resolve_inputs(args.input, staging_base / "before")
            staged_inputs_after = _resolve_inputs(args.input, staging_base / "after")

            def _restore_env(key: str, previous: str | None) -> None:
                if previous is None:
                    os.environ.pop(key, None)
                else:
                    os.environ[key] = previous

            prev_before = os.environ.get("CHEMBL_PERF_INPUT_BEFORE")
            prev_after = os.environ.get("CHEMBL_PERF_INPUT_AFTER")
            os.environ["CHEMBL_PERF_INPUT_BEFORE"] = ":".join(staged_inputs_before)
            os.environ["CHEMBL_PERF_INPUT_AFTER"] = ":".join(staged_inputs_after)
            stack.callback(lambda: _restore_env("CHEMBL_PERF_INPUT_BEFORE", prev_before))
            stack.callback(lambda: _restore_env("CHEMBL_PERF_INPUT_AFTER", prev_after))

        before = _profile_once(
            script=args.script,
            args=args.script_args,
            target_root=before_root,
            enable_line_profiler=args.line_profiler,
            top_n=top_n,
            artefact_prefix=f"baseline_{artefact_prefix_base}",
        )
        after = _profile_once(
            script=args.script,
            args=args.script_args,
            target_root=after_root,
            enable_line_profiler=args.line_profiler,
            top_n=top_n,
            artefact_prefix=f"optimised_{artefact_prefix_base}",
        )

    optimisations = args.optimisation or []
    date_str = (args.date or dt.datetime.utcnow().date().isoformat()).replace("-", "")
    report_name = DEFAULT_OPT_REPORT_TEMPLATE.format(date=date_str)
    report_path = DEFAULT_REPORTS_DIR / report_name
    _render_benchmark_report(
        script=args.script,
        args=args.script_args,
        before=before,
        after=after,
        before_label=before_label,
        after_label=after_label,
        optimisations=optimisations,
        report_path=report_path,
        top_n=top_n,
        inputs_before=staged_inputs_before,
        inputs_after=staged_inputs_after,
    )
    print(f"Benchmark report written to {report_path}")
    return int(after.result.get("return_code", 0))


def _temporary_worktree(ref: str) -> contextlib.AbstractContextManager[Path]:
    @contextlib.contextmanager
    def manager() -> Iterable[Path]:
        workdir = Path(tempfile.mkdtemp(prefix="perf-before-"))
        try:
            subprocess.run(
                ["git", "worktree", "add", "--detach", str(workdir), ref],
                cwd=REPO_ROOT,
                check=True,
                capture_output=True,
                text=True,
            )
            yield workdir
        finally:
            subprocess.run(
                ["git", "worktree", "remove", "--force", str(workdir)],
                cwd=REPO_ROOT,
                check=True,
                capture_output=True,
                text=True,
            )
            shutil.rmtree(workdir, ignore_errors=True)

    return manager()


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Profile and benchmark data acquisition scripts")
    subparsers = parser.add_subparsers(dest="command", required=True)

    profile_parser = subparsers.add_parser("profile", help="Run cProfile against a script")
    profile_parser.add_argument("script", help="Script module or path under scripts/")
    profile_parser.add_argument("script_args", nargs=argparse.REMAINDER, help="Arguments for the script")
    profile_parser.add_argument("--top", type=int, default=20, help="Number of functions to include in the summary table")
    profile_parser.add_argument("--line-profiler", action="store_true", help="Attempt to capture line_profiler statistics")
    profile_parser.add_argument("--output", type=Path, help="Optional destination for the Markdown report")
    profile_parser.add_argument("--target-root", type=Path, help="Alternate repository root (defaults to current tree)")
    profile_parser.set_defaults(func=_profile_command)

    bench_parser = subparsers.add_parser("benchmark", help="Run before/after profiling and emit optimisation report")
    bench_parser.add_argument("script", help="Script module or path under scripts/")
    bench_parser.add_argument("script_args", nargs=argparse.REMAINDER, help="Arguments passed to the script")
    bench_parser.add_argument("--before-ref", help="Git ref used for the baseline profile (defaults to HEAD^)")
    bench_parser.add_argument("--before-path", type=Path, help="Explicit filesystem path used for the baseline run")
    bench_parser.add_argument("--after-path", type=Path, help="Alternate path used for the optimised run (defaults to current repo)")
    bench_parser.add_argument("--top", type=int, default=20, help="Number of functions to include in the summary tables")
    bench_parser.add_argument("--line-profiler", action="store_true", help="Attempt to capture line_profiler statistics for both runs")
    bench_parser.add_argument("--input", action="append", default=[], help="Input files or directories copied for reproducible runs")
    bench_parser.add_argument("--optimisation", action="append", help="Describe an optimisation applied in the optimised build")
    bench_parser.add_argument("--date", help="Override the YYYYMMDD stamp used in the report filename")
    bench_parser.set_defaults(func=_benchmark_command)

    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)
    try:
        return args.func(args)
    except ProfilingError as exc:  # pragma: no cover - CLI level error handler
        parser.error(str(exc))
        return 2


if __name__ == "__main__":  # pragma: no cover - CLI entry point
    raise SystemExit(main())
