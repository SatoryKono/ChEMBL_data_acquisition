"""Project-wide cleanup utility.

This helper consolidates the maintenance sweep that removes build artefacts,
pytest leftovers and transient CSV exports created by the local pipelines.
It replaces the ad-hoc shell commands previously embedded in the Makefile and
exposes a `--dry-run` flag so developers can preview the actions before they are
executed.
"""
from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from collections.abc import Callable, Sequence
from dataclasses import dataclass, field
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent

# Directories that must never be traversed or removed by pattern-based sweeps.
EXCLUDED_DIR_NAMES = {".git", ".venv"}


def _is_under_excluded(path: Path) -> bool:
    try:
        parts = path.relative_to(PROJECT_ROOT).parts
    except ValueError:
        # Outside the repository root; treat as excluded to avoid accidental
        # deletions.
        return True
    return any(part in EXCLUDED_DIR_NAMES for part in parts)


def _load_tracked_paths(root: Path) -> set[str]:
    """Return Git-tracked paths relative to *root*.

    The Git CLI is optional – if it is missing or the command fails the helper
    simply treats every file as untracked and proceeds with best-effort cleanup.
    """

    try:
        result = subprocess.run(
            ["git", "ls-files"],
            cwd=root,
            check=True,
            capture_output=True,
            text=True,
        )
    except (FileNotFoundError, subprocess.CalledProcessError):
        return set()
    return {line.strip() for line in result.stdout.splitlines() if line.strip()}


@dataclass
class CleanupContext:
    dry_run: bool
    tracked_paths: set[str]
    removed: list[Path] = field(default_factory=list)
    skipped: list[tuple[Path, str]] = field(default_factory=list)

    def _relative(self, path: Path) -> Path:
        return path.relative_to(PROJECT_ROOT)

    def _is_tracked(self, path: Path) -> bool:
        rel = self._relative(path).as_posix()
        if rel in self.tracked_paths:
            return True
        prefix = rel + "/"
        return any(candidate.startswith(prefix) for candidate in self.tracked_paths)

    def remove_path(self, path: Path) -> None:
        if _is_under_excluded(path):
            return
        if not path.exists() and not path.is_symlink():
            return
        try:
            rel = self._relative(path)
        except ValueError:
            self.skipped.append((path, "outside project root"))
            return
        if self._is_tracked(path):
            self.skipped.append((path, "tracked by git"))
            return
        description = f"{rel}{'/' if path.is_dir() and not path.is_symlink() else ''}"
        if self.dry_run:
            print(f"[dry-run] Would remove {description}")
            self.removed.append(path)
            return
        try:
            if path.is_dir() and not path.is_symlink():
                shutil.rmtree(path)
            else:
                path.unlink(missing_ok=True)  # Python 3.8+: guard symlink races.
        except OSError as exc:
            self.skipped.append((path, f"failed to remove ({exc})"))
            return
        print(f"Removed {description}")
        self.removed.append(path)

    def clean_directory(self, directory: Path, *, remove_root: bool = False) -> None:
        if _is_under_excluded(directory) or not directory.exists():
            return
        try:
            self._relative(directory)
        except ValueError:
            self.skipped.append((directory, "outside project root"))
            return
        if not directory.is_dir():
            self.remove_path(directory)
            return
        for child in sorted(directory.iterdir(), key=lambda item: item.name):
            self.remove_path(child)
        if remove_root and directory.exists():
            # Attempt to remove the now-empty directory; skip if tracked.
            self.remove_path(directory)

    def remove_glob(self, pattern: str, *, files_only: bool = False) -> None:
        for candidate in PROJECT_ROOT.rglob(pattern):
            if files_only and candidate.is_dir():
                continue
            self.remove_path(candidate)


CleanupTask = Callable[[CleanupContext], None]


def _cleanup_build(ctx: CleanupContext) -> None:
    targets = [PROJECT_ROOT / "build", PROJECT_ROOT / "dist", PROJECT_ROOT / "htmlcov"]
    targets.extend(PROJECT_ROOT.glob("*.egg-info"))
    for target in targets:
        ctx.remove_path(target)


def _cleanup_python_artifacts(ctx: CleanupContext) -> None:
    for pattern in ("__pycache__", "*.pyc", "*.pyo", "*.pyd"):
        ctx.remove_glob(pattern)


def _cleanup_tool_caches(ctx: CleanupContext) -> None:
    for name in (".pytest_cache", ".mypy_cache", ".ruff_cache", ".hypothesis", ".coverage", "coverage.xml"):
        for path in PROJECT_ROOT.glob(name):
            ctx.remove_path(path)
    for cache_dir in (".pytest_cache", ".mypy_cache", ".ruff_cache", ".hypothesis", ".tox"):
        for path in PROJECT_ROOT.rglob(cache_dir):
            ctx.remove_path(path)


def _cleanup_reports(ctx: CleanupContext) -> None:
    report_files = [
        PROJECT_ROOT / "reports" / "test_report.json",
        PROJECT_ROOT / "reports" / "test_summary.md",
        PROJECT_ROOT / "reports" / "pytest_raw_report.json",
    ]
    for path in report_files:
        ctx.remove_path(path)
    ctx.clean_directory(PROJECT_ROOT / "reports" / "coverage", remove_root=True)
    ctx.clean_directory(PROJECT_ROOT / "reports" / "tests" / "artifacts", remove_root=True)


def _cleanup_logs(ctx: CleanupContext) -> None:
    logs_dir = PROJECT_ROOT / "logs"
    ctx.clean_directory(logs_dir)
    for pattern in ("*.log", "*.log.json"):
        ctx.remove_glob(pattern)


def _cleanup_data_outputs(ctx: CleanupContext) -> None:
    ctx.clean_directory(PROJECT_ROOT / "data" / "output", remove_root=True)
    ctx.clean_directory(PROJECT_ROOT / "data" / "output-smoke", remove_root=True)
    ctx.clean_directory(PROJECT_ROOT / "data" / "reports", remove_root=True)
    ctx.clean_directory(PROJECT_ROOT / "data" / "logs", remove_root=True)
    ctx.remove_glob("data/**/*.tmp")


def _cleanup_test_outputs(ctx: CleanupContext) -> None:
    ctx.clean_directory(PROJECT_ROOT / "tests" / "data", remove_root=True)
    ctx.clean_directory(PROJECT_ROOT / "tests" / "output", remove_root=True)
    ctx.clean_directory(PROJECT_ROOT / "tests" / "reports", remove_root=True)


def _cleanup_tmp_files(ctx: CleanupContext) -> None:
    for pattern in ("*.tmp", "*.lock", "*.bak", "*.swp", "*.swo"):
        ctx.remove_glob(pattern, files_only=True)


def _cleanup_activity_logs_tmp(ctx: CleanupContext) -> None:
    ctx.clean_directory(PROJECT_ROOT / "scripts" / "activity_logs_tmp" / "output")


CLEANUP_TASKS: dict[str, CleanupTask] = {
    "build": _cleanup_build,
    "python": _cleanup_python_artifacts,
    "caches": _cleanup_tool_caches,
    "reports": _cleanup_reports,
    "logs": _cleanup_logs,
    "data": _cleanup_data_outputs,
    "tests": _cleanup_test_outputs,
    "tmp": _cleanup_tmp_files,
    "activity-logs": _cleanup_activity_logs_tmp,
}


def _parse_args(argv: Sequence[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Remove transient build artefacts and cache directories.")
    parser.add_argument("--dry-run", action="store_true", help="Preview the cleanup actions without deleting anything.")
    parser.add_argument(
        "--all",
        action="store_true",
        help="Execute every cleanup category (default when no specific categories are listed).",
    )
    parser.add_argument(
        "--categories",
        nargs="+",
        choices=sorted(CLEANUP_TASKS.keys()),
        help="Explicit list of cleanup categories to run (overrides --all).",
    )
    parser.add_argument(
        "--list",
        action="store_true",
        help="Show the available cleanup categories and exit.",
    )
    return parser.parse_args(argv)


def _resolve_categories(args: argparse.Namespace) -> list[str]:
    if args.list:
        for name in sorted(CLEANUP_TASKS.keys()):
            print(name)
        return []
    if args.categories:
        return list(dict.fromkeys(args.categories))
    if args.all or not args.categories:
        return list(CLEANUP_TASKS.keys())
    return []


def _print_summary(ctx: CleanupContext) -> None:
    print()
    unique_removed = {path.resolve() for path in ctx.removed}
    print(f"Cleanup complete: {len(unique_removed)} paths targeted.")
    if ctx.skipped:
        print("Skipped:")
        for path, reason in ctx.skipped:
            try:
                rel = path.relative_to(PROJECT_ROOT)
            except ValueError:
                rel = path
            print(f"  - {rel}: {reason}")


def main(argv: Sequence[str] | None = None) -> int:
    args = _parse_args(argv or sys.argv[1:])
    categories = _resolve_categories(args)
    if args.list:
        return 0
    if not categories:
        print("No cleanup categories selected.")
        return 0
    tracked = _load_tracked_paths(PROJECT_ROOT)
    context = CleanupContext(dry_run=args.dry_run, tracked_paths=tracked)
    for name in categories:
        task = CLEANUP_TASKS[name]
        print(f"Running cleanup category: {name}")
        task(context)
    _print_summary(context)
    failed = any("failed" in reason for _, reason in context.skipped)
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
