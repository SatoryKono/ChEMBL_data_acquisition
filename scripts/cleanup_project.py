"""Project-wide cleanup utility for the ChEMBL data acquisition pipeline.

The script removes transient artefacts while preserving the canonical ETL
outputs and caches required for deterministic runs.  It is safe to execute in
dry-run mode to preview changes prior to deleting any files.
"""

from __future__ import annotations

import argparse
import re
import shutil
import subprocess
import sys
import time
from collections import defaultdict
from collections.abc import Callable, Sequence
from dataclasses import dataclass, field
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
EXCLUDED_DIR_NAMES = {".git", ".venv"}

CACHE_RETENTION_DAYS = 7
LOG_RETENTION_DAYS = 14

OUTPUT_PURGE_PATTERNS = (
    "*.tmp",
    "*.lock",
    "*_intermediate*",
    "*_debug*",
    "output.*_raw.csv",
)
QC_KEYWORDS = ("qc", "quality", "correlation")
PERSISTENT_CACHE_KEYWORDS = ("chembl", "pubchem", "molecule_parent_catalog")
PIPELINE_TMP_DIRS = ("tmp", "raw", "interim")
TEST_GOLDEN_PREFIX = "expected_"
TEST_GOLDEN_SUFFIX = ".csv"


def _is_under_excluded(path: Path) -> bool:
    try:
        parts = path.relative_to(PROJECT_ROOT).parts
    except ValueError:
        # Outside the repository root; treat as excluded to avoid accidental
        # deletions.
        return True
    return any(part in EXCLUDED_DIR_NAMES for part in parts)


def _load_tracked_paths(root: Path) -> set[str]:
    """Return Git-tracked paths relative to *root*."""

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
                path.unlink(missing_ok=True)
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
            self.remove_path(directory)

    def remove_glob(self, pattern: str, *, base: Path | None = None, files_only: bool = False) -> None:
        search_root = base or PROJECT_ROOT
        for candidate in search_root.rglob(pattern):
            if files_only and candidate.is_dir():
                continue
            self.remove_path(candidate)


CleanupTask = Callable[[CleanupContext], None]


def _cleanup_build_artifacts(ctx: CleanupContext) -> None:
    targets = [PROJECT_ROOT / "build", PROJECT_ROOT / "dist", PROJECT_ROOT / "htmlcov"]
    targets.extend(PROJECT_ROOT.glob("*.egg-info"))
    for target in targets:
        ctx.remove_path(target)


def _cleanup_python_artifacts(ctx: CleanupContext) -> None:
    for pattern in ("__pycache__", "*.pyc", "*.pyo", "*.pyd"):
        ctx.remove_glob(pattern)


def _cleanup_tool_caches(ctx: CleanupContext) -> None:
    for name in (".pytest_cache", ".mypy_cache", ".ruff_cache", ".hypothesis", ".tox"):
        ctx.remove_glob(name)
    for path in (PROJECT_ROOT / ".coverage", PROJECT_ROOT / "coverage.xml"):
        ctx.remove_path(path)


def _cleanup_reports(ctx: CleanupContext) -> None:
    reports_dir = PROJECT_ROOT / "reports"
    keep_dirs = {reports_dir / "archive"}
    if not reports_dir.exists():
        return
    for entry in reports_dir.iterdir():
        if entry in keep_dirs:
            continue
        if entry.is_dir():
            ctx.clean_directory(entry, remove_root=True)
        else:
            ctx.remove_path(entry)


def _cleanup_data_outputs(ctx: CleanupContext) -> None:
    output_dir = PROJECT_ROOT / "data" / "output"
    if not output_dir.exists():
        return

    for pattern in OUTPUT_PURGE_PATTERNS:
        ctx.remove_glob(pattern, base=output_dir, files_only=True)

    dataset_outputs: dict[str, list[tuple[Path, str]]] = defaultdict(list)
    output_regex = re.compile(r"^output\.(?P<dataset>.+)_(?P<date>\d{8})\.csv$")
    for file in output_dir.rglob("output.*_*.csv"):
        if not file.is_file():
            continue
        match = output_regex.match(file.name)
        if not match:
            continue
        dataset_outputs[match.group("dataset")].append((file, match.group("date")))

    kept_output_files: set[Path] = set()
    latest_dates: set[str] = set()
    for dataset, entries in dataset_outputs.items():
        latest_date = max(date for _, date in entries)
        latest_dates.add(latest_date)
        for file, date in entries:
            if date == latest_date:
                kept_output_files.add(file)
            else:
                ctx.remove_path(file)

    if dataset_outputs:
        qc_regex = re.compile(r"(\d{8})")
        for file in output_dir.rglob("*.csv"):
            if not file.is_file() or file in kept_output_files:
                continue
            lower_name = file.name.lower()
            if not any(keyword in lower_name for keyword in QC_KEYWORDS):
                continue
            match = qc_regex.search(file.stem)
            if match and match.group(1) not in latest_dates:
                ctx.remove_path(file)

    # Remove stale smoke-test outputs and archived reports.
    for extra in ("output-smoke", "reports", "logs"):
        ctx.clean_directory(PROJECT_ROOT / "data" / extra, remove_root=True)


def _cleanup_pipeline_caches(ctx: CleanupContext) -> None:
    now = time.time()
    retention_seconds = CACHE_RETENTION_DAYS * 86400
    for base in (PROJECT_ROOT / "cache", PROJECT_ROOT / "data" / "cache"):
        if not base.exists():
            continue
        for pattern in ("*.pkl", "*.json", "*.lock", "*.tmp", "*.bak"):
            for file in base.rglob(pattern):
                if not file.is_file():
                    continue
                lower_name = file.name.lower()
                if lower_name.startswith("pubchem_") and file.suffix == ".pkl":
                    ctx.remove_path(file)
                    continue
                if any(keyword in lower_name for keyword in PERSISTENT_CACHE_KEYWORDS):
                    continue
                if now - file.stat().st_mtime > retention_seconds:
                    ctx.remove_path(file)


def _cleanup_logs(ctx: CleanupContext) -> None:
    logs_dir = PROJECT_ROOT / "logs"
    if not logs_dir.exists():
        return

    now = time.time()
    retention_seconds = LOG_RETENTION_DAYS * 86400
    run_latest = logs_dir / "run_latest.log"

    kept_logs: list[Path] = []
    for path in logs_dir.glob("*.log"):
        if path == run_latest:
            continue
        if now - path.stat().st_mtime > retention_seconds:
            ctx.remove_path(path)
        else:
            kept_logs.append(path)

    for pattern in ("*.old", "*.log.*"):
        for path in logs_dir.glob(pattern):
            if path.is_file():
                ctx.remove_path(path)

    if kept_logs:
        kept_logs.sort(key=lambda item: item.stat().st_mtime)
        rel_dest = None
        if ctx.dry_run:
            try:
                rel_dest = run_latest.relative_to(PROJECT_ROOT)
            except ValueError:
                rel_dest = run_latest
            print(f"[dry-run] Would consolidate {len(kept_logs)} log file(s) into {rel_dest}")
        else:
            run_latest.parent.mkdir(parents=True, exist_ok=True)
            with run_latest.open("w", encoding="utf-8") as destination:
                for index, path in enumerate(kept_logs):
                    with path.open("r", encoding="utf-8", errors="ignore") as source:
                        shutil.copyfileobj(source, destination)
                    if index != len(kept_logs) - 1:
                        destination.write("\n")
        for path in kept_logs:
            if path != run_latest:
                ctx.remove_path(path)


def _cleanup_test_outputs(ctx: CleanupContext) -> None:
    test_data_dir = PROJECT_ROOT / "tests" / "data"
    if not test_data_dir.exists():
        return
    for entry in test_data_dir.iterdir():
        if entry.is_file():
            name = entry.name
            if not (name.startswith(TEST_GOLDEN_PREFIX) and name.endswith(TEST_GOLDEN_SUFFIX)):
                ctx.remove_path(entry)
        else:
            ctx.remove_path(entry)


def _cleanup_tmp_files(ctx: CleanupContext) -> None:
    for pattern in ("*.tmp", "*.lock", "*.bak", "*.swp", "*.swo"):
        ctx.remove_glob(pattern, files_only=True)


def _cleanup_tmp_directories(ctx: CleanupContext) -> None:
    for name in PIPELINE_TMP_DIRS:
        ctx.remove_path(PROJECT_ROOT / name)
        ctx.remove_path(PROJECT_ROOT / "data" / name)


def _cleanup_activity_logs_tmp(ctx: CleanupContext) -> None:
    ctx.clean_directory(PROJECT_ROOT / "scripts" / "activity_logs_tmp" / "output")


CLEANUP_TASKS: dict[str, CleanupTask] = {
    "build": _cleanup_build_artifacts,
    "python": _cleanup_python_artifacts,
    "tooling": _cleanup_tool_caches,
    "reports": _cleanup_reports,
    "data": _cleanup_data_outputs,
    "pipeline-cache": _cleanup_pipeline_caches,
    "logs": _cleanup_logs,
    "tests": _cleanup_test_outputs,
    "tmp": _cleanup_tmp_files,
    "tmp-dirs": _cleanup_tmp_directories,
    "activity-logs": _cleanup_activity_logs_tmp,
}


def _parse_args(argv: Sequence[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Remove transient build artefacts and stale pipeline caches."
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Preview the cleanup actions without deleting anything.",
    )
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
