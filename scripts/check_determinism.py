#!/usr/bin/env python3
"""Best-effort determinism smoke test for data acquisition pipelines.

The helper runs the activity export twice and compares the artefacts. When
invoked with ``--dry-run`` it hashes the combined stdout/stderr logs instead of
requiring CSV outputs, making it possible to validate deterministic planning in
CI environments where writes are forbidden.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import TYPE_CHECKING

# ruff: noqa: E402  # bootstrap alters import order for script compatibility
if TYPE_CHECKING:
    from . import _bootstrap as _bootstrap_module
elif __package__ in {None, ""}:
    import _bootstrap as _bootstrap_module  # pragma: no cover - CLI fallback
else:  # pragma: no cover - executed when imported as a package module
    from . import _bootstrap as _bootstrap_module

bootstrap_cli = _bootstrap_module.bootstrap_cli
ensure_project_root = _bootstrap_module.ensure_project_root

bootstrap_cli(__package__, __file__)
del bootstrap_cli
del _bootstrap_module


def _hash_file(path: Path) -> str:
    hasher = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            hasher.update(chunk)
    return hasher.hexdigest()


def _hash_process_output(result: subprocess.CompletedProcess[str]) -> str:
    hasher = hashlib.sha256()
    for stream in (result.stdout or "", result.stderr or ""):
        hasher.update(stream.encode("utf-8"))
        hasher.update(b"\0")
    return hasher.hexdigest()


def _metadata_path(csv_path: Path) -> Path:
    return csv_path.with_suffix(csv_path.suffix + ".meta.yaml")


def _run_activity(
    limit: int,
    destination: Path,
    input_csv: Path,
    *,
    dry_run: bool,
    timeout: float | None = None,
    offline: bool = False,
    fixtures_dir: Path | None = None,
) -> subprocess.CompletedProcess[str]:
    env = os.environ.copy()
    env.setdefault("PYTHONHASHSEED", "0")
    repo_root = Path(__file__).resolve().parents[1]
    if offline:
        if fixtures_dir is None:
            raise ValueError("fixtures_dir must be provided when offline is enabled")
        env["CHEMBL_DA_OFFLINE"] = "1"
        cmd = [
            sys.executable,
            "-m",
            "tests.helpers.activity_offline_cli",
            "--fixtures-dir",
            str(fixtures_dir),
            "--destination",
            str(destination),
            "--input",
            str(input_csv),
            "--limit",
            str(limit),
        ]
    else:
        cmd = [
            sys.executable,
            "-m",
            "scripts.get_activity_data",
            "--limit",
            str(limit),
            "--final-out",
            str(destination),
            "--input",
            str(input_csv),
        ]
        if dry_run:
            cmd.append("--dry-run")
    return subprocess.run(
        cmd,
        text=True,
        capture_output=True,
        env=env,
        cwd=str(repo_root),
        timeout=timeout,
    )


def _default_input_csv(tmp_dir: Path) -> Path:
    repo_root = Path(__file__).resolve().parents[1]
    candidate = repo_root / "data" / "input" / "activity.csv"
    if candidate.exists():
        return candidate

    fallback = tmp_dir / "activity.csv"
    with fallback.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["activity_chembl_id"])
        writer.writerow(["ACT1"])
        writer.writerow(["ACT2"])
    return fallback


def _report_process_failure(
    label: str, result: subprocess.CompletedProcess[str]
) -> None:
    sys.stderr.write(f"{label} failed with exit code {result.returncode}\n")
    if result.stdout:
        sys.stderr.write("stdout:\n")
        sys.stderr.write(result.stdout)
        if not result.stdout.endswith("\n"):
            sys.stderr.write("\n")
    if result.stderr:
        sys.stderr.write("stderr:\n")
        sys.stderr.write(result.stderr)
        if not result.stderr.endswith("\n"):
            sys.stderr.write("\n")


def _report_process_timeout(label: str, exc: subprocess.TimeoutExpired) -> None:
    timeout_value = exc.timeout
    sys.stderr.write(f"{label} timed out after {timeout_value} seconds\n")

    output = getattr(exc, "output", None)
    if output:
        sys.stderr.write("partial stdout:\n")
        sys.stderr.write(output)
        if not output.endswith("\n"):
            sys.stderr.write("\n")

    stderr = getattr(exc, "stderr", None)
    if stderr:
        sys.stderr.write("partial stderr:\n")
        sys.stderr.write(stderr)
        if not stderr.endswith("\n"):
            sys.stderr.write("\n")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--limit", type=int, default=10, help="Number of IDs to process"
    )
    parser.add_argument(
        "--input",
        dest="input_csv",
        type=Path,
        default=None,
        help=(
            "Path to the identifiers CSV forwarded to get_activity_data. "
            "Defaults to data/input/activity.csv when available."
        ),
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help=(
            "Forward the --dry-run flag to get_activity_data. "
            "Enable this option to perform a dry run; real writes are performed by default."
        ),
    )
    parser.add_argument(
        "--no-dry-run",
        dest="dry_run",
        action="store_false",
        help="Explicitly disable dry-run mode (default).",
    )
    parser.add_argument(
        "--timeout",
        type=float,
        default=600.0,
        help="Timeout in seconds for each pipeline run (default: 600).",
    )
    parser.add_argument(
        "--offline",
        action="store_true",
        help=(
            "Run the determinism check in offline mode using prepared fixtures. "
            "This replaces the get_activity_data CLI with a lightweight stub and "
            "avoids network calls."
        ),
    )
    parser.add_argument(
        "--fixtures-dir",
        type=Path,
        default=None,
        help=(
            "Path to a directory with offline fixtures. When --offline is set and "
            "no path is provided, defaults to tests/resources/expected_get_data."
        ),
    )
    parser.set_defaults(dry_run=False)

    args = parser.parse_args(argv)

    ensure_project_root(__file__)

    tmp_dir = Path(tempfile.mkdtemp(prefix="determinism_check_"))
    try:
        first = tmp_dir / "first.csv"
        second = tmp_dir / "second.csv"

        if args.offline and args.dry_run:
            sys.stderr.write(
                "Offline mode does not support --dry-run; rerun without --dry-run.\n"
            )
            return 2

        repo_root = Path(__file__).resolve().parents[1]
        fixtures_dir = args.fixtures_dir
        if args.offline and fixtures_dir is None:
            fixtures_dir = repo_root / "tests" / "resources" / "expected_get_data"
        if args.offline:
            if fixtures_dir is None or not fixtures_dir.exists():
                sys.stderr.write(
                    "Offline mode requested but fixtures directory is missing.\n"
                )
                return 2

        if args.input_csv is not None:
            input_csv = args.input_csv
            if not input_csv.exists():
                sys.stderr.write(
                    f"Input CSV '{input_csv}' does not exist; determinism check aborted.\n"
                )
                return 1
        else:
            input_csv = _default_input_csv(tmp_dir)

        try:
            first_run = _run_activity(
                args.limit,
                first,
                input_csv,
                dry_run=args.dry_run,
                timeout=args.timeout,
                offline=args.offline,
                fixtures_dir=fixtures_dir,
            )
        except subprocess.TimeoutExpired as exc:
            _report_process_timeout("first run", exc)
            return 1
        if first_run.returncode != 0:
            _report_process_failure("first run", first_run)
            return 1

        try:
            second_run = _run_activity(
                args.limit,
                second,
                input_csv,
                dry_run=args.dry_run,
                timeout=args.timeout,
                offline=args.offline,
                fixtures_dir=fixtures_dir,
            )
        except subprocess.TimeoutExpired as exc:
            _report_process_timeout("second run", exc)
            return 1
        if second_run.returncode != 0:
            _report_process_failure("second run", second_run)
            return 1

        if not first.exists() or not second.exists():
            if args.dry_run:
                first_logs_hash = _hash_process_output(first_run)
                second_logs_hash = _hash_process_output(second_run)

                if first_logs_hash != second_logs_hash:
                    print("Dry-run log hash check: mismatch")
                    print(
                        f"  first stdout/stderr SHA256:  {first_logs_hash}"
                    )
                    print(
                        f"  second stdout/stderr SHA256: {second_logs_hash}"
                    )
                    print(
                        "Dry-run outputs diverged; inspect the captured logs for"
                        " non-deterministic behaviour."
                    )
                    return 1

                print("Dry-run log hash check: matched")
                print(f"stdout/stderr SHA256: {first_logs_hash}")
                print("Deterministic dry-run output confirmed")
                return 0

            sys.stderr.write(
                "Determinism check failed: the pipeline exited without writing output files.\n"
            )
            sys.stderr.write(
                "Inspect the pipeline logs for warnings or errors before rerunning this check.\n"
            )
            return 1

        first_hash = _hash_file(first)
        second_hash = _hash_file(second)

        if first_hash != second_hash:
            print("Mismatch detected:")
            print(f"  first:  {first_hash}")
            print(f"  second: {second_hash}")
            return 1

        first_meta = _metadata_path(first)
        second_meta = _metadata_path(second)

        if first_meta.exists() and second_meta.exists():
            first_meta_hash = _hash_file(first_meta)
            second_meta_hash = _hash_file(second_meta)

            if first_meta_hash != second_meta_hash:
                print("Metadata hash check: mismatch")
                print("WARNING: Metadata hashes differ")
                print(f"  first:  {first_meta} ({first_meta_hash})")
                print(f"  second: {second_meta} ({second_meta_hash})")
                return 1

            print("Metadata hash check: matched")
            print(f"  first metadata SHA256:  {first_meta_hash}")
            print(f"  second metadata SHA256: {second_meta_hash}")
        else:
            print(
                "Metadata hash check: skipped (sidecars missing for one or both runs)"
            )

        print("Deterministic output confirmed")
        print(f"SHA256: {first_hash}")
        return 0
    finally:
        shutil.rmtree(tmp_dir, ignore_errors=True)


if __name__ == "__main__":
    raise SystemExit(main())
