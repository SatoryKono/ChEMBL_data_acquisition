#!/usr/bin/env python3
"""Best-effort determinism smoke test for data acquisition pipelines."""

from __future__ import annotations

if __package__ in {None, ""}:
    from _bootstrap import bootstrap_cli, ensure_project_root
else:  # pragma: no cover - executed when imported as a package module
    from ._bootstrap import bootstrap_cli, ensure_project_root

bootstrap_cli(__package__, __file__)
del bootstrap_cli

import argparse
import csv
import hashlib
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


def _hash_file(path: Path) -> str:
    hasher = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            hasher.update(chunk)
    return hasher.hexdigest()


def _metadata_path(csv_path: Path) -> Path:
    return csv_path.with_suffix(csv_path.suffix + ".meta.yaml")


def _run_activity(
    limit: int, destination: Path, input_csv: Path, *, dry_run: bool
) -> subprocess.CompletedProcess[str]:
    env = os.environ.copy()
    env.setdefault("PYTHONHASHSEED", "0")
    repo_root = Path(__file__).resolve().parents[1]
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


def _report_process_failure(label: str, result: subprocess.CompletedProcess[str]) -> None:
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
    parser.set_defaults(dry_run=False)

    args = parser.parse_args(argv)

    ensure_project_root(__file__)

    tmp_dir = Path(tempfile.mkdtemp(prefix="determinism_check_"))
    try:
        first = tmp_dir / "first.csv"
        second = tmp_dir / "second.csv"

        if args.input_csv is not None:
            input_csv = args.input_csv
            if not input_csv.exists():
                sys.stderr.write(
                    f"Input CSV '{input_csv}' does not exist; determinism check aborted.\n"
                )
                return 1
        else:
            input_csv = _default_input_csv(tmp_dir)

        first_run = _run_activity(args.limit, first, input_csv, dry_run=args.dry_run)
        if first_run.returncode != 0:
            _report_process_failure("first run", first_run)
            return 1

        second_run = _run_activity(args.limit, second, input_csv, dry_run=args.dry_run)
        if second_run.returncode != 0:
            _report_process_failure("second run", second_run)
            return 1

        if not first.exists() or not second.exists():
            if args.dry_run:
                sys.stderr.write(
                    "Determinism check failed: --dry-run prevents creating output files.\n"
                )
                sys.stderr.write(
                    "Re-run with --no-dry-run to verify that the pipeline produces stable results.\n"
                )
                return 2

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
