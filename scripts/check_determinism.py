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

from typing import Any

import yaml


def _hash_file(path: Path) -> str:
    hasher = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            hasher.update(chunk)
    return hasher.hexdigest()


_NON_DETERMINISTIC_META_KEYS = {"generated_at"}


def _metadata_path(csv_path: Path) -> Path:
    return csv_path.with_name(csv_path.name + ".meta.yaml")


def _load_metadata(path: Path) -> dict[str, Any] | None:
    if not path.exists():
        return None

    try:
        with path.open("r", encoding="utf-8") as handle:
            loaded = yaml.safe_load(handle)
    except (OSError, yaml.YAMLError):
        return {}

    if loaded is None:
        return {}

    if isinstance(loaded, dict):
        return dict(loaded)

    return {}


def _normalise_metadata(value: Any) -> Any:
    if isinstance(value, dict):
        items = []
        for key in sorted(value):
            if key in _NON_DETERMINISTIC_META_KEYS:
                continue
            items.append((key, _normalise_metadata(value[key])))
        return tuple(items)
    if isinstance(value, list):
        return tuple(_normalise_metadata(item) for item in value)
    return value


def _compare_metadata(first_csv: Path, second_csv: Path) -> bool:
    first_meta_path = _metadata_path(first_csv)
    second_meta_path = _metadata_path(second_csv)

    first_meta = _load_metadata(first_meta_path)
    second_meta = _load_metadata(second_meta_path)

    if first_meta is None and second_meta is None:
        return True

    if first_meta is None or second_meta is None:
        print("Metadata sidecar missing for one of the runs:")
        print(f"  first:  {first_meta_path}")
        print(f"  second: {second_meta_path}")
        return False

    return _normalise_metadata(first_meta) == _normalise_metadata(second_meta)


def _run_activity(
    limit: int, destination: Path, input_csv: Path, *, dry_run: bool
) -> subprocess.CompletedProcess[str]:
    env = os.environ.copy()
    env.setdefault("PYTHONHASHSEED", "0")
    cmd = [
        sys.executable,
        str(Path(__file__).resolve().parents[0] / "get_activity_data.py"),
        "--limit",
        str(limit),
        "--final-out",
        str(destination),
        "--input",
        str(input_csv),
    ]
    if dry_run:
        cmd.append("--dry-run")
    return subprocess.run(cmd, text=True, capture_output=True, env=env)


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
        action=argparse.BooleanOptionalAction,
        default=True,
        help=(
            "Forward the --dry-run flag to get_activity_data. "
            "Disable with --no-dry-run to perform a real write."
        ),
    )

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
            print(
                "Outputs not created; determinism check inconclusive", file=sys.stderr
            )
            return 1

        first_hash = _hash_file(first)
        second_hash = _hash_file(second)

        if first_hash != second_hash:
            print("Mismatch detected:")
            print(f"  first:  {first_hash}")
            print(f"  second: {second_hash}")
            return 1

        if not _compare_metadata(first, second):
            print("Metadata mismatch detected:")
            print(f"  first:  {_metadata_path(first)}")
            print(f"  second: {_metadata_path(second)}")
            return 1

        print("Deterministic output confirmed")
        print(f"SHA256: {first_hash}")
        return 0
    finally:
        shutil.rmtree(tmp_dir, ignore_errors=True)


if __name__ == "__main__":
    raise SystemExit(main())
