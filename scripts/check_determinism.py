#!/usr/bin/env python3
"""Best-effort determinism smoke test for data acquisition pipelines."""

from __future__ import annotations

# ruff: noqa: E402

if __package__ in {None, ""}:
    from _bootstrap import bootstrap_cli, ensure_project_root
else:  # pragma: no cover - executed when imported as a package module
    from ._bootstrap import bootstrap_cli, ensure_project_root

bootstrap_cli(__package__, __file__)
del bootstrap_cli

import argparse
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


def _run_activity(limit: int, destination: Path) -> subprocess.CompletedProcess[str]:
    env = os.environ.copy()
    env.setdefault("PYTHONHASHSEED", "0")
    cmd = [
        sys.executable,
        str(Path(__file__).resolve().parents[0] / "get_activity_data.py"),
        "--limit",
        str(limit),
        "--final-out",
        str(destination),
    ]
    return subprocess.run(cmd, text=True, capture_output=True, env=env)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--limit", type=int, default=10, help="Number of IDs to process"
    )
    args = parser.parse_args(argv)

    ensure_project_root(__file__)

    tmp_dir = Path(tempfile.mkdtemp(prefix="determinism_check_"))
    try:
        first = tmp_dir / "first.csv"
        second = tmp_dir / "second.csv"

        first_run = _run_activity(args.limit, first)
        if first_run.returncode != 0:
            sys.stderr.write("first run failed\n")
            sys.stderr.write(first_run.stderr)
            return 1

        second_run = _run_activity(args.limit, second)
        if second_run.returncode != 0:
            sys.stderr.write("second run failed\n")
            sys.stderr.write(second_run.stderr)
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

        print("Deterministic output confirmed")
        print(f"SHA256: {first_hash}")
        return 0
    finally:
        shutil.rmtree(tmp_dir, ignore_errors=True)


if __name__ == "__main__":
    raise SystemExit(main())
