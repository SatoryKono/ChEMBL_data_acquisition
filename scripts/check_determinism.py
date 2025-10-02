"""Deterministic CSV regression checker.

The helper executes a configured CLI command twice, normalises the
produced CSV outputs and compares their cryptographic hashes. It is
intended for lightweight smoke verification in offline environments
where the default network-heavy pipelines cannot be exercised.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
import sys
import tempfile
import time
from pathlib import Path
from typing import Iterable

import shlex

import pandas as pd

DEFAULT_INPUT = Path("tests/data/input-smoke/testitem.csv")
DEFAULT_KEY_COLS = "molecule_chembl_id"
DEFAULT_COMMAND = (
    f"{{python}} -m library.utils.cli_tools.csv_utils_main "
    f"--input {{input}} --key-cols {{key_cols}} --output {{output}} "
    f"--log-level ERROR"
)


class DeterminismError(RuntimeError):
    """Raised when generated outputs diverge."""


def _build_command(
    template: Iterable[str] | None,
    output: Path,
    input_path: Path,
    key_cols: str,
) -> list[str]:
    """Materialise the command line for a single run."""

    placeholder_map = {
        "python": sys.executable,
        "output": str(output),
        "input": str(input_path),
        "key_cols": key_cols,
    }
    if not template:
        command = DEFAULT_COMMAND.format_map(placeholder_map)
        return shlex.split(command)

    if isinstance(template, str):
        command = template.format_map(placeholder_map)
        return shlex.split(command)

    return [token.format_map(placeholder_map) for token in template]


def _normalise_csv(path: Path, na_token: str) -> tuple[str, Path]:
    """Return a SHA256 hash of a canonicalised CSV file."""

    df = pd.read_csv(path, dtype="string").sort_index(axis=1)
    df = df.fillna(na_token)
    sort_cols = list(df.columns)
    if sort_cols:
        df = df.sort_values(by=sort_cols, kind="mergesort").reset_index(drop=True)
    csv_bytes = df.to_csv(index=False, na_rep=na_token).encode("utf-8")
    digest = hashlib.sha256(csv_bytes).hexdigest()
    normalised_path = path.with_suffix(".normalized.csv")
    normalised_path.write_bytes(csv_bytes)
    return digest, normalised_path


def run(argv: Iterable[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "command",
        nargs=argparse.REMAINDER,
        help=(
            "Command template executed per run. Placeholders {python}, {input}, "
            "{output} and {key_cols} are expanded. Defaults to a csv-utils job "
            "reading tests/data/input-smoke/testitem.csv."
        ),
    )
    parser.add_argument("--runs", type=int, default=2, help="How many executions to compare")
    parser.add_argument(
        "--input",
        type=Path,
        default=DEFAULT_INPUT,
        help="Input CSV for the default command template",
    )
    parser.add_argument(
        "--key-cols",
        default=DEFAULT_KEY_COLS,
        help="Key columns for the default command template",
    )
    parser.add_argument(
        "--keep-temp",
        action="store_true",
        help="Keep temporary outputs for inspection",
    )
    parser.add_argument(
        "--na-token",
        default="",
        help="String used to fill NA values before hashing",
    )
    args = parser.parse_args(argv)

    if args.runs < 2:
        parser.error("--runs must be at least 2")

    workspace = Path(tempfile.mkdtemp(prefix="determinism_"))
    results: list[dict[str, object]] = []
    try:
        for idx in range(args.runs):
            output = workspace / f"run_{idx}.csv"
            command = _build_command(args.command, output, args.input, args.key_cols)
            start = time.perf_counter()
            env = os.environ.copy()
            env.setdefault("PYTHONHASHSEED", "0")
            proc = subprocess.run(command, env=env, check=False, capture_output=True, text=True)
            elapsed = time.perf_counter() - start
            run_result = {
                "run": idx,
                "returncode": proc.returncode,
                "elapsed": elapsed,
                "stdout": proc.stdout.strip(),
                "stderr": proc.stderr.strip(),
                "output": str(output),
            }
            if proc.returncode != 0:
                raise DeterminismError(json.dumps(run_result, indent=2))
            digest, normalised = _normalise_csv(output, args.na_token)
            run_result.update({"sha256": digest, "normalised": str(normalised)})
            results.append(run_result)

        first_hash = results[0]["sha256"]
        divergent = [r for r in results if r["sha256"] != first_hash]
        summary = {"runs": results, "hash": first_hash}
        print(json.dumps(summary, indent=2))
        if divergent:
            raise DeterminismError(json.dumps(divergent, indent=2))
        return 0
    except DeterminismError as exc:
        print(f"Determinism check failed: {exc}", file=sys.stderr)
        return 1
    finally:
        if args.keep_temp:
            print(f"Temporary data preserved at {workspace}")
        else:
            for path in workspace.glob("*"):
                path.unlink(missing_ok=True)
            workspace.rmdir()


if __name__ == "__main__":  # pragma: no cover - CLI entry point
    raise SystemExit(run())
