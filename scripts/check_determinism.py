#!/usr/bin/env python3
"""Verify deterministic CSV output.

This script writes a small test :class:`pandas.DataFrame` twice using
:func:`library.csv_utils.write_csv_deterministic` and compares the
SHA-256 hashes of the resulting files. A non-zero exit code is returned if the
hashes differ.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from tempfile import TemporaryDirectory
from time import perf_counter

# Allow running the script directly via ``python scripts/check_determinism.py``
# by ensuring the repository root is on ``sys.path`` when the module is executed
# outside of the ``scripts`` package. This mirrors the behaviour of installing
# the project in editable mode.
if __package__ in {None, ""}:
    sys.path.append(str(Path(__file__).resolve().parents[1]))

try:
    import pandas as pd
except ImportError as exc:  # pragma: no cover - import-time check
    raise SystemExit(
        "pandas is required to run scripts/check_determinism.py."
        " Install it with 'pip install pandas'."
    ) from exc

from library.cli import LoggerConfig, configure_logger
from library.csv_utils import sha256_file, write_csv_deterministic
from library.log import logger
from library.timing import log_duration


def run_check(tmp_dir: Path) -> bool:
    """Write the same CSV twice and compare their SHA-256 hashes.

    Parameters
    ----------
    tmp_dir:
        Directory where temporary CSV files will be created.

    Returns
    -------
    bool
        ``True`` if both files produce identical hashes, ``False`` otherwise.
    """

    df = pd.DataFrame({"a": [1, 2], "b": ["x", "y"]})

    first = tmp_dir / "first.csv"
    second = tmp_dir / "second.csv"

    # Write the DataFrame twice using the deterministic writer
    write_csv_deterministic(df, first)
    hash1 = sha256_file(first)

    write_csv_deterministic(df, second)
    hash2 = sha256_file(second)

    logger.debug("First hash: %s", hash1)
    logger.debug("Second hash: %s", hash2)

    return hash1 == hash2


def main() -> int:
    """CLI entry point.

    Returns
    -------
    int
        ``0`` on success, ``1`` when hashes differ.
    """

    start = perf_counter()
    try:
        parser = argparse.ArgumentParser(
            description="Check deterministic CSV writing",
        )
        parser.add_argument(
            "--log-level",
            default="INFO",
            help="Logging level (default: INFO).",
        )
        args = parser.parse_args()

        configure_logger(LoggerConfig(level=args.log_level))

        with TemporaryDirectory() as tmp:
            ok = run_check(Path(tmp))

        if ok:
            logger.info("Hashes match; output is deterministic")
            return 0

        logger.error("Hashes differ; output is not deterministic")
        return 1
    finally:
        log_duration(start)


if __name__ == "__main__":
    raise SystemExit(main())
