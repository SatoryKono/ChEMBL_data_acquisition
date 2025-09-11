#!/usr/bin/env python3
"""Verify deterministic CSV output.

This script writes a small test :class:`pandas.DataFrame` twice using
:func:`library.csv_utils.write_csv_deterministic` and compares the
SHA-256 hashes of the resulting files. A non-zero exit code is returned if the
hashes differ.
"""

from __future__ import annotations

import sys
from pathlib import Path
import argparse
from tempfile import TemporaryDirectory

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

try:
    import pandas as pd
except ImportError as exc:  # pragma: no cover - import-time check
    raise SystemExit(
        "pandas is required to run scripts/check_determinism.py."
        " Install it with 'pip install pandas'."
    ) from exc


sys.path.append(str(Path(__file__).resolve().parents[1]))
from library.csv_utils import sha256_file, write_csv_deterministic  # noqa: E402
from library.log import logger  # noqa: E402
from library.logging_setup import LoggerConfig, configure_logger  # noqa: E402


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

    parser = argparse.ArgumentParser(description="Check deterministic CSV writing")
    parser.add_argument(
        "--log-level", default="INFO", help="Logging level (default: INFO)."
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


if __name__ == "__main__":
    raise SystemExit(main())
