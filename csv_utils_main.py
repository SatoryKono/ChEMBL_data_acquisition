"""CLI wrapper for :func:`write_csv_deterministic`.

This script reads an input CSV file and re-serialises it deterministically
using :func:`library.csv_utils.write_csv_deterministic`.
"""

from __future__ import annotations

import time
from pathlib import Path
from collections.abc import Sequence

import pandas as pd

from library.csv_utils import write_csv_deterministic
from library.cli_utils import build_parser
from library.log import logger
from library.cli import LoggerConfig, configure_logger


def main(argv: Sequence[str] | None = None) -> int:
    """CLI entry point.

    Parameters
    ----------
    argv:
        Optional sequence of command-line arguments.

    Returns
    -------
    int
        Zero on success.
    """

    parser = build_parser()
    args = parser.parse_args(argv)
    configure_logger(LoggerConfig(level=args.log_level))

    start = time.perf_counter()
    df = pd.read_csv(args.input_csv, sep=args.sep, encoding=args.encoding)
    output = args.output_csv or Path(args.input_csv).with_name(
        f"output_{Path(args.input_csv).stem}.csv"
    )
    write_csv_deterministic(
        df,
        output,
        col_order=args.col_order or None,
        key_cols=args.key_cols or None,
    )
    elapsed = time.perf_counter() - start
    logger.info("written %s", output)
    logger.info("completed in %.3f seconds", elapsed)
    return 0


if __name__ == "__main__":  # pragma: no cover - CLI entry point
    raise SystemExit(main())
