"""CLI wrapper for :func:`write_csv_deterministic`.

This script reads an input CSV file and re-serialises it deterministically
using :func:`library.csv_utils.write_csv_deterministic`.

If ``--output`` is omitted, a file named
``output.<input-stem>_<YYYYMMDD>.csv`` is created next to the input.
"""

from __future__ import annotations

# ruff: noqa: E402
import time
from collections.abc import Sequence
from datetime import date
from pathlib import Path

import pandas as pd

from library import cli
from library.cli import LoggerConfig, configure_logger
from library.cli_utils import build_parser
from library.csv_utils import write_csv_chunks_deterministic
from library.log import logger
from library.parser_schema import CSVExportArgs


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
    ns = parser.parse_args(argv)
    if ns.output_csv is None:
        ns.output_csv = Path(ns.input_csv).with_name(
            f"output.{Path(ns.input_csv).stem}_{date.today():%Y%m%d}.csv"
        )
    cfg = cli.apply_config_overrides(
        ns,
        parser,
        ns.config,
        mapping={"chunk_size": "io.csv_chunksize"},
    )
    arg_data = {field: getattr(ns, field) for field in CSVExportArgs.model_fields}
    args = CSVExportArgs.model_validate(arg_data)
    if not args.key_cols:
        parser.error("--key-cols must be provided")
    configure_logger(LoggerConfig(level=args.log_level))

    start = time.perf_counter()

    sep = args.sep or cfg.io.csv_sep
    encoding = args.encoding or cfg.io.csv_encoding
    chunk_size = args.chunk_size or cfg.io.csv_chunksize

    output = args.output_csv or Path(args.input_csv).with_name(
        f"output.{Path(args.input_csv).stem}.csv"
    )
    with pd.read_csv(
        args.input_csv,
        sep=sep,
        encoding=encoding,
        chunksize=chunk_size,
    ) as reader:
        write_csv_chunks_deterministic(
            reader,
            output,
            col_order=args.col_order or None,
            key_cols=args.key_cols,
            chunksize=chunk_size,
            merge_chunksize=args.merge_chunk_size,
            sep=cfg.io.csv_sep,
            encoding=cfg.io.csv_encoding,
            cfg=cfg,
            drop_unexpected_cols=True,
        )
    elapsed = time.perf_counter() - start
    logger.info("write_done", path=str(args.output_csv))
    logger.info("run_completed", elapsed=elapsed)
    return 0


if __name__ == "__main__":  # pragma: no cover - CLI entry point
    raise SystemExit(main())
