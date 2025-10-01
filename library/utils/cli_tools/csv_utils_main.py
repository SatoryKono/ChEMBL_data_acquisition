"""CLI wrapper for :func:`write_csv_deterministic`.

This script reads an input CSV file and re-serialises it deterministically
using :func:`library.csv_utils.write_csv_deterministic`.

The command understands ``--base-path``, ``--input-dir`` and ``--output-dir``
for resolving relative paths. When ``--output`` is omitted, a file named
``output.<input-stem>_<YYYYMMDD>.csv`` is created next to the resolved input or
inside the selected output directory.
"""

from __future__ import annotations

# ruff: noqa: E402
import time
from collections.abc import Sequence
from datetime import date
from pathlib import Path

import pandas as pd

from library import cli
from library.cli import configure_logger, create_logger_config
from library.cli_utils import build_parser
from library.csv_utils import write_csv_chunks_deterministic
from library.config import print_config
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
    input_path = getattr(ns, "input_csv", None)
    output_stem = Path(input_path).stem if input_path else None
    cli.prepare_io_paths(ns, output_stem=output_stem)
    if ns.output_csv is None and output_stem is not None:
        target_dir = ns.output_dir or ns.base_path or Path(ns.input_csv).parent
        date_value = getattr(ns, "date", None)
        if date_value is None:
            date_value = f"{date.today():%Y%m%d}"
            setattr(ns, "date", date_value)
        ns.output_csv = (target_dir / f"output.{output_stem}_{date_value}.csv").resolve()
    cfg = cli.apply_config_overrides(
        ns,
        parser,
        ns.config,
        mapping={"chunk_size": "io.csv_chunksize"},
    )
    if getattr(ns, "print_config", False):
        print_config(cfg)
        configure_logger(LoggerConfig(level=ns.log_level))
        return 0
    arg_data = {field: getattr(ns, field) for field in CSVExportArgs.model_fields}
    args = CSVExportArgs.model_validate(arg_data)
    if not args.key_cols:
        parser.error("--key-cols must be provided")
    log_cfg = create_logger_config(args.log_level)
    configure_logger(log_cfg)

    start = time.perf_counter()

    sep = args.sep or cfg.io.csv_sep
    encoding = args.encoding or cfg.io.csv_encoding
    chunk_size = args.chunk_size or cfg.io.csv_chunksize
    merge_chunk_size = args.merge_chunk_size
    if merge_chunk_size is None:
        merge_chunk_size = max(chunk_size, 1)

    output = args.output_csv or Path(args.input_csv).with_name(
        f"output.{Path(args.input_csv).stem}.csv"
    )
    ns.output_csv = output
    if args.output_csv != output:
        args = args.model_copy(update={"output_csv": output})
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
            merge_chunksize=merge_chunk_size,
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
