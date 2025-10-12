"""CLI wrapper for :func:`write_csv_deterministic`.

This script reads an input CSV file and re-serialises it deterministically
using :func:`library.common.csv_utils.write_csv_deterministic`.

The command understands ``--base-path``, ``--input-dir`` and ``--output-dir``
for resolving relative paths. When ``--final-out`` is omitted, a file named
``output.<input-stem>_<YYYYMMDD>.csv`` is created next to the resolved input or
inside the selected output directory.
"""

from __future__ import annotations

# ruff: noqa: E402
import time
from collections.abc import Sequence
from datetime import date
from pathlib import Path
from typing import Any

import pandas as pd

from library import cli
from library.cli import LoggerConfig, configure_logger, create_logger_config
from library.cli_utils import build_parser, ensure_run_id
from library.common.csv_utils import write_csv_chunks_deterministic
from library.common.log import logger
from library.config import ConfigError, ensure_dirs, print_config
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
    provided_output = getattr(ns, "output_csv", None)
    provided_final = getattr(ns, "final_out", None)
    provided_date = getattr(ns, "date", None)
    cli.prepare_io_paths(ns, output_stem=output_stem)
    if output_stem is not None and provided_output is None and provided_final is None:
        target_dir = ns.output_dir or ns.base_path or Path(ns.input_csv).parent
        date_value = getattr(ns, "date", None)
        if provided_date is None:
            date_value = f"{date.today():%Y%m%d}"
            ns.date = date_value
        default_output = (
            target_dir / f"output.{output_stem}_{date_value}.csv"
        ).resolve()
        ns.output_csv = default_output
        ns.final_out = default_output
    try:
        cfg = cli.apply_config_overrides(
            ns,
            parser,
            ns.config,
            mapping={"chunk_size": "io.csv_chunksize"},
        )
    except (ConfigError, FileNotFoundError, ValueError) as exc:
        logger.error(
            "config_error",
            error=str(exc),
            config=str(ns.config),
        )
        return 1
    run_id_value = getattr(ns, "run_id", None)
    if isinstance(run_id_value, str):
        run_id_value = run_id_value.strip() or None
    log_cfg: LoggerConfig = create_logger_config(ns.log_level, run_id=run_id_value)
    ensure_run_id(ns, parser, log_cfg)
    if getattr(ns, "print_config", False):
        print_config(cfg)
        configure_logger(log_cfg)
        return 0
    arg_data = {field: getattr(ns, field) for field in CSVExportArgs.model_fields}
    args = CSVExportArgs.model_validate(arg_data)
    if not args.key_cols:
        parser.error("--key-cols must be provided")
    if args.log_level != log_cfg.level:
        log_cfg = LoggerConfig(
            level=args.log_level,
            run_id=log_cfg.run_id,
            redact_secrets=log_cfg.redact_secrets,
            stream=log_cfg.stream,
        )
    configure_logger(log_cfg)

    try:
        ensure_dirs(cfg)
    except (ValueError, TypeError) as exc:
        logger.error(
            "config_error",
            error=str(exc),
            config=str(ns.config),
        )
        return 1
    except OSError as exc:
        payload: dict[str, Any] = {"error": str(exc)}
        path = getattr(exc, "filename", None)
        if path:
            payload["path"] = str(path)
        logger.error("io_policy_violation", **payload)
        return 1

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
    csv_dtype = getattr(cfg.io, "csv_dtype", str)
    with pd.read_csv(
        args.input_csv,
        sep=sep,
        encoding=encoding,
        chunksize=chunk_size,
        dtype=csv_dtype,
    ) as reader:
        write_csv_chunks_deterministic(
            reader,
            output,
            col_order=args.col_order or None,
            key_cols=args.key_cols,
            chunksize=chunk_size,
            merge_chunksize=merge_chunk_size,
            sep=sep,
            encoding=encoding,
            cfg=cfg,
            drop_unexpected_cols=True,
        )
    elapsed = time.perf_counter() - start
    logger.info("write_done", path=str(args.output_csv))
    logger.info("run_completed", elapsed=elapsed)
    return 0


if __name__ == "__main__":  # pragma: no cover - CLI entry point
    raise SystemExit(main())
