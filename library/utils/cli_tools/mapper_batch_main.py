"""CLI for batch mapping ChEMBL IDs to UniProt accessions."""

from __future__ import annotations

import argparse
from collections.abc import Callable, Sequence
from pathlib import Path

import pandas as pd

from library import cli, io
from library.cli import (
    LoggerConfig,
    configure_logger,
)
from library.cli import build_parser as base_parser
from library.common.log import logger
from library.config import Config, ConfigError, ensure_dirs, print_config
from library.integration.mapper_batch_library import map_chembl_ids_to_uniprot


def run(cfg: Config, args: argparse.Namespace) -> int:
    """Run batch mapping of ChEMBL IDs to UniProt accessions."""
    try:
        marker_values = list(cfg.io.na_markers or ())
        keep_markers = bool(getattr(cfg.io, "keep_na_markers", False))
        na_values = marker_values if marker_values and not keep_markers else None
        try:
            df = io.read_csv(
                args.input_csv,
                cfg=cfg.io,
                sep=args.sep,
                encoding=args.encoding,
                dtype=str,
                na_values=na_values,
            )
        except io.CsvReadError as exc:
            logger.error(
                "read_fail",
                error=str(exc.original_error),
                path=str(args.input_csv),
                encoding=args.encoding,
                sep=args.sep,
            )
            return 1
        except (FileNotFoundError, OSError) as exc:
            logger.error(
                "read_fail",
                error=str(exc),
                path=str(args.input_csv),
                encoding=args.encoding,
                sep=args.sep,
            )
            return 1
        if args.column not in df.columns:
            logger.error(
                "missing_column",
                column=args.column,
                path=str(args.input_csv),
            )
            return 1

        marker_set = set(marker_values)
        ids = []
        for chembl_id in df[args.column]:
            if pd.isna(chembl_id):
                continue
            text = str(chembl_id).strip()
            if not text:
                continue
            if not keep_markers and text in marker_set:
                continue
            ids.append(text)
        mappings = map_chembl_ids_to_uniprot(
            ids,
            cfg.uniprot_mapping,
            batch_size=args.chunk_size,
            rps=args.rps,
            max_workers=args.workers,
        )
        for cid, uid in mappings.items():
            if uid:
                logger.info("mapped", chembl_id=cid, uniprot_id=uid)
            else:
                logger.warning("uniprot_id_missing", chembl_id=cid)
        df["mapping_uniprot_id"] = df[args.column].map(mappings).fillna("")
        output = args.output_csv or io.default_output_path(
            args.input_csv,
            cfg.io,
            date=getattr(args, "date", None),
        )
        try:
            io.write_csv(
                df,
                output,
                cfg=cfg,
                sep=args.sep,
                encoding=args.encoding,
            )
            logger.info("write_done", path=str(output))
        except OSError as exc:
            logger.error(
                "write_fail",
                error=str(exc),
                path=str(output),
            )
            return 1
        return 0
    except Exception as exc:  # pragma: no cover - defensive
        logger.exception("run_fail", exc=exc)
        return 1


def build_parser() -> tuple[argparse.ArgumentParser, LoggerConfig]:
    """Create CLI argument parser."""
    parser, log_cfg = base_parser(
        "Batch map ChEMBL IDs to UniProt accessions",
        column="chembl_id",
        chunk_size=10,
    )
    parser.add_argument(
        "--rps", type=float, default=5.0, help="Max requests per second"
    )
    parser.add_argument(
        "--workers", type=int, default=4, help="Number of worker threads"
    )
    parser.set_defaults(func=run)
    return parser, log_cfg


def main(argv: Sequence[str] | None = None) -> int:
    """Entry point for the batch mapper script.

    Notes
    -----
    Relative paths honour ``--base-path``, ``--input-dir`` and ``--output-dir``.
    """
    parser, log_cfg = build_parser()
    args = parser.parse_args(argv)
    input_path = getattr(args, "input_csv", None)
    output_stem = Path(input_path).stem if input_path else None
    cli.prepare_io_paths(args, output_stem=output_stem)
    run_id_value = getattr(args, "run_id", None)
    if isinstance(run_id_value, str):
        run_id_value = run_id_value.strip() or None
    if run_id_value is not None:
        log_cfg.run_id = run_id_value
    log_cfg.level = args.log_level
    logger = configure_logger(log_cfg)
    logger.info("pipeline_start", run_id=log_cfg.run_id)
    try:
        cfg = cli.apply_config_overrides(args, parser, args.config)
    except (ConfigError, FileNotFoundError, ValueError) as exc:
        logger.error(
            "config_error",
            error=str(exc),
            config=str(args.config),
        )
        logger.info("pipeline_end", exit_code=1)
        return 1

    try:
        ensure_dirs(cfg)
    except (ValueError, TypeError) as exc:
        logger.error(
            "config_error",
            error=str(exc),
            config=str(args.config),
        )
        logger.info("pipeline_end", exit_code=1)
        return 1
    except (FileNotFoundError, NotADirectoryError) as exc:
        logger.error("directory_setup_failed", error=str(exc))
        logger.info("pipeline_end", exit_code=1)
        return 1

    logger = configure_logger(log_cfg)
    if args.print_config:
        print_config(cfg)
        logger.info("pipeline_end", exit_code=0)
        return 0
    func: Callable[[Config, argparse.Namespace], int] = args.func
    exit_code = func(cfg, args)
    logger.info("pipeline_end", exit_code=exit_code)
    return exit_code


if __name__ == "__main__":  # pragma: no cover - CLI entry point
    raise SystemExit(main())
