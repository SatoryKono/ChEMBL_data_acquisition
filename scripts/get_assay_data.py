"""Command line interface for retrieving ChEMBL assay data."""

from __future__ import annotations

import argparse
import sys
from collections.abc import Sequence
from pathlib import Path

import requests
from pandera.errors import SchemaErrors

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from library import assay_postprocessing as ap  # noqa: E402
from library import chembl_library as cl  # noqa: E402
from library import io  # noqa: E402
from library.chembl_client import ChemblClient  # noqa: E402
from library.cli import (  # noqa: E402
    LoggerConfig,
    apply_config_overrides,
    configure_logger,
)
from library.cli import (  # noqa: E402
    build_parser as base_parser,
)
from library.config import (  # noqa: E402
    Config,
    _serialize_paths,
    ensure_dirs,
    print_config,
)
from library.log import logger  # noqa: E402
from library.metadata import Stats, file_sha256, write_meta_yaml  # noqa: E402
from library.sidecar import SidecarErrors  # noqa: E402
from library.table_quality import analyze_table_quality  # noqa: E402
from schemas import AssaysSchema, normalize_assays  # noqa: E402


def run_chembl(cfg: Config, args: argparse.Namespace) -> int:
    """Execute assay retrieval from the ChEMBL API.

    Parameters
    ----------
    cfg : Config
        Application configuration.
    args : argparse.Namespace
        Parsed command-line arguments.

    Returns
    -------
    int
        Zero on success, non-zero on failure.

    """
    # Prepare HTTP session for ChEMBL requests
    client = ChemblClient(cfg.api, cfg.retry, cfg.chembl)

    try:
        ids = io.read_ids(args.input_csv, column=cfg.assay.column, cfg=cfg.io)
    except (FileNotFoundError, ValueError) as exc:
        logger.error("%s", exc)
        return 1

    try:
        df = cl.get_assays(
            ids,
            cfg=cfg.api,
            client=client,
            chunk_size=cfg.assay.chunk_size,
            timeout=cfg.assay.timeout,
        )
    except (requests.RequestException, ValueError) as exc:
        logger.error("failed to retrieve assays: %s", exc)
        return 1
    df = ap.postprocess_assays(df)
    output = args.output_csv or io.default_output_path(args.input_csv, cfg.io)
    df = normalize_assays(df)
    rows_total = len(df)
    exit_code = 0
    required_cols = set(AssaysSchema.columns.keys())
    if required_cols.issubset(df.columns):
        try:
            df = AssaysSchema.validate(df, lazy=True)
        except SchemaErrors as exc:
            failure_path = output.with_name(f"{output.stem}_failure_cases.csv")
            errors = SidecarErrors()
            for row in exc.failure_cases.to_dict("records"):
                errors.add_error(row)
            errors.save(failure_path)
            logger.error(
                "validation failed; wrote %d failure cases to %s",
                len(exc.failure_cases),
                failure_path,
            )
            df = getattr(exc, "validated_data", df)
            exit_code = 1
    else:
        missing = required_cols.difference(df.columns)
        logger.warning("Skipping validation due to missing columns: %s", missing)
    rows_kept = len(df)
    rows_dropped = rows_total - rows_kept
    try:
        key_cols = [c for c in ["assay_chembl_id"] if c in df.columns]
        csv_path = io.write_csv(
            df,
            output,
            cfg=cfg,
            key_cols=key_cols or None,
        )
        logger.info("Wrote %d rows to %s", rows_kept, csv_path)
    except OSError as exc:
        logger.error("failed to write output CSV: %s", exc)
        return 1

    stats: Stats = {
        "rows_total": rows_total,
        "rows_kept": rows_kept,
        "rows_dropped": rows_dropped,
        "output_sha256": file_sha256(csv_path),
    }
    write_meta_yaml(
        csv_path=csv_path,
        command=" ".join(sys.argv),
        config_subset=_serialize_paths(cfg.to_dict()),
        inputs={"input_csv": str(args.input_csv)},
        stats=stats,
        schema="AssaysSchema",
    )

    try:
        analyze_table_quality(df, table_name=str(output.with_suffix("")))
    except ValueError as exc:
        logger.error("failed to generate quality report: %s", exc)
        return 1
    return exit_code


def build_parser() -> tuple[argparse.ArgumentParser, LoggerConfig]:
    """Create the command-line argument parser."""
    parser, log_cfg = base_parser(
        "ChEMBL assay data utilities", column="assay_chembl_id", chunk_size=10
    )
    parser.add_argument(
        "--timeout",
        type=float,
        default=30.0,
        help="Timeout in seconds for each HTTP request",
    )
    parser.set_defaults(func=run_chembl)
    return parser, log_cfg


def main(argv: Sequence[str] | None = None) -> int:
    """Command line entry point using :class:`Config` for defaults."""
    parser, log_cfg = build_parser()
    args = parser.parse_args(argv)
    log_cfg.level = args.log_level
    logger = configure_logger(log_cfg)
    logger.info("pipeline start run_id=%s", log_cfg.run_id, extra={"event": "start"})
    try:
        cfg: Config = apply_config_overrides(
            args,
            parser,
            args.config,
            mapping={
                "timeout": "assay.timeout",
                "column": "assay.column",
                "chunk_size": "assay.chunk_size",
            },
        )
        if args.print_config:
            print_config(cfg)
            configure_logger(log_cfg, fmt=cfg.log.format, datefmt=cfg.log.datefmt)
            logger.info(
                "pipeline done run_id=%s", log_cfg.run_id, extra={"event": "done"}
            )
            return 0
        ensure_dirs(cfg)
        logger = configure_logger(log_cfg, fmt=cfg.log.format, datefmt=cfg.log.datefmt)
    except (ValueError, TypeError) as exc:
        logger.error("%s", exc)
        logger.info("pipeline fail run_id=%s", log_cfg.run_id, extra={"event": "fail"})
        return 1
    except (FileNotFoundError, NotADirectoryError) as exc:
        logger.error("failed to set up directories: %s", exc)
        logger.info("pipeline fail run_id=%s", log_cfg.run_id, extra={"event": "fail"})
        return 1
    exit_code = args.func(cfg, args)
    if exit_code == 0:
        logger.info("pipeline done run_id=%s", log_cfg.run_id, extra={"event": "done"})
    else:
        logger.info("pipeline fail run_id=%s", log_cfg.run_id, extra={"event": "fail"})
    return exit_code


if __name__ == "__main__":  # pragma: no cover - CLI entry point
    raise SystemExit(main())
