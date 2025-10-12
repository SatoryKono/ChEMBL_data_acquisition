"""CLI for combining ChEMBL initialisation Excel sources.

Exports pair tables without merging:

- ``pairs_same_document.csv`` from sheet ``pairs_same_doc`` in ``--same-doc``
- ``pairs_independent.csv`` and ``pairs_non_independent.csv`` derived from
  ``step5_pairs`` in ``--all-doc`` based on the ``INDEPENDENT`` flag.

Additionally, for each pair table the corresponding ``activity``, ``assay``,

``document``, ``target``, ``testitem`` and ``system`` entries are exported with
matching suffixes, for example ``activity_independent.csv`` or
``assay_same_document.csv``.

"""

from __future__ import annotations

import argparse
import sys
from collections.abc import Sequence
from pathlib import Path

if __package__ in {None, ""}:
    project_root = Path(__file__).resolve().parents[3]
    project_root_str = str(project_root)
    if project_root_str not in sys.path:
        sys.path.insert(0, project_root_str)

from library.utils import bootstrap  # noqa: E402

bootstrap.ensure_project_root()

from library import cli  # noqa: E402
from library.cli import (  # noqa: E402
    LoggerConfig,
    configure_logger,
    path_argument,
)
from library.cli import (  # noqa: E402
    build_parser as base_parser,
)
from library.cli_utils import ensure_run_id  # noqa: E402
from library.common.log import logger  # noqa: E402
from library.config import Config, ConfigError, ensure_dirs, print_config  # noqa: E402
from library.integration import input_initialisation_library as lib  # noqa: E402
from library.qa.reporting import build_table_quality_hook  # noqa: E402


def run(cfg: Config, args: argparse.Namespace) -> int:
    """Execute table combination routine.

    Parameters
    ----------
    cfg : Config
        Application configuration.
    args:
        Parsed command line arguments.

    Returns
    -------
    int
        Zero on success, non-zero otherwise.

    """
    try:
        # ``same_doc`` and ``all_doc`` may be ``None`` when the caller bypasses
        # :func:`main` and invokes :func:`run` directly.  Guard against this to
        # avoid ``AttributeError`` when accessing ``Path.exists``.
        if args.same_doc is None or not args.same_doc.exists():
            raise FileNotFoundError(f"{args.same_doc} does not exist")
        if args.all_doc is None or not args.all_doc.exists():
            raise FileNotFoundError(f"{args.all_doc} does not exist")

        out_dir = Path(args.out_dir or cfg.init.output_dir)
        out_dir.mkdir(parents=True, exist_ok=True)

        logger.info("load_workbooks")
        same = lib.load_same_doc(args.same_doc)
        all_ = lib.load_all_doc(args.all_doc)

        logger.info("combine_tables")
        tables = lib.build_combined_tables(
            same,
            all_,
            dictionary_dir=args.dictionary_dir,
            targets_type_csv=cfg.resources.targets_type_csv,
        )
        logger.info("generate_pair_tables")
        tables = lib.generate_pair_entity_tables(
            tables,
            {
                "pairs_independent": "independent",
                "pairs_non_independent": "non_independent",
                "pairs_same_document": "same_document",
            },
        )
        logger.info("save_output")
        paths = lib.save_tables(tables, out_dir, cfg, fmt=args.format)
        # Ensure that files were actually written to disk
        missing = [str(p) for p in paths.values() if not p.exists()]
        if missing:
            raise RuntimeError("failed to write output files: " + ", ".join(missing))

        logger.info("generate_quality_reports")
        report_dir = out_dir / "data_validity_report"
        report_dir.mkdir(parents=True, exist_ok=True)
        doc_quality_cfg = cfg.system.doc_quality
        for entity, path in paths.items():
            logger.info("profiling", entity=entity)
            table_quality = build_table_quality_hook(
                doc_quality_cfg,
                table_name=path.with_suffix(""),
                destination=report_dir,
            )
            table_quality(path)

        logger.info("save_done", tables=len(paths), path=str(out_dir))
        return 0
    except KeyError as exc:
        msg = exc.args[0] if exc.args else str(exc)
        if any(
            indicator in msg
            for indicator in (
                "not in index",
                "table missing expected columns",
                "required column",
            )
        ):
            # Log a column-specific message when the error clearly refers to
            # missing columns rather than an absent table.
            logger.error("missing_columns", details=msg)
        else:
            logger.error("missing_table", table=msg)
        return 1


def build_parser() -> tuple[argparse.ArgumentParser, LoggerConfig]:
    """Create argument parser."""
    parser, log_cfg = base_parser(__doc__ or "Input initialisation", column="chembl_id")
    parser.add_argument(
        "--same-doc",
        type=path_argument,
        help="Path to same document workbook (default: config init.same_doc)",
    )
    parser.add_argument(
        "--all-doc",
        type=path_argument,
        help="Path to all document workbook (default: config init.all_doc)",
    )
    parser.add_argument(
        "--dictionary-dir",
        type=path_argument,
        default=None,
        help=(
            "Directory with _target/targets_type.csv and _Curation/citation_fraction.csv "
            "(default: config resources.dictionary_dir)"
        ),
    )
    parser.add_argument(
        "--out-dir",
        type=path_argument,
        help="Output directory (default: config init.output_dir)",
    )
    parser.add_argument(
        "--format", choices=["csv"], default="csv", help="Output format"
    )
    parser.set_defaults(func=run)
    return parser, log_cfg


def main(argv: Sequence[str] | None = None) -> int:
    """Entry point using :class:`Config` for defaults.

    Notes
    -----
    Relative paths honour ``--base-path``, ``--input-dir`` and ``--output-dir``.
    """
    parser, log_cfg = build_parser()

    args = parser.parse_args(argv)
    input_path = getattr(args, "input_csv", None)
    output_stem = Path(input_path).stem if input_path else None
    cli.prepare_io_paths(args, output_stem=output_stem)
    ensure_run_id(args, parser, log_cfg)
    log_cfg.level = args.log_level
    logger = configure_logger(log_cfg)
    logger.info("pipeline_start", run_id=log_cfg.run_id)
    try:
        cfg: Config = cli.apply_config_overrides(
            args,
            parser,
            args.config,
            mapping={
                "same_doc": "init.same_doc",
                "all_doc": "init.all_doc",
                "out_dir": "init.output_dir",
                "dictionary_dir": "resources.dictionary_dir",
            },
        )
    except (ConfigError, FileNotFoundError, ValueError) as exc:
        logger.error(
            "config_error",
            error=str(exc),
            config=str(args.config),
        )
        logger.info("pipeline_fail", run_id=log_cfg.run_id)
        return 1

    try:
        if args.print_config:
            print_config(cfg)
            configure_logger(log_cfg)
            logger.info("pipeline_done", run_id=log_cfg.run_id)
            return 0
        ensure_dirs(cfg)

        # ``dictionary_dir`` is optional and may be undefined in the
        # configuration.  Convert provided paths to :class:`Path` objects only
        # when values are available to avoid ``TypeError`` on ``None``.
        args.same_doc = Path(args.same_doc) if args.same_doc is not None else None
        args.all_doc = Path(args.all_doc) if args.all_doc is not None else None
        # ``out_dir`` is mandatory but may originate from the configuration
        # file.  Resolve to :class:`Path` lazily to avoid ``TypeError`` when both
        # sources are missing.
        args.out_dir = (
            Path(args.out_dir)
            if args.out_dir is not None
            else Path(cfg.init.output_dir)
        )
        args.dictionary_dir = (
            Path(args.dictionary_dir) if args.dictionary_dir is not None else None
        )
        if args.same_doc is None or args.all_doc is None:
            raise ValueError("same_doc and all_doc paths must be provided")

        logger = configure_logger(log_cfg)
    except (ValueError, TypeError) as exc:
        logger.error(
            "config_error",
            error=str(exc),
            config=str(args.config),
        )
        logger.info("pipeline_fail", run_id=log_cfg.run_id)
        return 1
    except (FileNotFoundError, NotADirectoryError) as exc:
        logger.error(
            "directory_setup_failed",
            error=str(exc),
            output=str(args.out_dir),
        )
        logger.info("pipeline_fail", run_id=log_cfg.run_id)
        return 1
    exit_code = int(args.func(cfg, args))
    if exit_code == 0:
        logger.info("pipeline_done", run_id=log_cfg.run_id)
    else:
        logger.info("pipeline_fail", run_id=log_cfg.run_id)
    return exit_code


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
