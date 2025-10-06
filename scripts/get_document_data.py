"""Command line interface for retrieving document metadata from external sources.

The tool integrates :mod:`library.integration.pubmed_library` and
:mod:`library.integration.chembl_library` to collect information about publications from
several public APIs.  The interface mirrors :mod:`scripts.get_target_data` and exposes a
single entry point configured via ``--mode``:

``--mode pubmed``
    Query PubMed, Semantic Scholar, OpenAlex and CrossRef for a list of PMIDs.
``--mode chembl``
    Retrieve document information from the ChEMBL API.
``--mode all``
    Run the ChEMBL and PubMed pipelines and merge the results.

Example
-------
Fetch PubMed metadata for identifiers listed in ``pmids.csv``::

    python scripts/get_document_data.py --mode pubmed --config config/config.yaml --input pmids.csv --final-out output.csv

The input file must contain a ``PMID`` column.

"""

from __future__ import annotations

if __package__ in {None, ""}:
    from _bootstrap import bootstrap_cli
else:  # pragma: no cover - executed when imported as a package module
    from ._bootstrap import bootstrap_cli

bootstrap_cli(__package__, __file__)
del bootstrap_cli

import argparse
import os
import sys
from collections.abc import Sequence
from pathlib import Path

from library.document_defaults import ALL_DEFAULTS, CHEMBL_DEFAULTS, PUBMED_DEFAULTS
from library.cli import (
    LoggerConfig,
    build_root_parser,
    configure_logger,
    path_argument,
    positive_int,
    prepare_io_paths,
)
from library.cli.logging import setup_cli_logging
from library.cli.utils import run_cli_command
from library.common.log import logger
from library.pipelines.document import runner as document_runner

DEFAULT_INPUT_NAME = "document.csv"
DEFAULT_OUTPUT_STEM = "documents"

run_pubmed = document_runner.run_pubmed
run_chembl = document_runner.run_chembl
run_all = document_runner.run_all
run = document_runner.run
MODE_HANDLERS = document_runner.MODE_HANDLERS


class _FallbackPathAction(argparse.Action):
    """Store fallback DOI CSV path under both legacy and new attribute names."""

    def __call__(
        self,
        parser: argparse.ArgumentParser,
        namespace: argparse.Namespace,
        values: object,
        option_string: str | None = None,
    ) -> None:
        setattr(namespace, self.dest, values)
        setattr(namespace, "fallback_doi_csv", values)


def build_parser() -> tuple[argparse.ArgumentParser, LoggerConfig]:
    """Create the argument parser for document utilities."""

    root, _, log_cfg = build_root_parser()
    root.set_defaults(input_csv=Path(DEFAULT_INPUT_NAME))
    parser = argparse.ArgumentParser(
        description="Document data utilities",
        parents=[root],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.set_defaults(input_csv=Path(DEFAULT_INPUT_NAME))

    pipeline_group = parser.add_argument_group("Pipeline selection")
    pipeline_group.add_argument(
        "--mode",
        choices=("chembl", "pubmed", "all"),
        required=False,
        default=None,
        help="Document pipeline to execute",
    )
    parser.add_argument(
        "command",
        nargs="?",
        choices=("chembl", "pubmed", "all"),
        help=argparse.SUPPRESS,
    )
    pipeline_group.add_argument(
        "--column",
        default=PUBMED_DEFAULTS.column,
        help=(
            "Input column containing identifiers (defaults: "
            f"pubmed={PUBMED_DEFAULTS.column}, chembl={CHEMBL_DEFAULTS.column}, "
            f"all={ALL_DEFAULTS.column})"
        ),
    )
    pipeline_group.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Maximum number of identifiers to process (default: no limit)",
    )
    pipeline_group.add_argument(
        "--offset",
        type=int,
        default=0,
        help="Number of identifiers to skip before processing",
    )
    pipeline_group.add_argument(
        "--openalex-rps",
        type=float,
        default=None,
        help="Requests per second limit for OpenAlex lookups",
    )
    pipeline_group.add_argument(
        "--crossref-rps",
        type=float,
        default=None,
        help="Requests per second limit for CrossRef lookups",
    )

    single_group = parser.add_argument_group("Single pipeline options")
    single_group.add_argument(
        "--batch-size",
        type=positive_int,
        default=PUBMED_DEFAULTS.batch_size,
        help=(
            "Maximum PMIDs per PubMed request when running in pubmed mode "
            f"(default: {PUBMED_DEFAULTS.batch_size})"
        ),
    )
    single_group.add_argument(
        "--sleep",
        type=float,
        default=PUBMED_DEFAULTS.sleep,
        help=(
            "Seconds to sleep between PubMed requests when running in pubmed mode "
            f"(default: {PUBMED_DEFAULTS.sleep})"
        ),
    )
    single_group.add_argument(
        "--workers",
        type=int,
        default=PUBMED_DEFAULTS.workers,
        help=(
            "Number of concurrent PubMed requests when running in pubmed mode "
            f"(default: {PUBMED_DEFAULTS.workers})"
        ),
    )
    single_group.add_argument(
        "--chunk-size",
        type=positive_int,
        default=CHEMBL_DEFAULTS.chunk_size,
        help=(
            "Maximum identifiers per ChEMBL request when running in chembl mode "
            f"(default: {CHEMBL_DEFAULTS.chunk_size})"
        ),
    )
    single_group.add_argument(
        "--timeout",
        type=float,
        default=CHEMBL_DEFAULTS.timeout,
        help=(
            "HTTP read timeout in seconds (defaults: chembl/all="
            f"{CHEMBL_DEFAULTS.timeout}, pubmed={PUBMED_DEFAULTS.timeout})"
        ),
    )

    combined_group = parser.add_argument_group("Combined pipeline overrides")
    combined_group.add_argument(
        "--chembl-chunk-size",
        "--chembl-batch-size",
        dest="chembl_chunk_size",
        type=positive_int,
        default=ALL_DEFAULTS.chunk_size,
        help=(
            "Maximum identifiers per ChEMBL request when running in all mode "
            f"(default: {ALL_DEFAULTS.chunk_size})"
        ),
    )
    combined_group.add_argument(
        "--chembl-timeout",
        dest="chembl_timeout",
        type=float,
        default=ALL_DEFAULTS.timeout,
        help=(
            "Timeout in seconds for ChEMBL requests when running in all mode "
            f"(default: {ALL_DEFAULTS.timeout})"
        ),
    )
    combined_group.add_argument(
        "--pubmed-sleep",
        "--chembl-sleep",
        dest="pubmed_sleep",
        type=float,
        default=ALL_DEFAULTS.sleep,
        help=(
            "Seconds to sleep between PubMed requests when running in all mode "
            f"(default: {ALL_DEFAULTS.sleep})"
        ),
    )
    combined_group.add_argument(
        "--pubmed-workers",
        "--chembl-workers",
        dest="pubmed_workers",
        type=int,
        default=ALL_DEFAULTS.workers,
        help=(
            "Number of concurrent PubMed requests when running in all mode "
            f"(default: {ALL_DEFAULTS.workers})"
        ),
    )
    combined_group.add_argument(
        "--pubmed-batch-size",
        "--pubmed-chunk-size",
        dest="pubmed_batch_size",
        type=positive_int,
        default=ALL_DEFAULTS.batch_size,
        help=(
            "Maximum PMIDs per PubMed request when running in all mode "
            f"(default: {ALL_DEFAULTS.batch_size})"
        ),
    )
    combined_group.add_argument(
        "--pubmed-timeout",
        dest="pubmed_timeout",
        type=float,
        default=PUBMED_DEFAULTS.timeout,
        help=(
            "Timeout in seconds for PubMed requests when running in all mode "
            f"(default: {PUBMED_DEFAULTS.timeout})"
        ),
    )

    fallback_group = parser.add_argument_group("Fallback DOI overrides")
    fallback_group.add_argument(
        "--fallback-doi-enabled",
        action="store_true",
        help="Enable lookup of DOI overrides from a CSV file",
    )
    parser.set_defaults(fallback_doi_csv=None)
    fallback_group.add_argument(
        "--fallback-doi-path",
        "--fallback-doi-csv",
        dest="fallback_doi_path",
        action=_FallbackPathAction,
        type=path_argument,
        default=None,
        help="CSV file containing DOI overrides keyed by PMID",
    )
    fallback_group.add_argument(
        "--fallback-doi-col-pmid",
        default="PMID",
        help="Column containing PubMed identifiers in the fallback CSV",
    )
    fallback_group.add_argument(
        "--fallback-doi-col-doi",
        default="DOI",
        help="Column containing DOI values in the fallback CSV",
    )
    fallback_group.add_argument(
        "--fallback-doi-delimiter",
        default=None,
        help="Delimiter used when reading the fallback CSV (default: io.csv_sep)",
    )
    fallback_group.add_argument(
        "--fallback-doi-encoding",
        default=None,
        help="Encoding used for the fallback CSV (default: io.csv_encoding)",
    )
    fallback_group.add_argument(
        "--fallback-doi-overwrite",
        action="store_true",
        help="Allow replacing existing DOIs with fallback values",
    )

    return parser, log_cfg


def main(argv: Sequence[str] | None = None) -> int:
    """Command line entry point using :class:`Config` for defaults.

    Parameters
    ----------
    argv : Sequence[str], optional
        Command line arguments to parse.  When ``None`` the arguments from
        :data:`sys.argv` are used instead.

    Returns
    -------
    int
        Exit status to return from :func:`sys.exit`.

    Raises
    ------
    SystemExit
        Raised when argument validation fails and ``argparse`` aborts
        execution.
    """

    parser, log_cfg = build_parser()
    if argv is None:
        argv_list: list[str] = list(sys.argv[1:])
    else:
        argv_list = [str(item) for item in argv]
    args = parser.parse_args(argv_list)
    prepare_io_paths(
        args,
        input_default=DEFAULT_INPUT_NAME,
        output_stem=DEFAULT_OUTPUT_STEM,
    )
    limit_value = getattr(args, "limit", None)
    if limit_value == 0:
        logger.info("pipeline_skip_limit", limit=limit_value)
        return 0
    mode = getattr(args, "mode", None)
    command_value = getattr(args, "command", None)
    if not mode and command_value:
        mode = command_value
        setattr(args, "mode", mode)
    if not mode:
        parser.error("--mode is required")
    if command_value is None and mode is not None:
        setattr(args, "command", mode)
    mode = str(mode)
    if limit_value is not None and limit_value < 0:
        parser.error("--limit must be zero or a positive integer")
    offset_value = getattr(args, "offset", 0)
    if offset_value < 0:
        parser.error("--offset must be zero or a positive integer")
    if mode in {"pubmed", "all"}:
        fallback_enabled_cli = getattr(args, "fallback_doi_enabled", False)
        fallback_path_cli = getattr(args, "fallback_doi_path", None)
        if fallback_path_cli is not None and not fallback_enabled_cli:
            fallback_enabled_cli = True
            setattr(args, "fallback_doi_enabled", True)
        if fallback_enabled_cli:
            if fallback_path_cli is None:
                parser.error(
                    "--fallback-doi-path is required when fallback DOI overrides are enabled"
                )
            fallback_path = (
                fallback_path_cli
                if isinstance(fallback_path_cli, Path)
                else Path(str(fallback_path_cli))
            )
            if not fallback_path.exists():
                parser.error("--fallback-doi-path must point to an existing file")
            if not fallback_path.is_file():
                parser.error("--fallback-doi-path must be a file")
            if not os.access(fallback_path, os.R_OK):
                parser.error("--fallback-doi-path must be readable")
            delimiter_cli = getattr(args, "fallback_doi_delimiter", None)
            if delimiter_cli is not None:
                delimiter_text = str(delimiter_cli)
                if not delimiter_text:
                    parser.error("--fallback-doi-delimiter must not be empty")
                if len(delimiter_text) > 1:
                    parser.error("--fallback-doi-delimiter must be a single character")
            encoding_cli = getattr(args, "fallback_doi_encoding", None)
            if encoding_cli is not None and not str(encoding_cli).strip():
                parser.error("--fallback-doi-encoding must not be empty")
            for attr_name in ("fallback_doi_col_pmid", "fallback_doi_col_doi"):
                attr_value = getattr(args, attr_name, None)
                if attr_value is None or not str(attr_value).strip():
                    option = attr_name.replace("_", "-")
                    parser.error(f"--{option} must not be empty")
    mapping = {
        "column": f"document.{mode}.column",
        "limit": f"document.{mode}.limit",
    }
    if mode == "pubmed":
        mapping.update(
            {
                "sleep": "document.pubmed.sleep",
                "workers": "document.pubmed.workers",
                "batch_size": "document.pubmed.batch_size",
                "timeout": "sources.pubmed.timeout_read",
            }
        )
    elif mode == "chembl":
        mapping.update(
            {
                "chunk_size": "document.chembl.chunk_size",
                "timeout": "document.chembl.timeout",
            }
        )
    elif mode == "all":
        mapping.update(
            {
                "chunk_size": "document.all.chunk_size",
                "timeout": "document.all.timeout",
                "chembl_chunk_size": "document.all.chunk_size",
                "pubmed_sleep": "document.all.sleep",
                "pubmed_workers": "document.all.workers",
                "pubmed_batch_size": "document.all.batch_size",
                "chembl_timeout": "document.all.timeout",
                "pubmed_timeout": "sources.pubmed.timeout_read",
            }
        )
    mapping |= {
        "openalex_rps": "openalex.rps",
        "crossref_rps": "crossref.rps",
    }
    with setup_cli_logging(
        Path(__file__).with_suffix("").name, log_cfg, getattr(args, "date", None)
    ) as logging_ctx:
        exit_code = run_cli_command(
            args=args,
            parser=parser,
            log_cfg=logging_ctx.log_cfg,
            mapping=mapping,
            run=run,
            logger=logger,
        )
    configure_logger(log_cfg)
    return exit_code


if __name__ == "__main__":  # pragma: no cover - CLI entry point
    raise SystemExit(main())
