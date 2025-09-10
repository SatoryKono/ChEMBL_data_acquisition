"""Command line interface for retrieving document metadata from external sources.

The tool integrates :mod:`library.pubmed_library` and
:mod:`library.chembl_library` to collect information about publications from
several public APIs.  The interface mirrors :mod:`get_target_data.py` and
provides three sub-commands:

``pubmed``
    Query PubMed, Semantic Scholar, OpenAlex and CrossRef for a list of PMIDs.
``chembl``
    Retrieve document information from the ChEMBL API.
``all``
    Run the ChEMBL and PubMed pipelines and merge the results.

Example:
-------
Fetch PubMed metadata for identifiers listed in ``pmids.csv``::

    python get_document_data.py pubmed pmids.csv output.csv

The input file must contain a ``PMID`` column.

"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path
from typing import Sequence

import pandas as pd

from library.config import Config, ensure_dirs, OpenAlexCfg, CrossRefCfg

import requests
from concurrent.futures import ThreadPoolExecutor, as_completed

from library import chembl_library as cl
from library import pubmed_library as pl
from library import semantic_scholar_library as ssl
from library import openalex_crossref_library as ocl
from library import io
from library import document_postprocessing as dp
from library.table_quality import analyze_table_quality
from library.cli import (
    apply_config_overrides,
    build_parser as base_parser,
    configure_logging,
)

logger = logging.getLogger(__name__)


def fetch_pubmed_records(
    pmids: list[str],
    sleep: float,


    openalex_cfg: OpenAlexCfg,
    crossref_cfg: CrossRefCfg,
    max_workers: int = 1,
    batch_size: int = 100,
) -> pd.DataFrame:
    """Retrieve metadata for a list of PubMed identifiers.

    Parameters
    ----------
    pmids:
        Identifiers to query.
    sleep:

        Seconds to pause between PubMed and Semantic Scholar requests.
    openalex_cfg:
        Configuration for OpenAlex API access.
    crossref_cfg:
        Configuration for CrossRef API access.

    max_workers:
        Maximum number of concurrent threads.
    batch_size:
        Maximum number of PMIDs per PubMed request.

    Returns
    -------
    pandas.DataFrame
        Combined metadata from the different sources.

    """

    def _fetch_batch(batch: list[str]) -> list[dict[str, str]]:
        """Fetch metadata for a batch of PMIDs.

        Each worker opens its own :class:`requests.Session` and retrieves PubMed
        entries for all PMIDs in ``batch`` using a single request. Metadata from
        Semantic Scholar, OpenAlex and CrossRef are then fetched individually
        for each PMID. Exceptions are logged so a failure in one batch does not
        abort the whole process.
        """
        try:
            with requests.Session() as session:
                pubmed_list = pl.fetch_pubmed_batch(session, batch, sleep)
                pmids_in_batch = [p.get("PubMed.PMID", "") for p in pubmed_list]

                # Fetch Semantic Scholar data in a single batch
                semsch_list = ssl.fetch_semantic_scholar_batch(
                    session, pmids_in_batch, sleep
                )

                # Create a map for easy lookup
                semsch_map = {s.get("scholar.PMID"): s for s in semsch_list}

                combined_records: list[dict[str, str]] = []
                for pubmed in pubmed_list:
                    pmid = pubmed.get("PubMed.PMID", "")
                    semsch = semsch_map.get(pmid, {})

                    # Still fetching these individually for now

                    openalex = ocl.fetch_openalex(session, pmid, openalex_cfg)
                    doi = pubmed.get("PubMed.DOI") or semsch.get("scholar.DOI") or ""
                    crossref = ocl.fetch_crossref(session, doi, crossref_cfg)


                    combined: dict[str, str] = {}
                    combined.update(pubmed)
                    combined.update(semsch)
                    combined.update(openalex)
                    combined.update(crossref)
                    combined_records.append(combined)
                return combined_records
        except requests.RequestException as exc:  # pragma: no cover - network errors
            logger.warning("failed to fetch PMIDs %s: %s", batch, exc)
            return [{} for _ in batch]

    if not pmids:
        return pd.DataFrame()

    records: list[dict[str, str]] = []
    batches = [pmids[i : i + batch_size] for i in range(0, len(pmids), batch_size)]
    total = len(pmids)
    processed = 0
    with ThreadPoolExecutor(max_workers=max_workers) as ex:
        futures = {ex.submit(_fetch_batch, batch): len(batch) for batch in batches}
        for future in as_completed(futures):
            batch_len = futures[future]
            records.extend(future.result())
            processed += batch_len
            percent = processed / total * 100
            logger.info("Processed %d/%d documents (%.1f%%)", processed, total, percent)
    if not records:
        return pd.DataFrame()
    return pd.DataFrame(records)


def run_pubmed(cfg: Config, args: argparse.Namespace) -> int:
    """Execute the ``pubmed`` sub-command.

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
    try:
        pmids = io.read_ids(
            args.input_csv,
            column=args.column,
            cfg=cfg.io,
            sep=args.sep,
            encoding=args.encoding,
        )
        df = fetch_pubmed_records(
            pmids,
            args.sleep,

            cfg.openalex,
            cfg.crossref,
            args.workers,
            args.batch_size,

        )
        output = args.output_csv or io.default_output_path(args.input_csv)
        io.write_csv(df, output, cfg=cfg, sep=args.sep, encoding=args.encoding)
        logger.info("Wrote %d rows to %s", len(df), output)
    except (FileNotFoundError, ValueError, OSError) as exc:
        logger.error("%s", exc)
        return 1
    try:
        analyze_table_quality(df, table_name=str(output.with_suffix("")))
    except ValueError as exc:
        logger.error("failed to generate quality report: %s", exc)
        return 1
    return 0


def run_chembl(cfg: Config, args: argparse.Namespace) -> int:
    """Execute the ``chembl`` sub-command.

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
    try:
        ids = io.read_ids(
            args.input_csv,
            column=args.column,
            cfg=cfg.io,
            sep=args.sep,
            encoding=args.encoding,
        )
    except (FileNotFoundError, ValueError) as exc:
        logger.error("%s", exc)
        return 1

    try:
        df = cl.get_documents(ids, chunk_size=args.chunk_size)  # type: ignore[attr-defined]
    except (requests.RequestException, ValueError) as exc:
        logger.error("failed to retrieve documents: %s", exc)
        return 1
    output = args.output_csv or io.default_output_path(args.input_csv)
    try:
        io.write_csv(df, output, cfg=cfg, sep=args.sep, encoding=args.encoding)
        logger.info("Wrote %d rows to %s", len(df), output)
    except OSError as exc:
        logger.error("failed to write output CSV: %s", exc)
        return 1
    try:
        analyze_table_quality(df, table_name=str(output.with_suffix("")))
    except ValueError as exc:
        logger.error("failed to generate quality report: %s", exc)
        return 1
    return 0


def run_all(cfg: Config, args: argparse.Namespace) -> int:
    """Run ChEMBL and PubMed pipelines and merge their outputs.

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
    try:
        ids = io.read_ids(
            args.input_csv,
            column=args.column,
            cfg=cfg.io,
            sep=args.sep,
            encoding=args.encoding,
        )
    except (FileNotFoundError, ValueError) as exc:
        logger.error("%s", exc)
        return 1

    try:
        doc_df = cl.get_documents(ids, chunk_size=args.chunk_size)  # type: ignore[attr-defined]
    except (requests.RequestException, ValueError) as exc:
        logger.error("failed to retrieve documents: %s", exc)
        return 1
    output = args.output_csv or io.default_output_path(args.input_csv)
    if doc_df.empty or "pubmed_id" not in doc_df:
        processed = dp.postprocess_documents(doc_df)
        # Merge any columns not covered by ``postprocess_documents`` back
        # into the result so the output retains all original fields.
        extra_cols = [c for c in doc_df.columns if c not in processed.columns]
        if extra_cols:
            processed = processed.merge(
                doc_df[["document_chembl_id"] + extra_cols],
                on="document_chembl_id",
                how="left",
            )
        try:
            io.write_csv(
                processed, output, cfg=cfg, sep=args.sep, encoding=args.encoding
            )
            logger.info("Wrote %d rows to %s", len(processed), output)
        except OSError as exc:
            logger.error("failed to write output CSV: %s", exc)
            return 1
        try:
            analyze_table_quality(processed, table_name=str(output.with_suffix("")))
        except ValueError as exc:
            logger.error("failed to generate quality report: %s", exc)
            return 1
        return 0

    # Normalise PubMed identifiers to strings to avoid dtype mismatches
    pubmed_ids = pd.to_numeric(doc_df["pubmed_id"], errors="coerce").astype("Int64")
    pmids = pubmed_ids.dropna().astype(str).tolist()
    pub_df = fetch_pubmed_records(
        pmids,

        cfg.openalex,
        cfg.crossref,
        args.workers,
        args.batch_size,

    )
    doc_df["pubmed_id"] = pubmed_ids.astype(str)
    if not pub_df.empty and "PubMed.PMID" in pub_df.columns:
        pub_df["PubMed.PMID"] = (
            pd.to_numeric(pub_df["PubMed.PMID"], errors="coerce")
            .astype("Int64")
            .astype(str)
        )
        merged = doc_df.merge(
            pub_df, how="left", left_on="pubmed_id", right_on="PubMed.PMID"
        )
    else:
        merged = doc_df
    processed = dp.postprocess_documents(merged)
    # Append any additional columns from the merged table that were not
    # included in the post-processing result.
    extra_cols = [c for c in merged.columns if c not in processed.columns]
    if extra_cols:
        processed = processed.merge(
            merged[["document_chembl_id"] + extra_cols],
            on="document_chembl_id",
            how="left",
        )
    try:
        io.write_csv(processed, output, cfg=cfg, sep=args.sep, encoding=args.encoding)
        logger.info("Wrote %d rows to %s", len(processed), output)
    except OSError as exc:
        logger.error("failed to write output CSV: %s", exc)
        return 1
    try:
        analyze_table_quality(processed, table_name=str(output.with_suffix("")))
    except ValueError as exc:
        logger.error("failed to generate quality report: %s", exc)
        return 1
    return 0


def build_parser() -> argparse.ArgumentParser:
    """Create the argument parser for document utilities.

    Returns
    -------
    argparse.ArgumentParser
        Parser populated with all sub-commands.

    """
    parser = base_parser("Document data utilities", column="chembl_id", chunk_size=5)
    sub = parser.add_subparsers(dest="command", required=True)

    pubmed = sub.add_parser("pubmed", help="Fetch data from PubMed and related APIs")
    pubmed.add_argument(
        "--input",
        dest="input_csv",
        type=Path,
        default=Path("input.csv"),
        help="CSV with a PMID column",
    )
    pubmed.add_argument(
        "--output",
        dest="output_csv",
        type=Path,
        default=None,
        help="Destination CSV file (default: auto-generate)",
    )
    pubmed.add_argument(
        "--column", default="PMID", help="Column name containing identifiers"
    )
    pubmed.add_argument("--sep", default=",", help="CSV delimiter")
    pubmed.add_argument("--encoding", default="utf8", help="File encoding")
    pubmed.add_argument(
        "--sleep", type=float, default=5.0, help="Seconds to sleep between requests"
    )
    pubmed.add_argument(
        "--workers", type=int, default=1, help="Number of concurrent requests"
    )
    pubmed.add_argument(
        "--batch-size",
        type=int,
        default=100,
        help="Maximum PMIDs per PubMed request",
    )
    pubmed.set_defaults(func=run_pubmed)

    chembl = sub.add_parser("chembl", help="Fetch document information from ChEMBL")
    chembl.add_argument(
        "--input",
        dest="input_csv",
        type=Path,
        default=Path("input.csv"),
        help="CSV with document_chembl_id column",
    )
    chembl.add_argument(
        "--output",
        dest="output_csv",
        type=Path,
        default=None,
        help="Destination CSV file (default: auto-generate)",
    )
    chembl.add_argument(
        "--column", default="chembl_id", help="Column name containing identifiers"
    )
    chembl.add_argument("--sep", default=",", help="CSV delimiter")
    chembl.add_argument("--encoding", default="utf8", help="File encoding")
    chembl.add_argument(
        "--chunk-size", type=int, default=5, help="Maximum number of IDs per request"
    )
    chembl.set_defaults(func=run_chembl)

    all_cmd = sub.add_parser("all", help="Run both ChEMBL and PubMed pipelines")
    all_cmd.add_argument(
        "--input",
        dest="input_csv",
        type=Path,
        default=Path("input.csv"),
        help="CSV with document_chembl_id column",
    )
    all_cmd.add_argument(
        "--output",
        dest="output_csv",
        type=Path,
        default=None,
        help="Destination CSV file (default: auto-generate)",
    )
    all_cmd.add_argument(
        "--column", default="chembl_id", help="Column in the input CSV"
    )
    all_cmd.add_argument("--sep", default=",", help="CSV delimiter")
    all_cmd.add_argument("--encoding", default="utf8", help="File encoding")
    all_cmd.add_argument(
        "--chunk-size", type=int, default=5, help="Maximum IDs per request"
    )
    all_cmd.add_argument(
        "--sleep",
        type=float,
        default=5.0,
        help="Seconds to sleep between PubMed requests",
    )
    all_cmd.add_argument(
        "--workers", type=int, default=1, help="Number of concurrent PubMed requests"
    )
    all_cmd.add_argument(
        "--batch-size",
        type=int,
        default=50,
        help="Maximum PMIDs per PubMed request",
    )
    all_cmd.set_defaults(func=run_all)

    return parser


def main(argv: Sequence[str] | None = None) -> int:
    """Command line entry point using :class:`Config` for defaults."""
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        cfg: Config = apply_config_overrides(args, parser, args.config)
        ensure_dirs(cfg)
        configure_logging(args.log_level, fmt=cfg.log.format, datefmt=cfg.log.datefmt)
    except (ValueError, TypeError) as exc:
        logger.error("%s", exc)
        return 1
    return args.func(cfg, args)


if __name__ == "__main__":  # pragma: no cover - CLI entry point
    raise SystemExit(main())
