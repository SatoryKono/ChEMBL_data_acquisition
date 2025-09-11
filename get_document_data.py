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

Example
-------
Fetch PubMed metadata for identifiers listed in ``pmids.csv``::

    python get_document_data.py pubmed --config config.yaml --input pmids.csv --output output.csv

The input file must contain a ``PMID`` column.

"""

from __future__ import annotations

import argparse
import sys
from typing import Sequence

import pandas as pd

from library.config import (
    Config,
    OpenAlexCfg,
    CrossRefCfg,
    RetryCfg,
    ensure_dirs,
    print_config,
    _serialize_paths,
)
from library.chembl_client import ChemblClient

import requests
from concurrent.futures import ThreadPoolExecutor, as_completed

from library import chembl_library as cl
from library import pubmed_library as pl
from library import semantic_scholar_library as ssl
from library import openalex_crossref_library as ocl
from library import io
from library import document_postprocessing as dp
from library.metadata import Stats, file_sha256, write_meta_yaml
from library.sidecar import SidecarErrors
from library.table_quality import analyze_table_quality
from library.cli import (
    apply_config_overrides,
    build_root_parser,
    configure_logger,
    LoggerConfig,
)
from library.log import logger
from pandera.errors import SchemaErrors
from schemas import DocumentsSchema, normalize_documents

from library import write_csv_deterministic


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
        output = args.output_csv or io.default_output_path(args.input_csv, cfg.io)
        df = normalize_documents(df)
        rows_total = len(df)
        exit_code = 0
        try:
            df = DocumentsSchema.validate(df, lazy=True)
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
        rows_kept = len(df)
        rows_dropped = rows_total - rows_kept
        key_cols = [c for c in ["document_chembl_id"] if c in df.columns]
        csv_path = write_csv_deterministic(
            df,
            output,
            col_order=[
                "document_chembl_id",
                "doi",
                "title",
                "year",
                "month",
                "day",
                "citation",
            ],
            key_cols=key_cols or None,
        )
        logger.info("Wrote %d rows to %s", rows_kept, csv_path)

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
            schema="DocumentsSchema",
        )
    except (FileNotFoundError, ValueError, OSError) as exc:
        logger.error("%s", exc)
        return 1
    try:
        analyze_table_quality(df, table_name=str(output.with_suffix("")))
    except ValueError as exc:
        logger.error("failed to generate quality report: %s", exc)
        return 1
    return exit_code


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
    # Configure session for ChEMBL requests
    client = ChemblClient(cfg.api, RetryCfg())

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
        df = cl.get_documents(  # type: ignore[attr-defined]
            ids,
            client=client,
            cfg=cfg.api,
            chunk_size=args.chunk_size,
            timeout=args.timeout,
        )
    except (requests.RequestException, ValueError) as exc:
        logger.error("failed to retrieve documents: %s", exc)
        return 1
    output = args.output_csv or io.default_output_path(args.input_csv, cfg.io)
    df = normalize_documents(df)
    rows_total = len(df)
    exit_code = 0
    required_cols = set(DocumentsSchema.columns.keys())
    if required_cols.issubset(df.columns):
        try:
            df = DocumentsSchema.validate(df, lazy=True)
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
        key_cols = [c for c in ["document_chembl_id"] if c in df.columns]
        csv_path = write_csv_deterministic(
            df,
            output,
            col_order=[
                "document_chembl_id",
                "doi",
                "title",
                "year",
                "month",
                "day",
                "citation",
            ],
            key_cols=key_cols or None,
        )
        logger.info("Wrote %d rows to %s", rows_kept, csv_path)
    except OSError as exc:
        logger.error("failed to write output CSV: %s", exc)
        return 1

    stats_all: Stats = {
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
        stats=stats_all,
        schema="DocumentsSchema",
    )
    try:
        analyze_table_quality(df, table_name=str(output.with_suffix("")))
    except ValueError as exc:
        logger.error("failed to generate quality report: %s", exc)
        return 1
    return exit_code


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
    # Prepare shared session before performing any API calls
    client = ChemblClient(cfg.api, RetryCfg())

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
        doc_df = cl.get_documents(  # type: ignore[attr-defined]
            ids,
            client=client,
            cfg=cfg.api,
            chunk_size=args.chunk_size,
            timeout=args.timeout,
        )
    except (requests.RequestException, ValueError) as exc:
        logger.error("failed to retrieve documents: %s", exc)
        return 1
    output = args.output_csv or io.default_output_path(args.input_csv, cfg.io)
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
        processed = normalize_documents(processed)
        rows_total = len(processed)
        exit_code = 0
        try:
            processed = DocumentsSchema.validate(processed, lazy=True)
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
            processed = getattr(exc, "validated_data", processed)
            exit_code = 1
        rows_kept = len(processed)
        rows_dropped = rows_total - rows_kept
        try:
            key_cols = [c for c in ["document_chembl_id"] if c in processed.columns]
            csv_path = write_csv_deterministic(
                processed,
                output,
                col_order=[
                    "document_chembl_id",
                    "doi",
                    "title",
                    "year",
                    "month",
                    "day",
                    "citation",
                ],
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
            schema="DocumentsSchema",
        )
        try:
            analyze_table_quality(processed, table_name=str(output.with_suffix("")))
        except ValueError as exc:
            logger.error("failed to generate quality report: %s", exc)
            return 1
        return exit_code

    # Normalise PubMed identifiers to strings to avoid dtype mismatches
    pubmed_ids = pd.to_numeric(doc_df["pubmed_id"], errors="coerce").astype("Int64")
    pmids = pubmed_ids.dropna().astype(str).tolist()
    pub_df = fetch_pubmed_records(
        pmids,
        args.sleep,
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
    processed = normalize_documents(processed)
    rows_total = len(processed)
    exit_code = 0
    try:
        processed = DocumentsSchema.validate(processed, lazy=True)
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
        processed = getattr(exc, "validated_data", processed)
        exit_code = 1
    rows_kept = len(processed)
    rows_dropped = rows_total - rows_kept
    try:
        key_cols = [c for c in ["document_chembl_id"] if c in processed.columns]
        csv_path = write_csv_deterministic(
            processed,
            output,
            col_order=[
                "document_chembl_id",
                "doi",
                "title",
                "year",
                "month",
                "day",
                "citation",
            ],
            key_cols=key_cols or None,
        )
        logger.info("Wrote %d rows to %s", rows_kept, csv_path)
    except OSError as exc:
        logger.error("failed to write output CSV: %s", exc)
        return 1

    stats_all: Stats = {
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
        stats=stats_all,
        schema="DocumentsSchema",
    )
    try:
        analyze_table_quality(processed, table_name=str(output.with_suffix("")))
    except ValueError as exc:
        logger.error("failed to generate quality report: %s", exc)
        return 1
    return exit_code


def build_parser() -> tuple[argparse.ArgumentParser, LoggerConfig]:
    """Create the argument parser for document utilities.

    Returns
    -------
    tuple[argparse.ArgumentParser, LoggerConfig]
        Parser populated with all sub-commands and logging configuration.

    """
    root, log_cfg = build_root_parser()
    parser = argparse.ArgumentParser(
        description="Document data utilities", parents=[root]
    )
    sub = parser.add_subparsers(dest="command", required=True)

    pubmed = sub.add_parser(
        "pubmed", parents=[root], help="Fetch data from PubMed and related APIs"
    )
    pubmed.add_argument(
        "--column", default="PMID", help="Column name containing identifiers"
    )
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

    chembl = sub.add_parser(
        "chembl", parents=[root], help="Fetch document information from ChEMBL"
    )
    chembl.add_argument(
        "--column", default="chembl_id", help="Column name containing identifiers"
    )
    chembl.add_argument(
        "--chunk-size", type=int, default=5, help="Maximum number of IDs per request"
    )
    chembl.add_argument(
        "--timeout",
        type=float,
        default=30.0,
        help="Timeout in seconds for each HTTP request",
    )
    chembl.set_defaults(func=run_chembl)

    all_cmd = sub.add_parser(
        "all", parents=[root], help="Run both ChEMBL and PubMed pipelines"
    )
    all_cmd.add_argument(
        "--column", default="chembl_id", help="Column in the input CSV"
    )
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
    all_cmd.add_argument(
        "--timeout",
        type=float,
        default=30.0,
        help="Timeout in seconds for each HTTP request",
    )
    all_cmd.set_defaults(func=run_all)

    setattr(
        parser,
        "subparsers_map",
        {"pubmed": pubmed, "chembl": chembl, "all": all_cmd},
    )

    return parser, log_cfg


def main(argv: Sequence[str] | None = None) -> int:
    """Command line entry point using :class:`Config` for defaults."""
    parser, log_cfg = build_parser()
    args = parser.parse_args(argv)
    log_cfg.level = args.log_level
    logger = configure_logger(log_cfg)
    logger.info("pipeline start run_id=%s", log_cfg.run_id, extra={"event": "start"})
    subparser_map = getattr(parser, "subparsers_map", {})
    subparser = subparser_map.get(args.command, parser)
    try:
        cfg: Config = apply_config_overrides(
            args, subparser, args.config, mapping={"timeout": "api.timeout_read"}
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
