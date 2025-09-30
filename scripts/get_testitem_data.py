"""Command line interface for retrieving ChEMBL and PubChem compound data."""

from __future__ import annotations

from pathlib import Path
import sys

import argparse
from collections import OrderedDict, deque
from collections.abc import Iterable, Iterator, Sequence
from dataclasses import dataclass
from itertools import islice, tee
from typing import Any

import pandas as pd
import requests

_REPO_ROOT = Path(__file__).resolve().parents[1]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from library import chembl_library as cl
from library import cli
from library import io

from library.clients import pubchem as pc
from library.chembl_client import ChemblClient

from library.cli import (
    LoggerConfig,
    configure_logger,
)
from library.cli import (
    build_parser as base_parser,
)
from library.config import ApiCfg, Config, IoCfg, _serialize_paths, ensure_dirs, print_config
from library.log import logger
from schemas import TestitemsSchema
from library.testitem_pipeline import (
    PARENT_LOOKUP_SOURCE_CACHE,
    PARENT_LOOKUP_SOURCE_LOOKUP,
    PARENT_LOOKUP_SOURCE_PARTIAL,
    PARENT_LOOKUP_SOURCE_SKIPPED,
    PARENT_LOOKUP_SOURCE_SYNC,
    ParentEnrichmentPreparation,
    ParentEnrichmentResult,
    ParentLookupPreparedData,
    ParentLookupStats,
    augment_pubchem,
    apply_testitem_enrichment,
    finalize_output,
    prepare_parent_enrichment,
    run_parent_enrichment,
)


# ===== Parameters =====

_FETCH_ERROR_SAMPLE_SIZE = 10

@dataclass
class ReadInputIdsResult:
    """Container holding the identifier iterator and a diagnostic sample."""

    ids_iter: Iterator[str]
    sample_ids: tuple[str, ...]
    requested_ids: tuple[str, ...]


def read_input_ids(
    input_csv: Path,
    *,
    column: str,
    io_cfg: IoCfg,
    limit: int | None,
    offset: int = 0,
) -> tuple[int, ReadInputIdsResult | None]:
    """Load identifiers from ``input_csv`` honouring ``limit`` when provided."""

    try:
        ids_iter = io.read_ids(
            input_csv,
            column=column,
            cfg=io_cfg,
            keep_na_markers=io_cfg.keep_na_markers,
        )
        if offset:
            ids_iter = islice(ids_iter, offset, None)
            logger.info("process_offset", offset=offset)
        if limit is not None:
            ids_iter = islice(ids_iter, limit)
        ids_iter, sample_iter, capture_iter = tee(ids_iter, 3)
        sample_ids = tuple(islice(sample_iter, _FETCH_ERROR_SAMPLE_SIZE))
        requested_ids = tuple(capture_iter)
    except (FileNotFoundError, ValueError) as exc:
        logger.error(
            "read_fail",
            error=str(exc),
            path=str(input_csv),
        )
        return 1, None

    return 0, ReadInputIdsResult(
        ids_iter=ids_iter,
        sample_ids=sample_ids,
        requested_ids=requested_ids,
    )


def _integrate_missing_identifiers(
    df: pd.DataFrame,
    *,
    missing_ids: Sequence[str],
    requested_ids: Sequence[str],
) -> pd.DataFrame:
    """Attach placeholder rows for ``missing_ids`` and restore input ordering."""

    if missing_ids:
        columns = list(df.columns)
        if "molecule_chembl_id" not in columns:
            columns.append("molecule_chembl_id")
        missing_frame = pd.DataFrame(
            pd.NA,
            index=range(len(missing_ids)),
            columns=columns,
        )
        missing_frame["molecule_chembl_id"] = list(missing_ids)
        df = pd.concat([df, missing_frame], ignore_index=True, sort=False)
    else:
        df = df.reset_index(drop=True)

    if not requested_ids or "molecule_chembl_id" not in df.columns:
        return df

    order_map: dict[str, deque[int]] = {}
    for index, identifier in enumerate(requested_ids):
        if identifier is None:
            continue
        key = str(identifier).strip()
        if not key:
            continue
        upper_key = key.upper()
        order_map.setdefault(upper_key, deque()).append(index)

    if not order_map:
        return df

    fallback_index = len(requested_ids)

    def _pop_order(value: Any) -> int:
        if value is None or pd.isna(value):
            return fallback_index
        key = str(value).strip()
        if not key:
            return fallback_index
        upper_key = key.upper()
        queue = order_map.get(upper_key)
        if queue:
            return queue.popleft()
        return fallback_index

    order_series = df["molecule_chembl_id"].map(_pop_order)
    ordered = df.assign(_input_order=order_series).sort_values(
        "_input_order", kind="stable"
    )
    ordered = ordered.drop(columns="_input_order").reset_index(drop=True)
    return ordered


def fetch_testitems(
    ids_iter: Iterable[str],
    *,
    api_cfg: ApiCfg,
    batch_size: int,
    timeout: float,
    client: ChemblClient,
    sample_ids: Sequence[str],
    fields: Sequence[str] | None,
    page_limit: int,
) -> tuple[int, pd.DataFrame | None]:
    """Retrieve ChEMBL test item records for ``ids_iter``."""

    logger.info("chembl_fetch_start", batch_size=batch_size)
    ids_iter, capture_iter = tee(ids_iter)
    requested_ids_raw = list(capture_iter)
    try:
        df = cl.get_testitem(
            ids_iter,
            cfg=api_cfg,
            client=client,
            chunk_size=batch_size,
            timeout=timeout,
            fields=fields,
            page_limit=page_limit,
        )
    except (requests.RequestException, ValueError) as exc:
        logger.error(
            "testitem_fetch_failed",
            error=str(exc),
            batch_size=batch_size,
            timeout=timeout,
            sample_ids=list(sample_ids),
        )
        return 1, None

    rows = len(df)
    logger.info("chembl_fetch_done", rows=rows)
    logger.info("identifiers_retrieved", count=rows)

    original_id_lookup: dict[str, str] = {}
    requested_ids: list[str] = []
    for identifier in requested_ids_raw:
        normalised = str(identifier).strip().upper()
        requested_ids.append(normalised)
        original_id_lookup[normalised] = str(identifier)

    if "molecule_chembl_id" not in df.columns:
        df["molecule_chembl_id"] = pd.Series(dtype="string")
    df["molecule_chembl_id"] = df["molecule_chembl_id"].astype("string").str.upper()

    retrieved_ids = set(df["molecule_chembl_id"].dropna())
    missing_ids = [identifier for identifier in requested_ids if identifier not in retrieved_ids]
    if missing_ids:
        missing_ids_original = [
            original_id_lookup.get(identifier, identifier) for identifier in missing_ids
        ]
        logger.warning(
            "chembl_missing_identifiers",
            missing_count=len(missing_ids),
            missing_ids=missing_ids,
            missing_ids_original=missing_ids_original,
        )
    else:
        logger.info("chembl_all_identifiers_retrieved")

    duplicates_mask = df["molecule_chembl_id"].duplicated(keep=False)
    if duplicates_mask.any():
        duplicate_ids = (
            df.loc[duplicates_mask, "molecule_chembl_id"].dropna().unique().tolist()
        )
        logger.warning(
            "chembl_duplicate_identifiers",
            duplicate_count=len(duplicate_ids),
            duplicate_ids=duplicate_ids,
        )
        df = df.drop_duplicates(subset="molecule_chembl_id", keep="first")

    index = pd.Index(requested_ids, name="molecule_chembl_id")
    df = df.set_index("molecule_chembl_id", drop=False).reindex(index)
    df["molecule_chembl_id"] = pd.Series(
        requested_ids,
        index=index,
        dtype="string",
    )
    df = df.where(pd.notna(df), pd.NA)
    df = df.convert_dtypes()
    df["molecule_chembl_id"] = df["molecule_chembl_id"].astype("string")
    df = df.reset_index(drop=True)

    return 0, df


def run_chembl(cfg: Config, args: argparse.Namespace) -> int:
    """Execute compound retrieval from the ChEMBL API and augment with PubChem data."""

    limit = cfg.testitem.limit
    if limit is not None and limit < 0:
        logger.error(
            "invalid_limit",
            section="testitem.limit",
            limit=limit,
        )
        return 1

    api_overrides: dict[str, Any] = {}
    if cfg.testitem.retries is not None:
        api_overrides["retries"] = cfg.testitem.retries
    if cfg.testitem.backoff_factor is not None:
        api_overrides["backoff_factor"] = cfg.testitem.backoff_factor
    api_cfg = cfg.api.model_copy(update=api_overrides) if api_overrides else cfg.api

    pc.init_session(api_cfg, cfg.retry)

    requested_ids: tuple[str, ...] = ()
    missing_ids: list[str] = []

    with ChemblClient(api_cfg, cfg.retry, cfg.chembl) as client:
        read_status, read_result = read_input_ids(
            args.input_csv,
            column=cfg.testitem.column,
            io_cfg=cfg.io,
            limit=limit,
            offset=getattr(args, "offset", cfg.testitem.offset),
        )
        if read_status != 0 or read_result is None:
            return read_status

        requested_ids = read_result.requested_ids
        fetch_status, df = fetch_testitems(
            read_result.ids_iter,
            api_cfg=api_cfg,
            batch_size=cfg.testitem.batch_size,
            timeout=cfg.testitem.timeout,
            client=client,
            sample_ids=read_result.sample_ids,
            fields=cfg.testitem.fields,
            page_limit=cfg.testitem.request_limit,
        )
        if fetch_status != 0 or df is None:
            return fetch_status

        requested_unique = OrderedDict()
        for identifier in requested_ids:
            key = identifier.strip()
            if not key:
                continue
            upper_key = key.upper()
            requested_unique.setdefault(upper_key, identifier)

        fetched_ids: set[str] = set()
        if "molecule_chembl_id" in df.columns:
            fetched_series = df["molecule_chembl_id"].dropna().astype(str).str.strip()
            fetched_ids = {value.upper() for value in fetched_series if value}

        missing_keys = [key for key in requested_unique if key not in fetched_ids]
        missing_ids = [requested_unique[key] for key in missing_keys]
        if missing_ids:
            logger.warning(
                "chembl_missing_identifiers",
                count=len(missing_ids),
                identifiers=missing_ids,
            )

        rows = len(df)
        if limit is not None:
            logger.info("process_limit", limit=min(limit, rows))

        enrichment_sources = getattr(cfg.testitem_molecule_enrichment, "sources", None)
        hierarchy_lookup_path = (
            getattr(enrichment_sources, "molecule_hierarchy_path", None)
            if enrichment_sources is not None
            else None
        )

        prep_status, prep = prepare_parent_enrichment(
            df,
            catalog_cfg=cfg.molecule_catalog,
            io_cfg=cfg.io,
            api_cfg=api_cfg,
            timeout=cfg.testitem.timeout,
            client=client,
            hierarchy_lookup_path=hierarchy_lookup_path,
        )
        if prep_status != 0 or prep is None:
            return prep_status

        parent_stats = prep.parent_stats
        parent_status, parent_result = run_parent_enrichment(
            prep,
            client=client,
            api_cfg=api_cfg,
            catalog_cfg=cfg.molecule_catalog,
            timeout=cfg.testitem.timeout,
        )
        if parent_status != 0 or parent_result is None:
            return parent_status

        df = parent_result.df
        parent_stats = parent_result.parent_stats

        df = augment_pubchem(
            df,
            pubchem_cfg=cfg.pubchem,
            api_cfg=api_cfg,
            timeout=cfg.testitem.timeout,
            client=client,
            fields=cfg.testitem.fields,
            request_limit=cfg.testitem.request_limit,
        )

        enrichment_status, enriched_df = apply_testitem_enrichment(
            df,
            enrichment_cfg=cfg.testitem_molecule_enrichment,
            io_cfg=cfg.io,
        )
        if enrichment_status != 0 or enriched_df is None:
            return enrichment_status

        df = enriched_df

    df = _integrate_missing_identifiers(
        df,
        missing_ids=missing_ids,
        requested_ids=requested_ids,
    )

    output_path = (
        Path(args.output_csv)
        if args.output_csv
        else io.default_output_path(args.input_csv, cfg.io)
    )
    rows_total = len(df)
    exit_code = finalize_output(
        df,
        cfg=cfg,
        output=output_path,
        parent_stats=parent_stats,
        input_csv=args.input_csv,
        rows_total=rows_total,
        missing_ids=missing_ids,
    )
    return exit_code


def build_parser() -> tuple[argparse.ArgumentParser, LoggerConfig]:
    """Create the command-line argument parser."""
    parser, log_cfg = base_parser(
        "ChEMBL and PubChem compound data utilities",
        column="molecule_chembl_id",
        chunk_size=1000,
        size_option="--batch-size",
        size_dest="batch_size",
    )
    parser.add_argument(
        "--timeout",
        type=float,
        default=30.0,
        help="Timeout in seconds for each HTTP request",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Maximum number of identifiers to process",
    )
    parser.add_argument(
        "--offset",
        type=int,
        default=0,
        help="Number of identifiers to skip before processing",
    )
    parser.set_defaults(func=run_chembl)
    return parser, log_cfg


def main(argv: Sequence[str] | None = None) -> int:
    """Command line entry point using :class:`Config` for defaults."""
    parser, log_cfg = build_parser()
    args = parser.parse_args(argv)
    if args.limit is not None and args.limit <= 0:
        parser.error("--limit must be a positive integer")
    if args.offset < 0:
        parser.error("--offset must be zero or a positive integer")
    log_cfg.level = args.log_level
    logger = configure_logger(log_cfg)
    logger.info("pipeline_start", run_id=log_cfg.run_id)
    try:
        cfg: Config = cli.apply_config_overrides(
            args,
            parser,
            args.config,
            mapping={
                "timeout": "testitem.timeout",
                "column": "testitem.column",
                "batch_size": "testitem.batch_size",
                "limit": "testitem.limit",
                "offset": "testitem.offset",
            },
        )
        if args.print_config:
            print_config(cfg)
            configure_logger(log_cfg, fmt=cfg.log.format, datefmt=cfg.log.datefmt)
            logger.info("pipeline_done", run_id=log_cfg.run_id)
            return 0
        ensure_dirs(cfg)
        logger = configure_logger(log_cfg, fmt=cfg.log.format, datefmt=cfg.log.datefmt)
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
        )
        logger.info("pipeline_fail", run_id=log_cfg.run_id)
        return 1
    exit_code: int = args.func(cfg, args)
    if exit_code == 0:
        logger.info("pipeline_done", run_id=log_cfg.run_id)
    else:
        logger.info("pipeline_fail", run_id=log_cfg.run_id)
    return exit_code


if __name__ == "__main__":  # pragma: no cover - CLI entry point
    raise SystemExit(main())
