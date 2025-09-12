"""Utilities for retrieving and merging publication metadata.

This module exposes a command line interface and re-exports the
implementation located in :mod:`library.pubmed`.
"""

from __future__ import annotations

import argparse
from collections.abc import Sequence
from datetime import date
from pathlib import Path

import pandas as pd

from .cli import LoggerConfig, configure_logger
from .cli import build_parser as base_parser
from .config import ApiCfg, Config, session_with_retry
from .csv_utils import write_csv_deterministic
from .log import logger
from .pubmed import (
    EMPTY_PUBMED,
    _do_request,
    _handle_response,
    _make_request,
    combine,
    fetch_crossref,
    fetch_openalex,
    fetch_pubmed,
    fetch_pubmed_batch,
    fetch_semantic_scholar,
    fetch_semantic_scholar_batch,
    find_all,
    find_one,
    merge_records,
    parse_pubmed_article,
    print_results,
    read_pmids,
    text_or_none,
)
from .rate_limiter import get_limiter

__all__ = [
    "read_pmids",
    "_make_request",
    "_handle_response",
    "_do_request",
    "fetch_pubmed_batch",
    "fetch_pubmed",
    "fetch_semantic_scholar",
    "fetch_semantic_scholar_batch",
    "fetch_openalex",
    "fetch_crossref",
    "text_or_none",
    "combine",
    "find_one",
    "find_all",
    "parse_pubmed_article",
    "EMPTY_PUBMED",
    "merge_records",
    "print_results",
    "parse_args",
    "main",
]


def parse_args(
    argv: Sequence[str] | None = None,
) -> tuple[argparse.Namespace, LoggerConfig]:
    """Parse command-line arguments."""
    parser, log_cfg = base_parser("Fetch publication metadata by PMID", column="PMID")
    parser.add_argument(
        "--input-csv",
        dest="input_csv",
        help="Input CSV path with PMID column",
    )
    args = parser.parse_args(argv)
    return args, log_cfg


def main(argv: Sequence[str] | None = None) -> int:
    """Command-line interface entry point."""
    args, log_cfg = parse_args(argv)
    log_cfg.level = args.log_level
    configure_logger(log_cfg)

    cfg = Config(api=ApiCfg(user_agent="chembl-da/0.1 (mailto:info@example.org)"))
    pubmed_cfg = cfg.pubmed
    semsch_cfg = cfg.semantic_scholar
    limiter = get_limiter("global", cfg.rate.global_rps, cfg.rate.global_burst)
    delay = 1.0 / cfg.rate.global_rps if cfg.rate.global_rps > 0 else 0.0

    pmid_df = read_pmids(args.input_csv, cfg=pubmed_cfg)
    pmids = pmid_df["PMID"].tolist()
    openalex_limiter = get_limiter("openalex", cfg.openalex.rps, cfg.openalex.burst)
    crossref_limiter = get_limiter("crossref", cfg.crossref.rps, cfg.crossref.burst)

    records: list[dict[str, str]] = []
    batch_size = 100
    with session_with_retry(cfg.api, cfg.retry) as session:
        for i in range(0, len(pmids), batch_size):
            batch_pmids = pmids[i : i + batch_size]
            limiter.acquire()
            pubmed_list = fetch_pubmed_batch(
                session, batch_pmids, delay, cfg=pubmed_cfg
            )
            limiter.acquire()
            semsch_list = fetch_semantic_scholar_batch(
                session, batch_pmids, delay, cfg=semsch_cfg
            )
            semsch_map = {s.get("scholar.PMID"): s for s in semsch_list}
            for pubmed in pubmed_list:
                pmid = pubmed.get("PubMed.PMID", "")
                semsch = semsch_map.get(pmid, {})

                openalex = fetch_openalex(
                    session, pmid, cfg=cfg.openalex, limiter=openalex_limiter
                )
                doi = pubmed.get("PubMed.DOI") or semsch.get("scholar.DOI") or ""
                crossref = fetch_crossref(
                    session, doi, cfg=cfg.crossref, limiter=crossref_limiter
                )

                combined = merge_records(pubmed, semsch, openalex, crossref)
                print_results([combined])
                records.append(combined)

    df = pd.DataFrame.from_records(records)
    output_path = (
        Path(args.output_csv)
        if args.output_csv
        else Path(f"output_{Path(args.input_csv).stem}_{date.today():%Y%m%d}.csv")
    )
    write_csv_deterministic(df, output_path, key_cols=sorted(df.columns))
    logger.info("file_written", path=str(output_path))
    return 0


if __name__ == "__main__":  # pragma: no cover - CLI entry point
    raise SystemExit(main())
