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

    python get_document_data.py pubmed pmids.csv output.csv

The input file must contain a ``PMID`` column.
"""

from __future__ import annotations

import argparse
import csv
import logging
from pathlib import Path
from typing import Sequence

import pandas as pd
import requests
from concurrent.futures import ThreadPoolExecutor, as_completed

from library import chembl_library as cl
from library import pubmed_library as pl
from library import semantic_scholar_library as ssl
from library import openalex_crossref_library as ocl

logger = logging.getLogger(__name__)


def read_ids(
    path: str | Path,
    column: str = "document_chembl_id",
    sep: str = ",",
    encoding: str = "utf8",
) -> list[str]:
    """Read identifier values from a CSV file.

    Parameters
    ----------
    path:
        CSV file path.
    column:
        Name of the column containing identifiers.
    sep:
        Field delimiter used in the CSV file.
    encoding:
        File encoding of the CSV file.

    Returns
    -------
    list[str]
        Identifier values in the order they appear. Empty strings and ``"#N/A"``
        entries are discarded.
    """

    try:
        with Path(path).open("r", encoding=encoding, newline="") as fh:
            reader = csv.DictReader(fh, delimiter=sep)
            if reader.fieldnames is None or column not in reader.fieldnames:
                raise ValueError(f"column '{column}' not found in {path}")
            ids: list[str] = []
            for row in reader:
                value = (row.get(column) or "").strip()
                if value and value != "#N/A":
                    ids.append(value)
            return ids
    except FileNotFoundError as exc:  # pragma: no cover - trivial
        raise FileNotFoundError(f"input file not found: {path}") from exc
    except csv.Error as exc:  # pragma: no cover - malformed CSV
        raise ValueError(f"malformed CSV in file: {path}: {exc}") from exc


def fetch_pubmed_records(
    pmids: list[str],
    sleep: float,
    max_workers: int = 1,
    batch_size: int = 100,
) -> pd.DataFrame:
    """Retrieve metadata for a list of PubMed identifiers.

    Parameters
    ----------
    pmids:
        Identifiers to query.
    sleep:
        Seconds to pause between API requests.
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
                semsch_list = ssl.fetch_semantic_scholar_batch(session, pmids_in_batch, sleep)
                
                # Create a map for easy lookup
                semsch_map = {s.get("scholar.PMID"): s for s in semsch_list}
                
                combined_records: list[dict[str, str]] = []
                for pubmed in pubmed_list:
                    pmid = pubmed.get("PubMed.PMID", "")
                    semsch = semsch_map.get(pmid, {})
                    
                    # Still fetching these individually for now
                    openalex = ocl.fetch_openalex(session, pmid, sleep)
                    doi = pubmed.get("PubMed.DOI") or semsch.get("scholar.DOI") or ""
                    crossref = ocl.fetch_crossref(session, doi, sleep)
                    
                    combined: dict[str, str] = {}
                    combined.update(pubmed)
                    combined.update(semsch)
                    combined.update(openalex)
                    combined.update(crossref)
                    combined_records.append(combined)
                return combined_records
        except Exception as exc:  # pragma: no cover - network errors
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
            logger.info(
                "Processed %d/%d documents (%.1f%%)", processed, total, percent
            )
    if not records:
        return pd.DataFrame()
    return pd.DataFrame(records)


def run_pubmed(args: argparse.Namespace) -> int:
    """Execute the ``pubmed`` sub-command."""

    try:
        pmids = pl.read_pmids(args.input_csv)
        df = fetch_pubmed_records(pmids, args.sleep, args.workers, args.batch_size)
        df.to_csv(args.output_csv, index=False)
        logger.info("Wrote %d rows to %s", len(df), args.output_csv)
        return 0
    except (FileNotFoundError, ValueError, OSError) as exc:
        logger.error("%s", exc)
        return 1


def run_chembl(args: argparse.Namespace) -> int:
    """Execute the ``chembl`` sub-command."""

    try:
        ids = read_ids(
            args.input_csv, column=args.column, sep=args.sep, encoding=args.encoding
        )
    except (FileNotFoundError, ValueError) as exc:
        logger.error("%s", exc)
        return 1

    df = cl.get_documents(ids, chunk_size=args.chunk_size)
    try:
        df.to_csv(args.output_csv, index=False)
        logger.info("Wrote %d rows to %s", len(df), args.output_csv)
        return 0
    except OSError as exc:
        logger.error("failed to write output CSV: %s", exc)
        return 1


def run_all(args: argparse.Namespace) -> int:
    """Run ChEMBL and PubMed pipelines and merge their outputs."""

    try:
        ids = read_ids(
            args.input_csv, column=args.column, sep=args.sep, encoding=args.encoding
        )
    except (FileNotFoundError, ValueError) as exc:
        logger.error("%s", exc)
        return 1

    doc_df = cl.get_documents(ids, chunk_size=args.chunk_size)
    if doc_df.empty or "pubmed_id" not in doc_df:
        try:
            doc_df.to_csv(args.output_csv, index=False)
            logger.info("Wrote %d rows to %s", len(doc_df), args.output_csv)
            return 0
        except OSError as exc:
            logger.error("failed to write output CSV: %s", exc)
            return 1

    # Normalise PubMed identifiers to strings to avoid dtype mismatches
    pubmed_ids = pd.to_numeric(doc_df["pubmed_id"], errors="coerce").astype("Int64")
    pmids = pubmed_ids.dropna().astype(str).tolist()
    pub_df = fetch_pubmed_records(pmids, args.sleep, args.workers, args.batch_size)
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
    try:
        merged.to_csv(args.output_csv, index=False)
        logger.info("Wrote %d rows to %s", len(merged), args.output_csv)
        return 0
    except OSError as exc:
        logger.error("failed to write output CSV: %s", exc)
        return 1


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Document data utilities")
    parser.add_argument("--log-level", default="INFO", help="Logging level")
    sub = parser.add_subparsers(dest="command", required=True)

    pubmed = sub.add_parser("pubmed", help="Fetch data from PubMed and related APIs")
    pubmed.add_argument("input_csv", type=Path, help="CSV with a PMID column")
    pubmed.add_argument("output_csv", type=Path, help="Destination CSV file")
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
        "input_csv", type=Path, help="CSV with document_chembl_id column"
    )
    chembl.add_argument("output_csv", type=Path, help="Destination CSV file")
    chembl.add_argument(
        "--column", default="document_chembl_id", help="Column name containing identifiers"
    )
    chembl.add_argument("--sep", default=",", help="CSV delimiter")
    chembl.add_argument("--encoding", default="utf8", help="File encoding")
    chembl.add_argument(
        "--chunk-size", type=int, default=5, help="Maximum number of IDs per request"
    )
    chembl.set_defaults(func=run_chembl)

    all_cmd = sub.add_parser("all", help="Run both ChEMBL and PubMed pipelines")
    all_cmd.add_argument(
        "input_csv", type=Path, help="CSV with document_chembl_id column"
    )
    all_cmd.add_argument("output_csv", type=Path, help="Destination CSV file")
    all_cmd.add_argument(
        "--column", default="document_chembl_id", help="Column in the input CSV"
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
    """Command line entry point."""

    parser = build_parser()
    args = parser.parse_args(argv)
    logging.basicConfig(level=getattr(logging, args.log_level.upper(), logging.INFO))
    return args.func(args)


if __name__ == "__main__":  # pragma: no cover - CLI entry point
    raise SystemExit(main())
