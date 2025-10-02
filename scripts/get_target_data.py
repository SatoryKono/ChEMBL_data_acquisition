"""Command line interface for retrieving target data from external sources.

Example
-------
Fetch ChEMBL target information for identifiers in ``targets.csv``::

    python scripts/get_target_data.py chembl --config config/config.yaml --input targets.csv
"""

from __future__ import annotations

# ruff: noqa: E402
import argparse
import sys
import shutil
from collections.abc import Iterator, Sequence
from contextlib import contextmanager
from dataclasses import dataclass
from functools import partial
from itertools import islice
from pathlib import Path
from typing import Any, cast

from datetime import datetime, timezone

from library.utils.bootstrap import ensure_project_root


if __package__ in {None, ""}:
    ensure_project_root()

import pandas as pd
import requests
from pandera.errors import SchemaErrors

import library.cli_utils as cli_utils_module
from library import chembl_library as cl
from library import cli
from library import io
from library import iuphar_library as ii
from library import protein_classification as pc
from library import target_postprocessing as tp
from library import uniprot_library as uu
from library.clients import ChemblClient
from library.cli_utils import PipelineError, run_cli_command, run_pipeline
from library.chembl_target import normalize_reaction_ec_numbers
from library.cli import (
    LoggerConfig,
    build_root_parser,
    configure_logger,
    path_argument,
    positive_int,
    prepare_io_paths,
)
from library.config import (
    Config,
    _serialize_paths,
)
from library.csv_utils import write_csv_deterministic
from library.log import logger
from library.metadata import Stats, file_sha256, write_meta_yaml
from library.pipeline_metadata import add_pipeline_metadata
from library.sidecar import SidecarErrors
from library.table_quality import analyze_table_quality
from library.validation import ValidationResult
from schemas import TargetsSchema, normalize_targets
from schemas.targets import TARGETS_COLUMN_ORDER


@contextmanager
def _override_cli_meta_writer() -> Iterator[None]:
    """Temporarily patch CLI metadata writer used by ``run_pipeline``."""

    original_cli_write_meta = cli_utils_module.write_meta_yaml
    cli_utils_module.write_meta_yaml = write_meta_yaml
    try:
        yield
    finally:
        cli_utils_module.write_meta_yaml = original_cli_write_meta


def _run_pipeline_with_meta(**kwargs: object) -> int:
    """Invoke :func:`run_pipeline` with project-specific metadata writer."""

    with _override_cli_meta_writer():
        return run_pipeline(**kwargs)

TARGETS_REQUIRED_COLUMNS: set[str] = {
    name for name, column in TargetsSchema.columns.items() if column.required
}

TARGETS_OPTIONAL_COLUMNS: list[str] = [
    column for column in TARGETS_COLUMN_ORDER if column not in TARGETS_REQUIRED_COLUMNS
]

TARGETS_OBJECT_COLUMNS: set[str] = {
    name
    for name, column in TargetsSchema.columns.items()
    if str(column.dtype) == "object"
}

UNIPROT_MISSING_VALUE = ""


DEFAULT_INPUT_NAME = "target.csv"
DEFAULT_OUTPUT_STEM = "targets"
RAW_SUFFIX = "_raw"
NORMALIZED_SUFFIX = "_normalized"


def _run_pipeline_with_meta(**kwargs: object) -> int:
    """Execute :func:`run_pipeline` with metadata writer override."""

    original_cli_write_meta = cli_utils_module.write_meta_yaml
    cli_utils_module.write_meta_yaml = write_meta_yaml
    try:
        return run_pipeline(**kwargs)
    finally:
        cli_utils_module.write_meta_yaml = original_cli_write_meta


@dataclass(frozen=True)
class _UniprotCandidate:
    """Container describing a UniProt identifier candidate for a target row."""

    value: str
    source: str
    original_id: str


@dataclass(frozen=True)
class _UniprotQueryPlan:
    """Deterministic mapping of ChEMBL rows to UniProt identifiers."""

    unique_records: list[dict[str, str]]
    row_candidates: list[list[_UniprotCandidate]]
    row_index: list[Any]


def _split_uniprot_tokens(value: str) -> Iterator[str]:
    """Yield cleaned UniProt identifiers from a pipe-delimited ``value``."""

    for token in value.split("|"):
        token = token.strip()
        if token:
            yield token


def _collect_uniprot_candidate_columns(df: pd.DataFrame, cfg: Config) -> list[str]:
    """Return ordered list of columns potentially holding UniProt accessions."""

    primary = cfg.target.all.uniprot_column
    ordered: list[str] = []
    if primary in df.columns:
        ordered.append(primary)

    preferred = ["uniprot_id", "mapping_uniprot_id"]
    for column in preferred:
        if column != primary and column in df.columns and column not in ordered:
            ordered.append(column)

    extra = [
        column
        for column in df.columns
        if column not in ordered
        and any(keyword in column.lower() for keyword in ("uniprot", "accession"))
    ]
    ordered.extend(extra)
    return ordered


def _build_uniprot_query_plan(df: pd.DataFrame, cfg: Config) -> _UniprotQueryPlan:
    """Create a deterministic plan for querying UniProt based on ``df``."""

    candidate_columns = _collect_uniprot_candidate_columns(df, cfg)
    unique_records: list[dict[str, str]] = []
    seen: set[str] = set()
    row_candidates: list[list[_UniprotCandidate]] = []
    row_index = list(df.index)

    if not candidate_columns:
        return _UniprotQueryPlan(unique_records, [[] for _ in row_index], row_index)

    positions = [df.columns.get_loc(column) for column in candidate_columns]
    for row in df.itertuples(index=False, name=None):
        candidates: list[_UniprotCandidate] = []
        row_seen: set[str] = set()
        for column, pos in zip(candidate_columns, positions):
            raw_value = row[pos]
            if not isinstance(raw_value, str) or not raw_value:
                continue
            for token in _split_uniprot_tokens(raw_value):
                if token in row_seen:
                    continue
                row_seen.add(token)
                candidate = _UniprotCandidate(
                    value=token, source=column, original_id=token
                )
                candidates.append(candidate)
                if token not in seen:
                    seen.add(token)
                    unique_records.append(
                        {"uniprot_id": candidate.value, "original_id": candidate.original_id}
                    )
        row_candidates.append(candidates)

    return _UniprotQueryPlan(unique_records, row_candidates, row_index)


def _resolve_uniprot_matches(
    plan: _UniprotQueryPlan, uniprot_df: pd.DataFrame
) -> pd.Series:
    """Return preferred UniProt identifier for each ChEMBL row in ``plan``."""

    if not plan.row_index:
        return pd.Series(dtype=object)

    lookup_column = "uniprot_id"
    if lookup_column not in uniprot_df.columns:
        return pd.Series(
            [UNIPROT_MISSING_VALUE for _ in plan.row_index],
            index=plan.row_index,
            dtype=object,
        )

    cleaned = uniprot_df[lookup_column].dropna().astype(str).map(str.strip)

    if "original_id" in uniprot_df.columns:
        original_series = (
            uniprot_df["original_id"].fillna("").astype(str).map(str.strip)
        )
        candidate_map: dict[str, str] = {}
        for canonical, original in zip(cleaned, original_series):
            if not canonical:
                continue
            candidate_map.setdefault(canonical, canonical)
            if original:
                for token in _split_uniprot_tokens(original):
                    if token and token not in candidate_map:
                        candidate_map[token] = canonical
    else:
        available = {value for value in cleaned if value}
        candidate_map = {value: value for value in available}
    resolved: list[str] = []
    for candidates in plan.row_candidates:
        match = UNIPROT_MISSING_VALUE
        for candidate in candidates:
            mapped = candidate_map.get(candidate.value)
            if mapped:
                match = mapped
                break
        resolved.append(match)

    return pd.Series(resolved, index=plan.row_index, dtype=object)


def _normalized_output_path(base: Path) -> Path:
    """Return deterministic path for the normalized export derived from ``base``."""

    suffix = base.suffix or ".csv"
    return base.with_name(f"{base.stem}{NORMALIZED_SUFFIX}{suffix}")


def _raw_output_path(base: Path) -> Path:
    """Return default path for the raw payload dump derived from ``base``."""

    suffix = base.suffix or ".csv"
    return base.with_name(f"{base.stem}{RAW_SUFFIX}{suffix}")


def _write_raw_dump(
    df: pd.DataFrame,
    destination: Path,
    *,
    cfg: Config,
    reindex_columns: bool,
) -> Path:
    """Persist ``df`` to ``destination`` honouring CSV or Parquet formats."""

    destination.parent.mkdir(parents=True, exist_ok=True)
    if reindex_columns and not df.empty:
        df = df.reindex(columns=sorted(df.columns))
    suffix = destination.suffix.lower()
    if suffix in {".parquet", ".pq"}:
        try:
            df.to_parquet(destination, index=False)
        except (ImportError, ValueError) as exc:
            raise OSError(f"failed to write parquet: {exc}") from exc
    else:
        df.to_csv(
            destination,
            index=False,
            sep=cfg.io.csv_sep,
            encoding=cfg.io.csv_encoding,
        )
    logger.info("raw_dump_written", rows=len(df), path=str(destination))
    return destination


def _pipe_merge(values: Sequence[str | None]) -> str:
    """Return a ``"|"``-joined string of unique, non-empty tokens.

    Parameters
    ----------
    values:
        Sequence of pipe-delimited strings to merge.

    Returns
    -------
    str
        Sorted, unique tokens separated by ``"|"``. Empty inputs yield
        an empty string.

    """
    tokens: set[str] = set()
    for value in values:
        if isinstance(value, str) and value:
            parts = [p.strip() for p in value.split("|") if p.strip()]
            tokens.update(parts)
    return "|".join(sorted(tokens))


def _first_token(value: str | None) -> str:
    """Return the first token from a pipe-delimited ``value``."""
    if isinstance(value, str) and value:
        return value.split("|")[0]
    return ""


def _prefer_primary(
    primary: pd.Series | None, secondary: pd.Series | None
) -> pd.Series:
    """Return a series prioritising ``primary`` values over ``secondary``.

    The function coalesces two series representing the same logical column
    coming from different data sources. Values from ``primary`` are preferred
    unless they are missing or empty, in which case the corresponding entry
    from ``secondary`` is used. Missing inputs yield an empty object-typed
    series to maintain downstream compatibility.
    """

    if primary is None and secondary is None:
        return pd.Series(dtype=object)

    if primary is None:
        return secondary.astype(object) if secondary is not None else pd.Series(dtype=object)

    if secondary is None:
        return primary.astype(object)

    primary = primary.astype(object)
    secondary = secondary.astype(object)
    result = primary.copy()
    if len(result) != len(secondary):
        secondary = secondary.reindex(result.index)
    mask = result.isna() | (result == "")
    if mask.any():
        result.loc[mask] = secondary.loc[mask]
    return result


def _prepare_targets_for_schema(
    df: pd.DataFrame,
) -> tuple[pd.DataFrame, set[str], set[str]]:
    """Return a copy of ``df`` aligned to :data:`TargetsSchema`.

    The returned frame contains only columns defined in
    :data:`TARGETS_COLUMN_ORDER`. Missing optional columns are created with
    ``"-"`` placeholders and schema fields expecting a generic ``object`` dtype
    are coerced accordingly.

    Returns
    -------
    tuple[pandas.DataFrame, set[str], set[str]]
        The prepared dataframe, names of missing required columns and optional
        columns that were injected.
    """

    missing_required = TARGETS_REQUIRED_COLUMNS - set(df.columns)
    missing_optional = {
        column for column in TARGETS_OPTIONAL_COLUMNS if column not in df.columns
    }

    if missing_optional:
        fill_values = {
            column: (
                pd.Series(["-"] * len(df), index=df.index, dtype=object)
                if len(df)
                else pd.Series(dtype=object)
            )
            for column in missing_optional
        }
        prepared = df.assign(**fill_values)
    else:
        prepared = df.copy()

    prepared = prepared.reindex(columns=TARGETS_COLUMN_ORDER)
    for column in TARGETS_OBJECT_COLUMNS & set(prepared.columns):
        prepared[column] = prepared[column].astype(object)
    return prepared, missing_required, missing_optional


def _save_snapshot(df: pd.DataFrame, base: Path, step: str, cfg: Config) -> Path:
    """Write ``df`` to a uniquely named snapshot CSV file with metadata.

    The file is created alongside ``base`` using the pattern
    ``<base>_<step>_<n>.csv`` where ``n`` increments to avoid overwriting
    existing files. Snapshot exports respect the configured CSV separator and
    encoding, and a ``.meta.yaml`` sidecar is generated containing minimal
    provenance details for reproducibility.

    Parameters
    ----------
    df:
        Data frame to serialise.
    base:
        Base path for the output file. Its stem and suffix determine the
        snapshot file name.
    step:
        Descriptive label inserted into the snapshot file name.
    cfg:
        Application configuration providing CSV export options.

    Returns
    -------
    Path
        Path to the written snapshot file.
    """
    stem = base.stem
    suffix = base.suffix or ".csv"
    index = 1
    while True:
        candidate = base.with_name(f"{stem}_{step}_{index}{suffix}")
        if not candidate.exists():
            work = df.copy()
            if work.columns.empty:
                key_cols: list[str] = []
            else:
                key_cols = [
                    column
                    for column in TARGETS_COLUMN_ORDER
                    if column in work.columns
                ]
                if not key_cols:
                    key_cols = list(work.columns)
            csv_path = write_csv_deterministic(
                work,
                candidate,
                col_order=list(df.columns),
                key_cols=key_cols,
                sep=cfg.io.csv_sep,
                encoding=cfg.io.csv_encoding,
                cfg=None,
            )
            stats: Stats = {
                "rows_total": len(df),
                "rows_kept": len(df),
                "rows_dropped": 0,
                "output_sha256": file_sha256(csv_path),
            }
            write_meta_yaml(
                csv_path=csv_path,
                command=" ".join(sys.argv),
                config_subset=_serialize_paths(cfg.to_dict()),
                inputs={"base": str(base), "step": step},
                stats=stats,
                schema="TargetSnapshot",
            )
            return csv_path
        index += 1


def build_parser() -> tuple[argparse.ArgumentParser, LoggerConfig]:
    """Create and return the top-level CLI argument parser.

    The command line interface is organised into sub-commands for retrieving
    data from individual sources (UniProt, ChEMBL and IUPHAR) as well as a
    convenience ``all`` command that runs all pipelines and merges their
    outputs.

    Returns
    -------
    tuple[argparse.ArgumentParser, LoggerConfig]
        Parser populated with every sub-command alongside the logging
        configuration used by :func:`main`.
    """

    root, shared, log_cfg = build_root_parser()

    def _add_output_arguments(
        parser_obj: argparse.ArgumentParser, *, defaults: bool
    ) -> None:
        final_default: object = None if defaults else argparse.SUPPRESS
        raw_default: object = None if defaults else argparse.SUPPRESS
        raw_format_default: object = "csv" if defaults else argparse.SUPPRESS
        id_cols_default: object = None if defaults else argparse.SUPPRESS
        no_reindex_default: object = False if defaults else argparse.SUPPRESS
        normalize_default: object = False if defaults else argparse.SUPPRESS
        parser_obj.add_argument(
            "--final-out",
            dest="final_out",
            type=path_argument,
            default=final_default,
            help=(
                "Destination for the validated target export "
                "(default: derived from input name)"
            ),
        )
        parser_obj.add_argument(
            "--raw-out",
            dest="raw_out",
            type=path_argument,
            default=raw_default,
            help=(
                "Location for the raw combined dataset prior to final "
                "normalisation"
            ),
        )
        parser_obj.add_argument(
            "--raw-format",
            dest="raw_format",
            choices=("csv", "parquet"),
            default=raw_format_default,
            help="Format used when writing the raw dataset (default: csv)",
        )
        parser_obj.add_argument(
            "--id-cols",
            dest="id_cols",
            nargs="+",
            default=id_cols_default,
            help="Identifier columns used for deterministic ordering",
        )
        parser_obj.add_argument(
            "--no-reindex-raw",
            dest="no_reindex_raw",
            action="store_true",
            default=no_reindex_default,
            help="Skip column reindexing when exporting the raw dataset",
        )
        parser_obj.add_argument(
            "--normalize-at-export",
            dest="normalize_at_export",
            action=argparse.BooleanOptionalAction,
            default=normalize_default,
            help=(
                "Apply normalisation immediately before writing the final output. "
                "Use --no-normalize-at-export to keep the raw payload."
            ),
        )
        has_out_alias = any(
            "--out" in action.option_strings for action in parser_obj._actions
        )
        if not has_out_alias:
            parser_obj.add_argument(
                "--out",
                dest="final_out_alias",
                type=path_argument,
                default=argparse.SUPPRESS,
                help=argparse.SUPPRESS,
            )

    _add_output_arguments(root, defaults=True)
    _add_output_arguments(shared, defaults=False)

    root.set_defaults(input_csv=Path(DEFAULT_INPUT_NAME))
    parser = argparse.ArgumentParser(
        description="Target data utilities", parents=[root]
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    # ----------------------------
    # UniProt sub-command
    # ----------------------------
    uniprot = subparsers.add_parser(
        "uniprot",
        parents=[shared],
        help="Extract information for UniProt accessions",
    )
    uniprot.add_argument(
        "--column",
        default="uniprot_id",
        choices=["uniprot_id", "mapping_uniprot_id"],
        help="Column in the input CSV containing UniProt accessions",
    )
    uniprot.add_argument(
        "--data-dir",
        type=path_argument,
        default=None,
        help=(
            "Directory containing '<uniprot_id>.json' files "
            "(default: config resources.uniprot_data_dir)"
        ),
    )
    uniprot.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Maximum number of identifiers to process",
    )
    uniprot.add_argument(
        "--offset",
        type=int,
        default=0,
        help="Number of identifiers to skip before processing",
    )
    uniprot.set_defaults(func=run_uniprot)

    # ----------------------------
    # ChEMBL sub-command
    # ----------------------------
    chembl = subparsers.add_parser(
        "chembl",
        parents=[shared],
        help="Retrieve target information from ChEMBL",
        conflict_handler="resolve",
    )
    # ``--id-cols`` and other export options are inherited from ``shared``.
    chembl.set_defaults(normalize_at_export=True)
    chembl.add_argument(
        "--column",
        default="target_chembl_id",
        help="Column name in the input CSV containing identifiers",
    )
    chembl.add_argument(
        "--chunk-size",
        type=positive_int,
        default=5,
        help="Maximum number of identifiers to request per call",
    )
    chembl.add_argument(
        "--timeout",
        type=float,
        default=30.0,
        help="Timeout in seconds for each HTTP request",
    )
    chembl.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Maximum number of identifiers to process",
    )
    chembl.add_argument(
        "--offset",
        type=int,
        default=0,
        help="Number of identifiers to skip before processing",
    )

    chembl.set_defaults(normalize_at_export=True)

    chembl.set_defaults(func=run_chembl)

    # ----------------------------
    # IUPHAR sub-command
    # ----------------------------
    iuphar = subparsers.add_parser(
        "iuphar",
        parents=[shared],
        help="Map UniProt accessions to IUPHAR classifications",
    )
    iuphar.add_argument(
        "--target-csv",
        type=path_argument,
        default=None,
        help=(
            "Path to the _IUPHAR_target.csv file "
            "(default: config resources.iuphar_target_csv)"
        ),
    )
    iuphar.add_argument(
        "--family-csv",
        type=path_argument,
        default=None,
        help=(
            "Path to the _IUPHAR_family.csv file "
            "(default: config resources.iuphar_family_csv)"
        ),
    )
    iuphar.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Maximum number of rows to process",
    )
    iuphar.add_argument(
        "--offset",
        type=int,
        default=0,
        help="Number of identifiers to skip before processing",
    )
    iuphar.set_defaults(func=run_iuphar)

    # ----------------------------
    # Combined pipeline
    # ----------------------------
    all_cmd = subparsers.add_parser(
        "all",
        parents=[shared],
        help="Run ChEMBL, UniProt and IUPHAR pipelines and merge results",
    )
    all_cmd.add_argument(
        "--chembl-out",
        dest="chembl_out",
        type=path_argument,
        help="Optional path to save intermediate ChEMBL data",
    )
    all_cmd.add_argument(
        "--uniprot-out",
        dest="uniprot_out",
        type=path_argument,
        help="Optional path to save intermediate UniProt data",
    )
    all_cmd.add_argument(
        "--iuphar-out",
        dest="iuphar_out",
        type=path_argument,
        help="Optional path to save intermediate IUPHAR data",
    )
    all_cmd.add_argument(
        "--data-dir",
        type=path_argument,
        default=None,
        help=(
            "Directory containing '<uniprot_id>.json' files "
            "(default: config resources.uniprot_data_dir)"
        ),
    )
    all_cmd.add_argument(
        "--target-csv",
        type=path_argument,
        default=None,
        help=(
            "Path to the _IUPHAR_target.csv file "
            "(default: config resources.iuphar_target_csv)"
        ),
    )
    all_cmd.add_argument(
        "--family-csv",
        type=path_argument,
        default=None,
        help=(
            "Path to the _IUPHAR_family.csv file "
            "(default: config resources.iuphar_family_csv)"
        ),
    )
    all_cmd.add_argument(
        "--timeout",
        type=float,
        default=30.0,
        help="Timeout in seconds for each HTTP request",
    )
    all_cmd.add_argument(
        "--uniprot-column",
        default="uniprot_id",
        choices=["uniprot_id", "mapping_uniprot_id"],
        help="Column from ChEMBL output to use for UniProt processing",
    )
    all_cmd.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Maximum number of identifiers to process",
    )
    all_cmd.add_argument(
        "--offset",
        type=int,
        default=0,
        help="Number of identifiers to skip before processing",
    )
    all_cmd.set_defaults(func=run_all)

    parser.subparsers_map = {  # type: ignore[attr-defined]
        "uniprot": uniprot,
        "chembl": chembl,
        "iuphar": iuphar,
        "all": all_cmd,
    }

    return parser, log_cfg


def run_uniprot(cfg: Config, args: argparse.Namespace) -> int:
    """Execute the ``uniprot`` sub-command.

    Parameters
    ----------
    cfg : Config
        Application configuration containing UniProt-specific overrides and
        shared IO settings.
    args : argparse.Namespace
        Parsed command-line arguments produced by :func:`build_parser` for the
        ``uniprot`` sub-command.

    Returns
    -------
    int
        ``0`` on success. Non-zero values indicate failures while reading the
        input CSV, interacting with the UniProt data directory or producing the
        derived artefacts. Input validation errors are logged and converted into
        a failure code.
    """
    limit = cfg.target.uniprot.limit
    if limit is not None and limit < 0:
        logger.error(
            "invalid_limit",
            section="target.uniprot.limit",
            limit=limit,
        )
        return 1

    try:
        df = pd.read_csv(
            args.input_csv, sep=cfg.io.csv_sep, encoding=cfg.io.csv_encoding, dtype=str
        )
        column = cfg.target.uniprot.column
        if column not in df.columns:
            raise ValueError(f"column '{column}' not found in {args.input_csv}")
        df = df.fillna("")
        df = df[(df[column].str.strip() != "") & (df[column] != "#N/A")].reset_index(
            drop=True
        )
        offset = getattr(args, "offset", 0)
        if offset:
            original_rows = len(df)
            df = df.iloc[offset:].reset_index(drop=True)
            logger.info("process_offset", offset=min(offset, original_rows))
        if limit is not None:
            df = df.head(limit)
            logger.info("process_limit", limit=len(df))
        ids = df[column].to_numpy(copy=False)
        rows_total = len(ids)

        from tempfile import NamedTemporaryFile

        with NamedTemporaryFile(
            "w", delete=False, encoding=cfg.io.csv_encoding, newline=""
        ) as tmp:
            tmp_path = Path(tmp.name)

        pd.DataFrame({"uniprot_id": ids}).to_csv(
            tmp_path,
            index=False,
            sep=cfg.io.csv_sep,
            encoding=cfg.io.csv_encoding,
        )

        output_candidate = getattr(args, "output_csv", None)
        output = output_candidate or io.default_output_path(args.input_csv, cfg.io)
        data_dir = cfg.target.uniprot.data_dir
        uu.init_session(cfg.api, cfg.retry)
        try:
            uu.process(
                input_csv=str(tmp_path),
                output_csv=str(output),
                data_dir=data_dir,
                cfg=cfg.uniprot,
                gtop_cfg=cfg.iuphar,
                sep=cfg.io.csv_sep,
                encoding=cfg.io.csv_encoding,
            )
        finally:
            tmp_path.unlink(missing_ok=True)

        out_df = pd.read_csv(
            output, sep=cfg.io.csv_sep, encoding=cfg.io.csv_encoding, dtype=str
        )
        rows_kept = len(out_df)
        if "mapping_uniprot_id" in df.columns:
            out_df.insert(
                1,
                "mapping_uniprot_id",
                df["mapping_uniprot_id"].to_numpy(copy=False),
            )
        csv_path = io.write_csv(
            out_df,
            output,
            cfg=cfg,
            sep=cfg.io.csv_sep,
            encoding=cfg.io.csv_encoding,
            key_cols=["uniprot_id"],
        )
        rows_dropped = max(rows_total - rows_kept, 0)
        stats: Stats = {
            "rows_total": rows_total,
            "rows_kept": rows_kept,
            "rows_dropped": rows_dropped,
            "output_sha256": file_sha256(csv_path),
        }
        inputs = {"input_csv": str(args.input_csv)}
        if cfg.target.uniprot.data_dir:
            inputs["data_dir"] = str(cfg.target.uniprot.data_dir)
        write_meta_yaml(
            csv_path=csv_path,
            command=" ".join(sys.argv),
            config_subset=_serialize_paths(cfg.to_dict()),
            inputs=inputs,
            stats=stats,
            schema="UniProtExport",
        )
    except (FileNotFoundError, ValueError, OSError) as exc:
        logger.error(
            "uniprot_processing_failed",
            error=str(exc),
            input=str(args.input_csv),
            output=str(output),
        )
        return 1
    try:
        analyze_table_quality(out_df, table_name=str(output.with_suffix("")))
    except ValueError as exc:
        logger.error(
            "quality_report_failed",
            error=str(exc),
            path=str(output),
        )
        return 1
    return 0


def _limited_ids(ids_iter: Iterator[str], limit: int) -> Iterator[str]:
    """Yield up to ``limit`` identifiers while logging the processed count."""

    count = 0
    try:
        for target_id in islice(ids_iter, limit):
            count += 1
            yield target_id
    finally:
        logger.info("process_limit", limit=count)


def run_chembl(cfg: Config, args: argparse.Namespace) -> int:
    """Execute the ``chembl`` sub-command.

    Parameters
    ----------
    cfg : Config
        Application configuration providing ChEMBL client settings, retry
        policy and CSV export behaviour.
    args : argparse.Namespace
        Parsed command-line arguments produced by :func:`build_parser`.

    Returns
    -------
    int
        ``0`` on success. Non-zero values indicate that reading identifiers,
        fetching data from the ChEMBL API, validating the result or writing the
        CSV artefact failed. All errors are logged with structured context.
    """
    limit = cfg.target.chembl.limit
    if limit is not None and limit < 0:
        logger.error(
            "invalid_limit",
            section="target.chembl.limit",
            limit=limit,
        )
        return 1

    try:
        ids_iter = io.read_ids(
            args.input_csv, column=cfg.target.chembl.column, cfg=cfg.io
        )
    except (FileNotFoundError, ValueError) as exc:
        logger.error(
            "read_fail",
            error=str(exc),
            path=str(args.input_csv),
        )
        return 1

    offset = getattr(args, "offset", 0)
    if offset:
        ids_iter = islice(ids_iter, offset, None)
        logger.info("process_offset", offset=offset)
    limited_ids = list(islice(ids_iter, limit)) if limit is not None else list(ids_iter)
    processed_ids = len(limited_ids)

    output_candidate = getattr(args, "output_csv", None)
    base_output = Path(
        output_candidate or io.default_output_path(args.input_csv, cfg.io)
    )

    raw_candidate = getattr(args, "raw_out", None)
    if raw_candidate in (None, argparse.SUPPRESS):
        legacy_raw_candidate = getattr(args, "raw_output", None)
        if legacy_raw_candidate not in (None, argparse.SUPPRESS):
            raw_candidate = legacy_raw_candidate
    raw_path_override = (
        Path(raw_candidate)
        if raw_candidate not in (None, argparse.SUPPRESS)
        else None
    )
    raw_destination = raw_path_override or _raw_output_path(base_output)
    normalized_output = _normalized_output_path(base_output)


    normalize_flag = getattr(args, "normalize_at_export", None)
    if normalize_flag in (None, argparse.SUPPRESS):
        normalize_at_export = True
    else:
        normalize_at_export = bool(normalize_flag)
    reindex_raw = not bool(getattr(args, "no_reindex_raw", False))


    raw_chunks: list[pd.DataFrame] = []
    parsed_chunks: list[pd.DataFrame] = []

    if limited_ids:
        try:
            with ChemblClient(cfg.api, cfg.retry, cfg.chembl) as client:
                for _, raw_chunk, parsed_chunk in cl.iter_target_batches(
                    limited_ids,
                    cfg=cfg.api,
                    client=client,
                    mapping_cfg=cfg.uniprot_mapping,
                    chunk_size=cfg.target.chembl.chunk_size,
                    timeout=cfg.target.chembl.timeout,
                ):
                    if not raw_chunk.empty:
                        raw_chunks.append(raw_chunk)
                    if not parsed_chunk.empty:
                        parsed_chunks.append(parsed_chunk)
        except (requests.RequestException, ValueError) as exc:
            logger.error(
                "chembl_fetch_failed",
                error=str(exc),
                chunk_size=cfg.target.chembl.chunk_size,
                timeout=cfg.target.chembl.timeout,
            )
            return 1

    raw_df = (
        pd.concat(raw_chunks, ignore_index=True) if raw_chunks else pd.DataFrame()
    )

    if normalize_at_export:
        try:
            _write_raw_dump(
                raw_df, raw_destination, cfg=cfg, reindex_columns=reindex_raw
            )
        except (OSError, ValueError) as exc:
            logger.error("raw_dump_failed", error=str(exc), path=str(raw_destination))
            return 1
    else:
        def _raw_fetcher() -> Iterator[pd.DataFrame]:
            yield raw_df

        def _raw_writer(
            chunks: Iterator[pd.DataFrame],
            destination: Path,
            col_order: Sequence[str] | None,
            key_cols: Sequence[str],
        ) -> Path:
            frames = list(chunks)
            if frames:
                merged = pd.concat(frames, ignore_index=True)
            else:
                merged = pd.DataFrame()
            _write_raw_dump(merged, destination, cfg=cfg, reindex_columns=reindex_raw)
            return destination

        failure_path = raw_destination.with_name(
            f"{raw_destination.stem}_failure_cases.csv"
        )
        exit_code = _run_pipeline_with_meta(
            fetcher=_raw_fetcher,
            schema=None,
            schema_name="raw_target_payload",
            validators=[],
            metadata_hooks=[],
            writer=_raw_writer,
            output_path=raw_destination,
            failure_path=failure_path,
            command=" ".join(sys.argv),
            config_snapshot=_serialize_paths(cfg.to_dict()),
            inputs={"input_csv": str(args.input_csv)},
            key_columns=[],
            table_quality=lambda _: None,
            cfg=cfg,
            logger=logger,
        )
        if exit_code != 0:
            return exit_code

        if not raw_destination.exists():
            logger.error("raw_dump_missing", path=str(raw_destination))
            return 1

        destination = base_output
        if raw_destination != destination:
            destination.parent.mkdir(parents=True, exist_ok=True)
            try:
                shutil.copy2(raw_destination, destination)
            except OSError as exc:
                logger.error(
                    "raw_to_final_copy_failed",
                    error=str(exc),
                    source=str(raw_destination),
                    destination=str(destination),
                )
                return 1

        if limit is not None:
            logger.info("process_limit", limit=processed_ids)
        logger.info(
            "chembl_normalization_skipped",
            output=str(destination),
            rows=len(raw_df),
        )
        return 0


    final_candidate = getattr(args, "final_out", None)
    if final_candidate in (None, argparse.SUPPRESS):
        output_candidate = getattr(args, "output_csv", None)
        if output_candidate not in (None, argparse.SUPPRESS):
            final_output = Path(output_candidate)
        else:
            final_output = Path(io.default_output_path(args.input_csv, cfg.io))
    else:
        final_output = Path(final_candidate)

    if raw_path_override is None:
        raw_output = final_output
    else:
        raw_output = raw_path_override

    raw_format = str(getattr(args, "raw_format", "csv") or "csv").lower()
    if raw_format not in {"csv", "parquet"}:
        logger.warning(
            "unsupported_raw_format", format=raw_format, fallback="csv"
        )
        raw_format = "csv"

    reindex_raw = not bool(getattr(args, "no_reindex_raw", False))
    # ``normalize_at_export`` already coerced above to avoid divergent values.

    id_cols_value = getattr(args, "id_cols", None)
    if id_cols_value in (None, argparse.SUPPRESS):
        key_columns = ["target_chembl_id"]
    elif isinstance(id_cols_value, str):
        key_columns = [id_cols_value]
    else:
        key_columns = list(id_cols_value) or ["target_chembl_id"]
    failure_path = final_output.with_name(f"{final_output.stem}_failure_cases.csv")


    missing_optional_columns: set[str] = set()
    placeholder_replacements = 0
    fetched_rows_total = 0
    raw_dump_rows_total = 0
    post_cleanup_rows_total = 0

    def _prepare_chunk(frame: pd.DataFrame) -> pd.DataFrame:
        nonlocal placeholder_replacements, post_cleanup_rows_total
        prepared, _, missing_optional = _prepare_targets_for_schema(frame)
        if missing_optional and not frame.empty:
            placeholder_replacements += len(frame) * len(missing_optional)
        post_cleanup_rows_total += len(prepared)
        missing_optional_columns.update(missing_optional)
        return prepared

    def _normalize_chunk(frame: pd.DataFrame) -> pd.DataFrame:
        nonlocal raw_dump_rows_total
        normalized_chunk = normalize_targets(frame)
        raw_dump_rows_total += len(normalized_chunk)
        return normalized_chunk

    def _validate_chunk(frame: pd.DataFrame) -> ValidationResult:
        try:
            validated = TargetsSchema.validate(frame, lazy=True)
        except SchemaErrors as exc:
            validated_subset = getattr(exc, "validated_data", frame)
            return ValidationResult(
                validated_subset,
                exc.failure_cases.copy(),
                "TargetsSchema",
            )
        return ValidationResult(validated, pd.DataFrame(), "TargetsSchema")

    def fetcher() -> Iterator[pd.DataFrame]:

        nonlocal fetched_rows_total
        if parsed_chunks:
            for chunk_df in parsed_chunks:
                fetched_rows_total += len(chunk_df)
                yield chunk_df
            return

        if not limited_ids:
            yield pd.DataFrame(columns=TARGETS_COLUMN_ORDER)
            return

        with ChemblClient(cfg.api, cfg.retry, cfg.chembl) as client:
            try:
                for _, _, parsed_chunk in cl.iter_target_batches(
                    limited_ids,
                    cfg=cfg.api,
                    client=client,
                    mapping_cfg=cfg.uniprot_mapping,
                    chunk_size=cfg.target.chembl.chunk_size,
                    timeout=cfg.target.chembl.timeout,
                ):
                    if parsed_chunk.empty:
                        continue
                    fetched_rows_total += len(parsed_chunk)
                    yield parsed_chunk
            except (requests.RequestException, ValueError) as exc:
                logger.error(
                    "chembl_fetch_failed",
                    error=str(exc),
                    chunk_size=cfg.target.chembl.chunk_size,
                    timeout=cfg.target.chembl.timeout,
                )
                raise PipelineError(str(exc)) from exc


    def writer(
        chunks: Iterator[pd.DataFrame],
        destination: Path,
        col_order: Sequence[str] | None,
        key_cols: Sequence[str],
    ) -> Path:
        resolved_keys = list(key_cols) if key_cols else key_columns
        column_order: Sequence[str] | None = col_order if reindex_raw else None

        if raw_format == "csv" and not normalize_at_export:
            raw_path = io.write_csv(
                chunks,
                raw_output,
                cfg=cfg,
                sep=cfg.io.csv_sep,
                encoding=cfg.io.csv_encoding,
                key_cols=resolved_keys or None,
                col_order=column_order,
                chunksize=cfg.io.csv_chunksize,
            )
            if final_output != raw_output:
                shutil.copy2(raw_path, final_output)
                return final_output
            return raw_path

        frames: list[pd.DataFrame] = []
        for chunk in chunks:
            working = chunk
            if column_order is not None:
                working = working.reindex(columns=column_order)
            frames.append(working)

        if frames:
            combined = pd.concat(frames, ignore_index=True)
        else:
            combined = pd.DataFrame(columns=column_order or [])

        if resolved_keys:
            missing_keys = [col for col in resolved_keys if col not in combined.columns]
            if not missing_keys:
                combined = combined.sort_values(by=resolved_keys).reset_index(drop=True)

        if raw_format == "parquet":
            try:
                combined.to_parquet(raw_output, index=False)
            except ImportError as exc:  # pragma: no cover - optional dependency
                raise ValueError(
                    "Parquet export requires optional pyarrow or fastparquet"
                ) from exc
            raw_path = raw_output
        else:
            raw_path = io.write_csv(
                combined,
                raw_output,
                cfg=cfg,
                sep=cfg.io.csv_sep,
                encoding=cfg.io.csv_encoding,
                key_cols=resolved_keys or None,
                col_order=column_order,
            )

        final_df = combined
        if normalize_at_export:
            final_df = normalize_targets(final_df)
            final_df, _, missing_optional = _prepare_targets_for_schema(final_df)
            missing_optional_columns.update(missing_optional)

        if final_output == raw_output and not normalize_at_export and raw_format == "parquet":
            final_path = raw_path
        else:
            final_path = io.write_csv(
                final_df,
                final_output,
                cfg=cfg,
                sep=cfg.io.csv_sep,
                encoding=cfg.io.csv_encoding,
                key_cols=resolved_keys or None,
                col_order=TARGETS_COLUMN_ORDER,
            )
        return final_path

    failure_path = normalized_output.with_name(
        f"{normalized_output.stem}_failure_cases.csv"
    )
    table_quality = partial(
        analyze_table_quality,

        table_name=str(final_output.with_suffix("")),
    )

    metadata_hooks = [add_pipeline_metadata, _prepare_chunk]
    if not normalize_at_export:
        metadata_hooks.insert(0, normalize_targets)


    exit_code = _run_pipeline_with_meta(
        fetcher=fetcher,
        schema=TargetsSchema,
        schema_name="TargetsSchema",
        validators=[_validate_chunk],
        metadata_hooks=metadata_hooks,
        writer=writer,
        output_path=raw_output,
        failure_path=failure_path,
        command=" ".join(sys.argv),
        config_snapshot=_serialize_paths(cfg.to_dict()),
        inputs={"input_csv": str(args.input_csv)},
        key_columns=key_columns,
        table_quality=table_quality,
        cfg=cfg,
        logger=logger,
    )

    if limit is not None:
        logger.info("process_limit", limit=processed_ids)

    if missing_optional_columns:
        logger.debug(
            "schema_optional_columns_missing",
            columns=sorted(missing_optional_columns),
        )

    logger.info("chembl_stage_rows", stage="fetch", rows=fetched_rows_total)
    logger.info("chembl_stage_rows", stage="raw_dump", rows=raw_dump_rows_total)
    logger.info(
        "chembl_stage_rows",
        stage="post_cleanup",
        rows=post_cleanup_rows_total,
    )
    logger.info(
        "chembl_placeholder_replacements",
        total=placeholder_replacements,
    )

    return exit_code


def run_iuphar(cfg: Config, args: argparse.Namespace) -> int:
    """Execute the ``iuphar`` sub-command.

    Parameters
    ----------
    cfg : Config
        Application configuration containing paths to the IUPHAR export files
        and shared IO settings.
    args : argparse.Namespace
        Parsed command-line arguments produced by :func:`build_parser`.

    Returns
    -------
    int
        ``0`` on success. Non-zero values indicate that the IUPHAR sources
        could not be read, the combined dataset failed validation or the CSV
        output could not be written.
    """
    limit = cfg.target.iuphar.limit
    if limit is not None and limit < 0:
        logger.error(
            "invalid_limit",
            section="target.iuphar.limit",
            limit=limit,
        )
        return 1

    tmp_path: Path | None = None
    source_csv = args.input_csv

    try:
        df_to_process: pd.DataFrame | None = None
        offset = getattr(args, "offset", 0)
        if limit is not None or offset:
            df_full = pd.read_csv(
                args.input_csv,
                sep=cfg.io.csv_sep,
                encoding=cfg.io.csv_encoding,
                dtype=str,
            )
            if offset:
                total_rows = len(df_full)
                df_full = df_full.iloc[offset:].reset_index(drop=True)
                logger.info("process_offset", offset=min(offset, total_rows))
            if limit is not None:
                df_full = df_full.head(limit)
                logger.info("process_limit", limit=len(df_full))
            df_to_process = df_full
        if df_to_process is not None:
            from tempfile import NamedTemporaryFile

            with NamedTemporaryFile(delete=False) as tmp:
                tmp_path = Path(tmp.name)
            df_to_process.to_csv(
                tmp_path,
                index=False,
                sep=cfg.io.csv_sep,
                encoding=cfg.io.csv_encoding,
            )
            source_csv = tmp_path

        data = ii.IUPHARData.from_files(
            target_path=cfg.target.iuphar.target_csv,
            family_path=cfg.target.iuphar.family_csv,
            encoding=cfg.io.csv_encoding,
        )
        output_candidate = getattr(args, "output_csv", None)
        output = output_candidate or io.default_output_path(args.input_csv, cfg.io)
        data.map_uniprot_file(
            input_path=source_csv,
            output_path=output,
            encoding=cfg.io.csv_encoding,
            sep=cfg.io.csv_sep,
        )
    except (FileNotFoundError, ValueError, OSError) as exc:
        logger.error(
            "iuphar_processing_failed",
            error=str(exc),
            input=str(source_csv),
            target_csv=str(cfg.target.iuphar.target_csv),
            family_csv=str(cfg.target.iuphar.family_csv),
        )
        return 1
    finally:
        if tmp_path is not None:
            tmp_path.unlink(missing_ok=True)
    try:
        analyze_table_quality(output, table_name=str(output.with_suffix("")))
    except ValueError as exc:
        logger.error(
            "quality_report_failed",
            error=str(exc),
            path=str(output),
        )
        return 1
    return 0


def fetch_chembl(
    cfg: Config,
    input_csv: Path,
    final_out: Path,
    limit: int | None = None,
    *,
    raw_out: Path | None = None,
    raw_format: str = "csv",
    id_cols: Sequence[str] | None = None,
    chunk_size: int | None = None,
    offset: int = 0,
    normalize_at_export: bool = False,
    no_reindex_raw: bool = False,
) -> pd.DataFrame:
    """Fetch target information from ChEMBL.

    Parameters
    ----------
    cfg : Config
        Application configuration used to drive the ChEMBL pipeline.
    input_csv : pathlib.Path
        Source CSV containing target identifiers.
    final_out : pathlib.Path
        Destination path used by :func:`run_chembl` to persist results.
    limit : int, optional
        Maximum number of identifiers to process. ``None`` processes all rows.
    chunk_size : int, optional
        Temporary override for the batch size used when calling the API.
    offset : int, optional
        Number of identifiers to skip before starting the retrieval.

    Returns
    -------
    pandas.DataFrame
        Retrieved ChEMBL data loaded from ``final_out``.

    Raises
    ------
    RuntimeError
        Raised when :func:`run_chembl` reports a non-zero exit code.
    """

    logger.info("fetch_chembl_start", input=str(input_csv), output=str(final_out))
    chembl_args = argparse.Namespace(
        input_csv=input_csv,
        final_out=final_out,
        output_csv=final_out,
        raw_out=raw_out,
        raw_format=raw_format,
        id_cols=list(id_cols) if id_cols is not None else None,
        offset=offset,
        normalize_at_export=normalize_at_export,
        no_reindex_raw=no_reindex_raw,
    )
    original_limit = cfg.target.chembl.limit
    original_chunk_size = cfg.target.chembl.chunk_size
    if limit is not None:
        cfg.target.chembl.limit = limit
    if chunk_size is not None:
        cfg.target.chembl.chunk_size = chunk_size
    try:
        if run_chembl(cfg, chembl_args) != 0:
            raise RuntimeError("ChEMBL retrieval failed")
    finally:
        if limit is not None:
            cfg.target.chembl.limit = original_limit
        if chunk_size is not None:
            cfg.target.chembl.chunk_size = original_chunk_size
    normalized_output = _normalized_output_path(final_out)
    df = pd.read_csv(

        final_out, sep=cfg.io.csv_sep, encoding=cfg.io.csv_encoding, dtype=str

    )
    logger.info("fetch_chembl_done", rows=len(df))
    return df


def fetch_uniprot(
    cfg: Config, chembl_df: pd.DataFrame, output_csv: Path
) -> pd.DataFrame:
    """Retrieve UniProt annotations for targets in ``chembl_df``.

    Parameters
    ----------
    cfg : Config
        Application configuration used to invoke :func:`run_uniprot`.
    chembl_df : pandas.DataFrame
        DataFrame containing ChEMBL target records with at least one UniProt
        identifier column.
    output_csv : pathlib.Path
        Destination path populated by the UniProt pipeline.

    Returns
    -------
    pandas.DataFrame
        UniProt records with an additional ``original_id`` column preserving
        the queried accessions.

    Raises
    ------
    RuntimeError
        Raised when :func:`run_uniprot` returns a non-zero exit code.
    """

    logger.info("fetch_uniprot_start", output=str(output_csv))
    plan = _build_uniprot_query_plan(chembl_df, cfg)
    if plan.unique_records:
        id_df = pd.DataFrame(plan.unique_records, dtype=object)
    else:
        id_df = pd.DataFrame(columns=["uniprot_id", "original_id"], dtype=object)

    id_df = id_df.copy()
    id_df["__query_order"] = range(len(id_df))
    query_input_df = id_df.drop(columns=["__query_order"], errors="ignore")

    from tempfile import NamedTemporaryFile

    with NamedTemporaryFile(
        "w", delete=False, encoding=cfg.io.csv_encoding, newline=""
    ) as tmp:
        tmp_path = Path(tmp.name)

    query_input_df.to_csv(
        tmp_path,
        index=False,
        sep=cfg.io.csv_sep,
        encoding=cfg.io.csv_encoding,
    )

    uniprot_args = argparse.Namespace(input_csv=tmp_path, output_csv=output_csv)
    orig_dir = cfg.target.uniprot.data_dir
    cfg.target.uniprot.data_dir = cfg.target.all.data_dir
    try:
        if run_uniprot(cfg, uniprot_args) != 0:
            raise RuntimeError("UniProt retrieval failed")
    finally:
        cfg.target.uniprot.data_dir = orig_dir
        tmp_path.unlink(missing_ok=True)

    fetched_df = pd.read_csv(
        output_csv, sep=cfg.io.csv_sep, encoding=cfg.io.csv_encoding, dtype=str
    )
    if "uniprot_id" in fetched_df.columns:
        df = id_df.merge(fetched_df, on="uniprot_id", how="left", sort=False)
    else:
        df = id_df.copy()

    if "__query_order" in df.columns:
        df = df.sort_values("__query_order").drop(columns=["__query_order"], errors="ignore")
        df = df.reset_index(drop=True)

    left_original = df.pop("original_id_x") if "original_id_x" in df.columns else None
    right_original = df.pop("original_id_y") if "original_id_y" in df.columns else None

    if left_original is not None or right_original is not None:
        df["original_id"] = _prefer_primary(left_original, right_original)
    elif "original_id" not in df.columns:
        df["original_id"] = pd.Series(dtype=object)
    logger.info("fetch_uniprot_done", rows=len(df))
    return df


def fetch_iuphar(
    cfg: Config,
    chembl_df: pd.DataFrame,
    uniprot_df: pd.DataFrame,
    output_csv: Path,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Retrieve IUPHAR classifications for combined target data.

    Parameters
    ----------
    cfg : Config
        Application configuration containing IUPHAR resource paths.
    chembl_df : pandas.DataFrame
        ChEMBL target data obtained from :func:`fetch_chembl`.
    uniprot_df : pandas.DataFrame
        UniProt annotations with an ``original_id`` column provided by
        :func:`fetch_uniprot`.
    output_csv : pathlib.Path
        Destination where :func:`run_iuphar` persists the retrieved
        classifications.

    Returns
    -------
    tuple[pandas.DataFrame, pandas.DataFrame]
        Two data frames: the merged ChEMBL/UniProt input and the IUPHAR
        classification results.

    Raises
    ------
    RuntimeError
        Raised when :func:`run_iuphar` reports a non-zero exit code.
    """

    logger.info("fetch_iuphar_start", output=str(output_csv))
    merge_column = cfg.target.all.uniprot_column
    plan = _build_uniprot_query_plan(chembl_df, cfg)
    lookup_series = _resolve_uniprot_matches(plan, uniprot_df)
    if lookup_series.empty and len(chembl_df.index):
        lookup_series = pd.Series(
            [UNIPROT_MISSING_VALUE] * len(chembl_df.index),
            index=chembl_df.index,
            dtype=object,
        )
    else:
        lookup_series = lookup_series.reindex(
            chembl_df.index, fill_value=UNIPROT_MISSING_VALUE
        )
    lookup_column = "__uniprot_lookup_id"

    if merge_column != "uniprot_id":
        chembl_for_merge = chembl_df.drop(columns=["uniprot_id"], errors="ignore")
    else:
        chembl_for_merge = chembl_df.copy()

    if merge_column not in chembl_for_merge.columns:
        logger.warning("missing_uniprot_column", column=merge_column)
        placeholder = pd.Series(
            UNIPROT_MISSING_VALUE,
            index=chembl_for_merge.index,
            dtype=object,
        )
        chembl_for_merge = chembl_for_merge.assign(**{merge_column: placeholder})

    chembl_for_merge[lookup_column] = lookup_series

    combined_df = pd.merge(
        chembl_for_merge,
        uniprot_df,
        left_on=lookup_column,
        right_on="uniprot_id",
        how="left",
        suffixes=("_chembl", "_uniprot"),
    )

    merge_suffix_priority = ("_uniprot", "_chembl")

    def _base_name(column: str) -> str | None:
        for suffix in merge_suffix_priority:
            if column.endswith(suffix):
                return column[: -len(suffix)]
        return None

    def _coalesce_column(df: pd.DataFrame, column: str) -> None:
        preferred: pd.Series | None = None
        for suffix in merge_suffix_priority:
            candidate_name = f"{column}{suffix}"
            if candidate_name not in df.columns:
                continue
            candidate = df.pop(candidate_name).astype(object)
            if preferred is None:
                preferred = candidate
            else:
                preferred = _prefer_primary(preferred, candidate)
        if preferred is not None:
            df[column] = preferred

    critical_columns = set(TARGETS_COLUMN_ORDER)
    critical_columns.update(
        {
            "gene",
            "gene_name",
            "synonyms",
            "component_description",
            "chembl_alternative_name",
            "names",
            "mapping_uniprot_id",
            "ec_numbers",
            "reaction_ec_numbers",
            "ec_code",
        }
    )

    suffixed_columns = [
        column
        for column in combined_df.columns
        if _base_name(column) is not None
    ]
    base_columns = {_base_name(column) for column in suffixed_columns if _base_name(column)}

    for column in sorted(base_columns & critical_columns):
        _coalesce_column(combined_df, column)

    for column in sorted(base_columns - critical_columns):
        _coalesce_column(combined_df, column)

    remaining_suffixes = [
        column
        for column in combined_df.columns
        if _base_name(column) is not None
    ]
    if remaining_suffixes:
        rename_map = {
            column: _base_name(column)
            for column in remaining_suffixes
            if _base_name(column)
        }
        combined_df = combined_df.rename(columns=rename_map)

    combined_df = combined_df.drop(columns=[lookup_column], errors="ignore")
    if "original_id" in combined_df.columns:
        combined_df = combined_df.drop(columns=["original_id"])

    overlap_columns = sorted(
        set(chembl_for_merge.columns)
        & (set(uniprot_df.columns) - {"original_id"})
    )
    for column in overlap_columns:
        chembl_col = f"{column}_chembl"
        uniprot_col = f"{column}_uniprot"
        chembl_series = (
            combined_df.pop(chembl_col)
            if chembl_col in combined_df.columns
            else None
        )
        uniprot_series = (
            combined_df.pop(uniprot_col)
            if uniprot_col in combined_df.columns
            else None
        )
        if column == "reaction_ec_numbers":
            chembl_values = (
                chembl_series.reindex(combined_df.index)
                if chembl_series is not None
                else pd.Series(index=combined_df.index, dtype=object)
            )
            uniprot_values = (
                uniprot_series.reindex(combined_df.index)
                if uniprot_series is not None
                else pd.Series(index=combined_df.index, dtype=object)
            )
            combined_df[column] = pd.Series(
                [
                    normalize_reaction_ec_numbers([u, c])
                    for u, c in zip(uniprot_values, chembl_values)
                ],
                index=combined_df.index,
                dtype=object,
            )
        else:
            combined_df[column] = _prefer_primary(uniprot_series, chembl_series)

    if "gene" not in combined_df.columns:
        combined_df["gene"] = pd.Series(
            UNIPROT_MISSING_VALUE,
            index=combined_df.index,
            dtype=object,
        )

    ec_number_columns = [
        column for column in combined_df.columns if column.startswith("ec_numbers")
    ]
    if ec_number_columns:
        combined_df["ec_numbers"] = combined_df.apply(
            lambda r: _pipe_merge([r.get(column) for column in ec_number_columns]),
            axis=1,
        )
        extra_ec_columns = [
            column for column in ec_number_columns if column != "ec_numbers"
        ]
        if extra_ec_columns:
            combined_df = combined_df.drop(columns=extra_ec_columns, errors="ignore")

    reaction_ec_columns = [
        column
        for column in combined_df.columns
        if column.startswith("reaction_ec_numbers")
    ]
    if reaction_ec_columns:
        combined_df["reaction_ec_numbers"] = combined_df.apply(
            lambda r: normalize_reaction_ec_numbers(
                [r.get(column) for column in reaction_ec_columns]
            ),
            axis=1,
        )
        extra_reaction_columns = [
            column for column in reaction_ec_columns if column != "reaction_ec_numbers"
        ]
        if extra_reaction_columns:
            combined_df = combined_df.drop(
                columns=extra_reaction_columns, errors="ignore"
            )

    combined_df["synonyms"] = combined_df.apply(
        lambda r: _pipe_merge(
            [
                r.get("pref_name"),
                r.get("component_description"),
                r.get("gene"),
                r.get("chembl_alternative_name"),
                r.get("names"),
                r.get("secondaryAccessionNames"),
            ]
        ),
        axis=1,
    )
    combined_df["ec_number"] = combined_df.apply(
        lambda r: _pipe_merge([r.get("ec_numbers"), r.get("reaction_ec_numbers")]),
        axis=1,
    )
    combined_df["gene_name"] = combined_df.get("gene", pd.Series(dtype=str)).apply(
        _first_token
    )
    combined_df = combined_df.drop(columns=["ec_numbers"], errors="ignore")

    if "mapping_uniprot_id" in combined_df.columns:
        combined_df["mapping_uniprot_id"] = (
            combined_df["mapping_uniprot_id"].fillna("").astype(str)
        )

    resolved_ids = lookup_series.reindex(
        chembl_df.index, fill_value=UNIPROT_MISSING_VALUE
    )
    if merge_column in combined_df.columns:
        combined_df[merge_column] = combined_df[merge_column].fillna(
            UNIPROT_MISSING_VALUE
        )
    if merge_column == "uniprot_id":
        combined_df[merge_column] = resolved_ids.reindex(
            combined_df.index, fill_value=UNIPROT_MISSING_VALUE
        )
        if "mapping_uniprot_id" in combined_df.columns:
            empty_mask = (
                combined_df[merge_column].astype(str).str.strip()
                == UNIPROT_MISSING_VALUE
            )
            combined_df.loc[empty_mask, merge_column] = combined_df.loc[
                empty_mask, "mapping_uniprot_id"
            ]

    chembl_reactions = chembl_df.get(
        "reaction_ec_numbers", pd.Series(dtype=object)
    ).reindex(combined_df.index, fill_value=UNIPROT_MISSING_VALUE)
    uniprot_reaction_map = (
        uniprot_df.set_index("uniprot_id")
        if "uniprot_id" in uniprot_df.columns
        else pd.DataFrame()
    )
    if "reaction_ec_numbers" in uniprot_reaction_map.columns:
        uniprot_reactions = resolved_ids.reindex(
            combined_df.index, fill_value=UNIPROT_MISSING_VALUE
        ).map(uniprot_reaction_map["reaction_ec_numbers"])
        uniprot_reactions = uniprot_reactions.fillna(UNIPROT_MISSING_VALUE)
    else:
        uniprot_reactions = pd.Series(
            UNIPROT_MISSING_VALUE,
            index=combined_df.index,
            dtype=object,
        )

    combined_df["reaction_ec_numbers"] = [
        normalize_reaction_ec_numbers([u, c])
        for u, c in zip(uniprot_reactions, chembl_reactions)
    ]

    from tempfile import NamedTemporaryFile

    with NamedTemporaryFile(
        "w", delete=False, encoding=cfg.io.csv_encoding, newline=""
    ) as tmp:
        combined_df.to_csv(
            tmp, index=False, sep=cfg.io.csv_sep, encoding=cfg.io.csv_encoding
        )

        iuphar_input = Path(tmp.name)

    iuphar_args = argparse.Namespace(input_csv=iuphar_input, output_csv=output_csv)
    orig_target = cfg.target.iuphar.target_csv
    orig_family = cfg.target.iuphar.family_csv
    cfg.target.iuphar.target_csv = cfg.target.all.target_csv
    cfg.target.iuphar.family_csv = cfg.target.all.family_csv
    try:
        if run_iuphar(cfg, iuphar_args) != 0:
            raise RuntimeError("IUPHAR retrieval failed")
    finally:
        cfg.target.iuphar.target_csv = orig_target
        cfg.target.iuphar.family_csv = orig_family
        iuphar_input.unlink(missing_ok=True)

    try:
        iuphar_df = pd.read_csv(
            output_csv, sep=cfg.io.csv_sep, encoding=cfg.io.csv_encoding, dtype=str
        )
    except pd.errors.EmptyDataError:
        logger.warning("empty_iuphar_output", path=str(output_csv))
        iuphar_df = pd.DataFrame({"uniprot_id": pd.Series(dtype=object)})
    else:
        if "uniprot_id" not in iuphar_df.columns:
            logger.warning(
                "missing_iuphar_uniprot_column",
                path=str(output_csv),
                columns=list(iuphar_df.columns),
            )
            placeholder = pd.Series(
                UNIPROT_MISSING_VALUE, index=iuphar_df.index, dtype=object
            )
            iuphar_df = iuphar_df.copy()
            iuphar_df.insert(0, "uniprot_id", placeholder)
    existing_cols = set(combined_df.columns)
    classification_cols = [c for c in iuphar_df.columns if c not in existing_cols]
    iuphar_df = iuphar_df[["uniprot_id", *classification_cols]].copy()
    logger.info("fetch_iuphar_done", rows=len(iuphar_df))
    return combined_df, iuphar_df


def merge_results(
    combined_df: pd.DataFrame,
    iuphar_df: pd.DataFrame,
    cfg: Config,
    *,
    classifier: ii.IUPHARClassifier | None = None,
) -> pd.DataFrame:
    """Merge ChEMBL, UniProt and IUPHAR data into a single table.

    Parameters
    ----------
    combined_df : pandas.DataFrame
        DataFrame containing ChEMBL and UniProt information produced by
        :func:`fetch_iuphar`.
    iuphar_df : pandas.DataFrame
        IUPHAR classification results aligned with ``combined_df`` on
        ``uniprot_id``.
    cfg : Config
        Application configuration providing classifier settings.
    classifier : library.iuphar_library.IUPHARClassifier, optional
        Pre-initialised classifier. When ``None`` a classifier is created from
        ``cfg``.

    Returns
    -------
    pandas.DataFrame
        Merged and post-processed target data ready for validation and export.
    """

    logger.info("merge_results_start")
    merged = combined_df.merge(iuphar_df, on="uniprot_id", how="left")
    if classifier is None:
        classifier = pc.classifier_from_config(cfg)
    merged = pc.append_protein_class_predictions(merged, classifier)
    processed = tp.postprocess_targets(merged)
    final_df = tp.finalise_targets(processed)
    logger.info("merge_results_done", rows=len(final_df))
    return final_df


def validate_and_write(
    df: pd.DataFrame,
    output: Path,
    cfg: Config,
    *,
    raw_out: Path | None = None,
    id_cols: Sequence[str] | None = None,
    raw_format: str = "csv",
    reindex_raw: bool = True,
) -> int:
    """Normalise, validate and export the target table.

    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame produced by :func:`merge_results`.
    output : pathlib.Path
        Destination CSV path.
    cfg : Config
        Application configuration providing schema definitions and IO settings.

    Returns
    -------
    int
        ``0`` on success. Non-zero values indicate validation errors or
        failures when generating table-quality reports.
    """

    logger.info("validate_write_start", output=str(output))

    key_columns = list(id_cols) if id_cols else ["target_chembl_id"]
    normalized_output = _normalized_output_path(output)
    input_rows = len(df)

    if raw_out is not None:
        raw_format_value = (raw_format or "csv").lower()
        raw_frame = df.copy()
        raw_order: Sequence[str] | None = (
            TARGETS_COLUMN_ORDER if reindex_raw else None
        )
        if raw_order is not None:
            raw_frame = raw_frame.reindex(columns=raw_order)
        if raw_format_value == "parquet":
            try:
                raw_frame.to_parquet(raw_out, index=False)
            except ImportError as exc:  # pragma: no cover - optional dependency
                raise ValueError(
                    "Parquet export requires optional pyarrow or fastparquet"
                ) from exc
        else:
            io.write_csv(
                raw_frame,
                raw_out,
                cfg=cfg,
                sep=cfg.io.csv_sep,
                encoding=cfg.io.csv_encoding,
                key_cols=key_columns or None,
                col_order=raw_order,
            )


    normalized = normalize_targets(df)
    normalized_rows = len(normalized)
    logger.info(
        "normalize_targets_counts",
        before=input_rows,
        after=normalized_rows,
    )
    normalized = add_pipeline_metadata(normalized)
    prepared, missing_required, missing_optional = _prepare_targets_for_schema(
        normalized
    )
    final_df = prepared
    logger.info("prepare_targets_rows", rows=len(final_df))

    drop_reasons: dict[str, set[Any]] = {
        "schema": set(),
        "fk": set(),
        "regex": set(),
        "not_null": set(),
    }

    def _categorize_failure(row: pd.Series) -> str:
        check_value = str(row.get("check", "")).lower()
        context_value = str(row.get("schema_context", "")).lower()
        failure_value = str(row.get("failure_case", "")).lower()
        if "fk" in check_value or "foreign" in check_value or "fk" in context_value:
            return "fk"
        if "regex" in check_value or "match" in check_value or "pattern" in check_value:
            return "regex"
        if (
            "notnull" in check_value
            or "not null" in check_value
            or "nullable" in check_value
            or "null" in failure_value
        ):
            return "not_null"
        return "schema"

    exit_code = 0
    if not missing_required:
        if missing_optional:
            logger.warning(
                "optional_columns_missing",
                columns=sorted(missing_optional),
            )
        logger.info("targets_schema_validate_start", rows=len(final_df))
        try:
            final_df = TargetsSchema.validate(final_df, lazy=True)
        except SchemaErrors as exc:
            failure_path = output.with_name(f"{output.stem}_failure_cases.csv")
            errors = SidecarErrors()
            for row in exc.failure_cases.to_dict("records"):
                errors.add_error(row)
            errors.save(failure_path, cfg=cfg)
            logger.error(
                "validation_failed",
                failures=len(exc.failure_cases),
                path=str(failure_path),
            )
            failure_cases = exc.failure_cases.copy()
            for _, failure_row in failure_cases.iterrows():
                reason = _categorize_failure(failure_row)
                identifier = failure_row.get("index")
                if pd.isna(identifier):
                    identifier = (
                        failure_row.get("column"),
                        failure_row.get("failure_case"),
                    )
                drop_reasons[reason].add(identifier)
            logger.info(
                "targets_schema_validate_result",
                status="failed",
                rows=len(final_df),
                failures=len(failure_cases),
            )
            final_df = getattr(exc, "validated_data", final_df)
            exit_code = 1
        else:
            logger.info(
                "targets_schema_validate_result",
                status="success",
                rows=len(final_df),
                failures=0,
            )
    else:

        logger.warning(
            "validation_skipped",
            missing_columns=sorted(missing_required),
        )


    total_dropped = sum(len(ids) for ids in drop_reasons.values())
    logger.info(
        "validation_drop_stats",
        total=total_dropped,
        schema=len(drop_reasons["schema"]),
        fk=len(drop_reasons["fk"]),
        regex=len(drop_reasons["regex"]),
        not_null=len(drop_reasons["not_null"]),
    )

    placeholders_before_fill = int(final_df.isna().sum().sum())
    logger.info("placeholder_fillna_pending", replacements=placeholders_before_fill)
    final_df = final_df.fillna("-")
    before_dedup = len(final_df)
    final_df = final_df.drop_duplicates()
    logger.info("deduplicated_rows", dropped=before_dedup - len(final_df))
    io.write_csv(
        final_df,
        normalized_output,
        cfg=cfg,
        sep=cfg.io.csv_sep,
        encoding=cfg.io.csv_encoding,
        col_order=TARGETS_COLUMN_ORDER,
        key_cols=key_columns or None,
    )
    if final_df.empty:
        logger.info(
            "quality_report_skipped",
            reason="empty_dataframe",
            table=str(normalized_output.with_suffix("")),
        )
    else:
        try:
            analyze_table_quality(
                final_df, table_name=str(normalized_output.with_suffix(""))
            )
        except ValueError as exc:
            logger.error(
                "quality_report_failed",
                error=str(exc),
                path=str(normalized_output),
            )
            return 1
    logger.info("validate_write_done", rows=len(final_df))
    return exit_code


def run_all(cfg: Config, args: argparse.Namespace) -> int:
    """Run the full target acquisition pipeline.

    Parameters
    ----------
    cfg : Config
        Application configuration combining defaults for every individual
        pipeline.
    args : argparse.Namespace
        Parsed command-line arguments produced by :func:`build_parser`.

    Returns
    -------
    int
        ``0`` on success. Non-zero values indicate failures while fetching
        data, validating the merged dataset or writing the resulting CSV. The
        offending step is logged in the ``pipeline_step_failed`` event.
    """

    limit = cfg.target.all.limit
    if limit is not None and limit < 0:
        logger.error(
            "invalid_limit",
            section="target.all.limit",
            limit=limit,
        )
        return 1

    try:
        final_candidate = getattr(args, "final_out", None)
        if final_candidate in (None, argparse.SUPPRESS):
            final_output = Path(io.default_output_path(args.input_csv, cfg.io))
        else:
            final_output = Path(final_candidate)

        raw_candidate = getattr(args, "raw_out", None)
        if raw_candidate in (None, argparse.SUPPRESS):
            raw_output: Path | None = None
        else:
            raw_output = Path(raw_candidate)
        chembl_out = cfg.target.all.chembl_out or final_output.with_name(
            final_output.stem + "_chembl.csv"
        )
        uniprot_out = cfg.target.all.uniprot_out or final_output.with_name(
            final_output.stem + "_uniprot.csv"
        )
        iuphar_out = cfg.target.all.iuphar_out or final_output.with_name(
            final_output.stem + "_iuphar.csv"
        )

        raw_format = str(getattr(args, "raw_format", "csv") or "csv").lower()
        reindex_raw = not bool(getattr(args, "no_reindex_raw", False))
        id_cols_value = getattr(args, "id_cols", None)
        if id_cols_value in (None, argparse.SUPPRESS):
            key_columns = ["target_chembl_id"]
        elif isinstance(id_cols_value, str):
            key_columns = [id_cols_value]
        else:
            key_columns = list(id_cols_value) or ["target_chembl_id"]

        chembl_df = fetch_chembl(
            cfg,
            args.input_csv,
            chembl_out,
            limit=limit,
            chunk_size=cfg.target.all.chunk_size,
            offset=getattr(args, "offset", 0),
            id_cols=key_columns,
        )
        uniprot_df = fetch_uniprot(cfg, chembl_df, uniprot_out)
        combined_df, iuphar_df = fetch_iuphar(cfg, chembl_df, uniprot_df, iuphar_out)
        merged = merge_results(combined_df, iuphar_df, cfg)
        exit_code = validate_and_write(
            merged,
            final_output,
            cfg,
            raw_out=raw_output,
            id_cols=key_columns,
            raw_format=raw_format,
            reindex_raw=reindex_raw,
        )
        return exit_code
    except (FileNotFoundError, ValueError, OSError, RuntimeError) as exc:
        logger.error(
            "pipeline_step_failed",
            error=str(exc),
            step="all",
            input=str(args.input_csv),
            output=str(final_output),
        )
        return 1


def run(cfg: Config, args: argparse.Namespace) -> int:
    """Execute the selected target pipeline applying CLI policies."""

    final_candidate = getattr(args, "final_out", None)
    if final_candidate in (None, argparse.SUPPRESS):
        final_output = Path(io.default_output_path(args.input_csv, cfg.io))
    else:
        final_output = Path(final_candidate)
    args.final_out = final_output
    args.output_csv = final_output
    if args.skip_existing and final_output.exists() and not args.force:
        logger.info("pipeline_skip_existing", output=str(final_output))
        return 0
    if getattr(args, "_out_alias_used", False):
        logger.warning(
            "deprecated_output_alias_used",
            alias="--out",
            replacement="--output",
        )
    func = getattr(args, "func", None)
    if func is None:
        logger.error(
            "missing_subcommand_handler", command=getattr(args, "command", "")
        )
        return 1
    result = func(cfg, args)
    return int(result)


def main(argv: Sequence[str] | None = None) -> int:
    """Command line entry point using :class:`Config` for defaults.

    Parameters
    ----------
    argv : Sequence[str] | None, optional
        Command-line arguments to parse. When ``None`` the values from
        :data:`sys.argv` are used.

    Returns
    -------
    int
        ``0`` when the selected target pipeline succeeds, non-zero otherwise.

    Raises
    ------
    SystemExit
        Raised when invalid command-line options are provided to the parser.
    """
    parser, log_cfg = build_parser()
    args = parser.parse_args(argv)
    prepare_io_paths(
        args,
        input_default=DEFAULT_INPUT_NAME,
        output_stem=DEFAULT_OUTPUT_STEM,
    )
    date_value = getattr(args, "date", None)
    if not isinstance(date_value, str) or not date_value:
        date_value = datetime.now(timezone.utc).strftime("%Y%m%d")
        setattr(args, "date", date_value)
    log_dir = Path("logs")
    log_dir.mkdir(parents=True, exist_ok=True)
    log_path = log_dir / f"get_target_data_{date_value}.log"
    with log_path.open("a", encoding="utf-8") as log_stream:
        log_cfg.stream = log_stream
        configure_logger(log_cfg)

        limit_value = getattr(args, "limit", None)
        if limit_value == 0:
            logger.info("pipeline_skip_limit", limit=limit_value)
            return 0
        subparser_map = getattr(parser, "subparsers_map", {})
        subparser = subparser_map.get(args.command, parser)
        if limit_value is not None and limit_value < 0:
            subparser.error("--limit must be zero or a positive integer")
        offset_value = getattr(args, "offset", 0)
        if offset_value < 0:
            subparser.error("--offset must be zero or a positive integer")
        mapping: dict[str, str] = {}
        if args.command == "uniprot":
            mapping = {
                "column": "target.uniprot.column",
                "data_dir": "target.uniprot.data_dir",
                "limit": "target.uniprot.limit",
            }
        elif args.command == "chembl":
            mapping = {
                "column": "target.chembl.column",
                "chunk_size": "target.chembl.chunk_size",
                "timeout": "target.chembl.timeout",
                "limit": "target.chembl.limit",
            }
        elif args.command == "iuphar":
            mapping = {
                "target_csv": "target.iuphar.target_csv",
                "family_csv": "target.iuphar.family_csv",
                "limit": "target.iuphar.limit",
            }
        elif args.command == "all":
            mapping = {
                "timeout": "target.all.timeout",
                "data_dir": "target.all.data_dir",
                "target_csv": "target.all.target_csv",
                "family_csv": "target.all.family_csv",
                "uniprot_column": "target.all.uniprot_column",
                "chembl_out": "target.all.chembl_out",
                "uniprot_out": "target.all.uniprot_out",
                "iuphar_out": "target.all.iuphar_out",
                "limit": "target.all.limit",
            }
        args_dict = vars(args).copy()
        output_candidate = args_dict.get("output_csv")
        if output_candidate in (
            None,
            argparse.SUPPRESS,
        ):
            final_value = args_dict.get("final_out")
            if final_value in (None, argparse.SUPPRESS):
                date_token = args_dict.get(
                    "date", datetime.now(timezone.utc).strftime("%Y%m%d")
                )
                inferred = Path(args_dict["input_csv"]).with_name(
                    f"output.{DEFAULT_OUTPUT_STEM}_{date_token}.csv"
                )
                args_dict["final_out"] = inferred
                args_dict["output_csv"] = inferred
            else:
                candidate = Path(final_value)
                args_dict["final_out"] = candidate
                args_dict["output_csv"] = candidate
        else:
            args_dict["output_csv"] = Path(output_candidate)
        exit_code = run_cli_command(
            args=argparse.Namespace(**args_dict),
            parser=subparser,
            base_parser=parser,
            log_cfg=log_cfg,
            mapping=mapping,
            run=run,
            logger=logger,
        )

    return exit_code


if __name__ == "__main__":  # pragma: no cover - CLI entry point
    raise SystemExit(main())
