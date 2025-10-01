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
from collections.abc import Iterator, Sequence
from functools import partial
from itertools import islice
from pathlib import Path
from typing import cast

from library.utils.bootstrap import ensure_project_root


if __package__ in {None, ""}:
    ensure_project_root()

import pandas as pd
import requests
from pandera.errors import SchemaErrors

from library import chembl_library as cl
from library import cli
from library import io
from library import iuphar_library as ii
from library import protein_classification as pc
from library import target_postprocessing as tp
from library import uniprot_library as uu
from library.clients import ChemblClient, _chunked
from library.cli_utils import PipelineError, run_pipeline
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
    ensure_dirs,
    print_config,
)
from library.log import logger
from library.metadata import Stats, file_sha256, write_meta_yaml
from library.pipeline_metadata import add_pipeline_metadata
from library.sidecar import SidecarErrors
from library.table_quality import analyze_table_quality
from library.validation import ValidationResult
from schemas import TargetsSchema, normalize_targets
from schemas.targets import TARGETS_COLUMN_ORDER

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
    """Write ``df`` to a uniquely named snapshot CSV file.

    The file is created alongside ``base`` using the pattern
    ``<base>_<step>_<n>.csv`` where ``n`` increments to avoid overwriting
    existing files.

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
        Application configuration (currently unused but kept for API
        compatibility).

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
            df.to_csv(candidate, index=False)
            return candidate
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
    )
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

        output = args.output_csv or io.default_output_path(args.input_csv, cfg.io)
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
    processed_ids = 0

    def _iter_ids() -> Iterator[str]:
        nonlocal processed_ids
        for identifier in ids_iter:
            processed_ids += 1
            yield identifier

    limited_ids: Iterator[str]
    if limit is not None:
        limited_ids = islice(_iter_ids(), limit)
    else:
        limited_ids = _iter_ids()

    output_path = args.output_csv or io.default_output_path(args.input_csv, cfg.io)
    output = Path(output_path)
    failure_path = output.with_name(f"{output.stem}_failure_cases.csv")

    missing_optional_columns: set[str] = set()

    def _prepare_chunk(frame: pd.DataFrame) -> pd.DataFrame:
        prepared, _, missing_optional = _prepare_targets_for_schema(frame)
        missing_optional_columns.update(missing_optional)
        return prepared

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
        with ChemblClient(cfg.api, cfg.retry, cfg.chembl) as client:
            for chunk_ids in _chunked(limited_ids, cfg.target.chembl.chunk_size):
                if not chunk_ids:
                    continue
                try:
                    chunk_df = cl.get_targets(
                        chunk_ids,
                        cfg=cfg.api,
                        client=client,
                        mapping_cfg=cfg.uniprot_mapping,
                        chunk_size=cfg.target.chembl.chunk_size,
                        timeout=cfg.target.chembl.timeout,
                    )
                except (requests.RequestException, ValueError) as exc:
                    logger.error(
                        "chembl_fetch_failed",
                        error=str(exc),
                        chunk_size=cfg.target.chembl.chunk_size,
                        timeout=cfg.target.chembl.timeout,
                    )
                    raise PipelineError(str(exc)) from exc
                yield chunk_df

    def writer(
        chunks: Iterator[pd.DataFrame],
        destination: Path,
        col_order: Sequence[str],
        key_cols: Sequence[str],
    ) -> Path:
        key_columns = list(key_cols)
        return io.write_csv(
            chunks,
            destination,
            cfg=cfg,
            sep=cfg.io.csv_sep,
            encoding=cfg.io.csv_encoding,
            key_cols=key_columns or None,
            col_order=col_order,
            chunksize=cfg.io.csv_chunksize,
        )

    table_quality = partial(
        analyze_table_quality,
        table_name=str(output.with_suffix("")),
    )

    metadata_hooks = [
        normalize_targets,
        add_pipeline_metadata,
        _prepare_chunk,
    ]

    exit_code = run_pipeline(
        fetcher=fetcher,
        schema=TargetsSchema,
        schema_name="TargetsSchema",
        validators=[_validate_chunk],
        metadata_hooks=metadata_hooks,
        writer=writer,
        output_path=output,
        failure_path=failure_path,
        command=" ".join(sys.argv),
        config_snapshot=_serialize_paths(cfg.to_dict()),
        inputs={"input_csv": str(args.input_csv)},
        key_columns=["target_chembl_id"],
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
        output = args.output_csv or io.default_output_path(args.input_csv, cfg.io)
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
    output_csv: Path,
    limit: int | None = None,
    *,
    chunk_size: int | None = None,
    offset: int = 0,
) -> pd.DataFrame:
    """Fetch target information from ChEMBL.

    Parameters
    ----------
    cfg : Config
        Application configuration used to drive the ChEMBL pipeline.
    input_csv : pathlib.Path
        Source CSV containing target identifiers.
    output_csv : pathlib.Path
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
        Retrieved ChEMBL data loaded from ``output_csv``.

    Raises
    ------
    RuntimeError
        Raised when :func:`run_chembl` reports a non-zero exit code.
    """

    logger.info("fetch_chembl_start", input=str(input_csv), output=str(output_csv))
    chembl_args = argparse.Namespace(
        input_csv=input_csv, output_csv=output_csv, offset=offset
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
    df = pd.read_csv(
        output_csv, sep=cfg.io.csv_sep, encoding=cfg.io.csv_encoding, dtype=str
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
    column = cfg.target.all.uniprot_column
    series = chembl_df.get(column)
    if series is None:
        uids: list[str] = []
    else:
        values = series.to_numpy(copy=False)
        uids = [uid for uid in values if isinstance(uid, str) and uid]
    from tempfile import NamedTemporaryFile

    with NamedTemporaryFile(
        "w", delete=False, encoding=cfg.io.csv_encoding, newline=""
    ) as tmp:
        tmp_path = Path(tmp.name)

    pd.DataFrame({"uniprot_id": uids}).to_csv(
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

    df = pd.read_csv(
        output_csv, sep=cfg.io.csv_sep, encoding=cfg.io.csv_encoding, dtype=str
    )
    df["original_id"] = pd.Series(uids, index=df.index, dtype=object)
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

    combined_df = pd.merge(
        chembl_for_merge,
        uniprot_df,
        left_on=merge_column,
        right_on="original_id",
        how="left",
    )

    if "original_id" in combined_df.columns:
        combined_df = combined_df.drop(columns=["original_id"])
    if merge_column == "uniprot_id":
        left_series = combined_df.pop("uniprot_id_x") if "uniprot_id_x" in combined_df else None
        right_series = combined_df.pop("uniprot_id_y") if "uniprot_id_y" in combined_df else None
        if left_series is not None and right_series is not None:
            combined_df["uniprot_id"] = right_series.fillna(left_series)
        elif left_series is not None:
            combined_df["uniprot_id"] = left_series
        elif right_series is not None:
            combined_df["uniprot_id"] = right_series

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


def validate_and_write(df: pd.DataFrame, output: Path, cfg: Config) -> int:
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
    normalized = normalize_targets(df)
    normalized = add_pipeline_metadata(normalized)
    final_df, missing_required, missing_optional = _prepare_targets_for_schema(
        normalized
    )
    exit_code = 0
    if not missing_required:
        if missing_optional:
            logger.warning(
                "optional_columns_missing",
                columns=sorted(missing_optional),
            )
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
            final_df = getattr(exc, "validated_data", final_df)
            exit_code = 1
    else:
        logger.warning(
            "validation_skipped",
            missing_columns=sorted(missing_required),
        )
    final_df = final_df.fillna("-")
    final_df = final_df.drop_duplicates()
    io.write_csv(
        final_df,
        output,
        cfg=cfg,
        sep=cfg.io.csv_sep,
        encoding=cfg.io.csv_encoding,
        col_order=TARGETS_COLUMN_ORDER,
    )
    if final_df.empty:
        logger.info(
            "quality_report_skipped",
            reason="empty_dataframe",
            table=str(output.with_suffix("")),
        )
    else:
        try:
            analyze_table_quality(final_df, table_name=str(output.with_suffix("")))
        except ValueError as exc:
            logger.error(
                "quality_report_failed",
                error=str(exc),
                path=str(output),
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
        output = args.output_csv or io.default_output_path(args.input_csv, cfg.io)
        chembl_out = cfg.target.all.chembl_out or output.with_name(
            output.stem + "_chembl.csv"
        )
        uniprot_out = cfg.target.all.uniprot_out or output.with_name(
            output.stem + "_uniprot.csv"
        )
        iuphar_out = cfg.target.all.iuphar_out or output.with_name(
            output.stem + "_iuphar.csv"
        )

        chembl_df = fetch_chembl(
            cfg,
            args.input_csv,
            chembl_out,
            limit=limit,
            chunk_size=cfg.target.all.chunk_size,
            offset=getattr(args, "offset", 0),
        )
        uniprot_df = fetch_uniprot(cfg, chembl_df, uniprot_out)
        combined_df, iuphar_df = fetch_iuphar(cfg, chembl_df, uniprot_df, iuphar_out)
        merged = merge_results(combined_df, iuphar_df, cfg)
        exit_code = validate_and_write(merged, output, cfg)
        return exit_code
    except (FileNotFoundError, ValueError, OSError, RuntimeError) as exc:
        logger.error(
            "pipeline_step_failed",
            error=str(exc),
            step="all",
            input=str(args.input_csv),
            output=str(args.output_csv or output),
        )
        return 1


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
    subparser_map = getattr(parser, "subparsers_map", {})
    subparser = subparser_map.get(args.command, parser)
    limit_value = getattr(args, "limit", None)
    if limit_value is not None and limit_value <= 0:
        subparser.error("--limit must be a positive integer")
    offset_value = getattr(args, "offset", 0)
    if offset_value < 0:
        subparser.error("--offset must be zero or a positive integer")
    log_cfg.level = args.log_level
    logger = configure_logger(log_cfg)
    logger.info("pipeline_start", run_id=log_cfg.run_id)
    try:
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
        cfg: Config = cli.apply_config_overrides(
            args, subparser, args.config, mapping=mapping, base_parser=parser
        )
        if args.print_config:
            print_config(cfg)
            configure_logger(log_cfg)
            logger.info("pipeline_done", run_id=log_cfg.run_id)
            return 0
        ensure_dirs(cfg)
        output_path = Path(
            args.output_csv or io.default_output_path(args.input_csv, cfg.io)
        )
        args.output_csv = output_path
        if args.skip_existing and output_path.exists() and not args.force:
            logger.info("pipeline_skip_existing", output=str(output_path))
            logger.info("pipeline_done", run_id=log_cfg.run_id)
            return 0
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
        )
        logger.info("pipeline_fail", run_id=log_cfg.run_id)
        return 1
    if hasattr(args, "func"):
        exit_code = cast(int, args.func(cfg, args))
        if exit_code == 0:
            logger.info("pipeline_done", run_id=log_cfg.run_id)
        else:
            logger.info("pipeline_fail", run_id=log_cfg.run_id)
        return exit_code
    parser.print_help()
    logger.info("pipeline_fail", run_id=log_cfg.run_id)
    return 1


if __name__ == "__main__":  # pragma: no cover - CLI entry point
    raise SystemExit(main())
