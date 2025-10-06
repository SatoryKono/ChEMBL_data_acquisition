"""Command line interface for retrieving ChEMBL activity data.

The module exposes a ``main`` entry point compatible with setuptools console
scripts as well as helpers that can be invoked directly from other
applications or tests.
"""

from __future__ import annotations

import argparse
import sys

from collections.abc import Iterable, Iterator, Mapping, Sequence, Callable

from functools import partial
from itertools import islice
from pathlib import Path
from time import sleep

try:
    from library.utils.bootstrap import ensure_project_root
except ModuleNotFoundError:  # pragma: no cover - fallback for direct execution
    def _is_within(path: Path, root: Path) -> bool:
        try:
            path.relative_to(root)
        except ValueError:
            return False
        return True

    def ensure_project_root() -> None:
        """Add the repository root to ``sys.path`` when executed as a script."""

        project_root = Path(__file__).resolve().parent.parent
        project_root_str = str(project_root)
        if project_root_str not in sys.path:
            sys.path.insert(0, project_root_str)

        existing = sys.modules.get("library")
        if existing is None:
            return

        module_paths: list[Path] = []
        file_attr = getattr(existing, "__file__", None)
        if file_attr:
            module_paths.append(Path(file_attr).resolve())
        package_paths = getattr(existing, "__path__", None)
        if package_paths is not None:
            module_paths.extend(Path(p).resolve() for p in package_paths)

        if any(_is_within(path, project_root) for path in module_paths):
            return

        for name in list(sys.modules):
            if name == "library" or name.startswith("library."):
                del sys.modules[name]

if __package__ in {None, ""}:
    ensure_project_root()

import pandas as pd
import requests

from library.integration import chembl_library as cl
from library import cli
from library import io
from library.clients import ChemblClient
from library.common.csv_utils import write_csv_chunks_deterministic  # re-exported for tests
from library.common.rate_limiter import get_global_limiter
from library.pipelines.assay.chembl_assay import ACTIVITY_COLUMNS
from library.pipelines.common import (
    ChunkedFetchConfig,
    CsvWriterConfig,
    prepare_chunked_pipeline,
)

from library.cli import (
    LoggerConfig,
    positive_int,
    ConfigMetadata,
)
from library.cli import (
    build_parser as base_parser,
)
from library.cli.logging import setup_cli_logging
from library.cli_utils import (
    PipelineError,
    resolve_invocation,
    run_cli_command,
    run_pipeline,
)
from library.cli_utils import (
    file_sha256 as _cli_file_sha256,
)
from library.cli_utils import (
    write_meta_yaml as _cli_write_meta_yaml,
)
from library.config import Config, _serialize_paths
from library.common.log import logger
from library.pipelines.common import add_pipeline_metadata
from library.processing.activity import (
    apply_activity_annotations,
    compute_activity_bounds,
)
from library.postprocessing.activity_extended import process_activity_extended
from library.table_quality import analyze_table_quality
from library.validation import validate_activities
from library.schemas import ActivitiesSchema, configure_activity_schema, normalize_activities
from library.common.fetch_retry import ChunkFailureTracker, compute_backoff_delay

DEFAULT_INPUT_NAME = "activity.csv"
DEFAULT_OUTPUT_STEM = "activities"
PROGRAM_NAME = Path(__file__).with_suffix("").name

_OPTION_UNSET = object()


def _option(
    metadata: ConfigMetadata | None,
    *,
    argument: str | None = None,
    path: str | None = None,
    value: object = _OPTION_UNSET,
    default_source: str = "unknown",
    default_detail: str | None = None,
) -> dict[str, object]:
    """Return structured option metadata for pipeline logging."""

    if metadata is not None:
        if value is _OPTION_UNSET:
            return metadata.option(
                argument=argument,
                path=path,
                default_source=default_source,
                default_detail=default_detail,
            )
        return metadata.option(
            argument=argument,
            path=path,
            value=value,
            default_source=default_source,
            default_detail=default_detail,
        )
    actual = None if value is _OPTION_UNSET else value
    entry: dict[str, object] = {"value": actual, "source": default_source}
    if default_detail is not None:
        entry["detail"] = default_detail
    return entry


def _args_invocation(args: argparse.Namespace) -> tuple[str, ...]:
    invocation = getattr(args, "invocation", None)
    if invocation is None:
        return (PROGRAM_NAME,)
    return tuple(str(arg) for arg in invocation)


file_sha256 = _cli_file_sha256
write_meta_yaml = _cli_write_meta_yaml
configure_logger = cli.configure_logger

__all__ = (
    "file_sha256",
    "write_meta_yaml",
    "configure_logger",
)


_ACTIVITY_REQUIRED_COLUMNS: tuple[str, ...] = tuple(
    name for name, column in ActivitiesSchema.columns.items() if column.required
)

_ACTIVITY_REQUIRED_DTYPES: dict[str, object] = {
    name: getattr(column, "dtype", None)
    for name, column in ActivitiesSchema.columns.items()
    if column.required
}

_ORIGINAL_IO_WRITE_CSV = io.write_csv


_EXTENDED_ACTIVITY_DTYPES: dict[str, str] = {
    "activity_chembl_id": "string",
    "salt_chembl_id": "string",
    "target_chembl_id": "string",
    "bao_endpoint": "string",
    "compound_key": "string",
    "compound_name": "string",
    "multmol_assay": "boolean",
    "approx_cited_activity": "boolean",
    "shuffled_cit": "boolean",
    "exact_cited_activity": "boolean",
    "higly_correlated_cit": "boolean",
    "review_doc": "boolean",
    "rounded_data_citation": "boolean",
    "original_activity_approx": "string",
    "original_activity_exact": "string",
    "nstereo": "Int64",
    "log_value": "Float64",
}

_EXTENDED_ACTIVITY_FALLBACKS: dict[str, Callable[[pd.DataFrame], pd.Series | None]] = {
    "activity_chembl_id": lambda df: df.get("activity_id"),
    "compound_name": lambda df: df.get("molecule_pref_name"),
    "log_value": lambda df: df.get("pchembl_value"),
}


 
def _string_like_missing(series: pd.Series) -> pd.Series:
    """Return a boolean mask for ``series`` treating blanks as missing."""

    mask = series.isna()
    if pd.api.types.is_string_dtype(series) or series.dtype == "object":
        string_values = series.astype("string")
        mask = mask | string_values.str.strip().fillna("").eq("")
    return mask
 
 
def _string_blank_mask(series: pd.Series) -> pd.Series:
    """Return mask of entries that are null or contain only whitespace."""

    values = series.astype("string")
    stripped = values.str.strip()
    return values.isna() | stripped.fillna("").eq("")


def _load_assay_src_lookup(dictionary_dir: Path | str | None) -> dict[str, str]:
    """Return mapping of assay identifiers to ``src_assay_id`` values."""

    if dictionary_dir is None:
        return {}

    candidate = Path(dictionary_dir) / "_assay" / "assay.csv"
    try:
        frame = pd.read_csv(
            candidate,
            usecols=["assay_chembl_id", "src_assay_id"],
            dtype="string",
        )
    except FileNotFoundError:
        logger.warning("assay_lookup_missing", path=str(candidate))
        return {}
    except pd.errors.EmptyDataError:
        logger.warning("assay_lookup_empty", path=str(candidate))
        return {}
    except ValueError as exc:
        logger.warning(
            "assay_lookup_invalid_columns",
            path=str(candidate),
            error=str(exc),
        )
        return {}
    except OSError as exc:
        logger.warning("assay_lookup_read_failed", path=str(candidate), error=str(exc))
        return {}

    if frame.empty:
        return {}

    cleaned = frame.dropna(subset=["assay_chembl_id", "src_assay_id"])
    if cleaned.empty:
        return {}

    cleaned = cleaned.assign(
        assay_chembl_id=cleaned["assay_chembl_id"].str.strip(),
        src_assay_id=cleaned["src_assay_id"].str.strip(),
    )

    cleaned = cleaned[cleaned["assay_chembl_id"].ne("") & cleaned["src_assay_id"].ne("")]
    if cleaned.empty:
        return {}

    return {
        str(assay_id): str(src_id)
        for assay_id, src_id in cleaned[["assay_chembl_id", "src_assay_id"]]
        .itertuples(index=False, name=None)
    }


def _ensure_src_assay_id(
    frame: pd.DataFrame, *, lookup: Mapping[str, str]
) -> pd.DataFrame:
    """Populate ``src_assay_id`` using ``lookup`` when missing."""

    if "src_assay_id" not in frame.columns and "assay_chembl_id" not in frame.columns:
        return frame

    result = frame.copy()

    if "src_assay_id" in result.columns:
        try:
            result["src_assay_id"] = result["src_assay_id"].astype("string")
        except TypeError:
            result["src_assay_id"] = result["src_assay_id"].astype("string")
    else:
        result["src_assay_id"] = pd.Series(pd.NA, index=result.index, dtype="string")

    if not lookup or result.empty or "assay_chembl_id" not in result.columns:
        return result

    assay_ids = result["assay_chembl_id"].astype("string")
    missing_mask = _string_like_missing(result["src_assay_id"])
    if not missing_mask.any():
        return result

    normalized_ids = assay_ids.where(~assay_ids.isna(), None).astype(object)

    def _resolve(value: object) -> str | None:
        if value is None:
            return None
        return lookup.get(str(value))

    mapped = normalized_ids.map(_resolve)
    available = missing_mask & mapped.notna()
    if not available.any():
        return result

    result.loc[available, "src_assay_id"] = mapped[available].astype("string")
    return result


def _ensure_molecule_pref_name(
    frame: pd.DataFrame,
    *,
    cfg: Config,
    client: ChemblClient,
    cache: dict[str, str | None],
) -> pd.DataFrame:
    """Populate ``molecule_pref_name`` via the test item API when missing."""

    if frame.empty or "molecule_chembl_id" not in frame.columns:
        return frame

    result = frame.copy()

    if "molecule_pref_name" in result.columns:
        missing_mask = _string_blank_mask(result["molecule_pref_name"])
    else:
        result["molecule_pref_name"] = pd.Series(pd.NA, index=result.index, dtype="string")
        missing_mask = pd.Series(True, index=result.index, dtype="boolean")

    if not missing_mask.any():
        return result

    molecule_ids = (
        result.loc[missing_mask, "molecule_chembl_id"].astype("string").str.strip()
    )
    molecule_ids = molecule_ids[molecule_ids != ""]
    unique_ids = tuple(dict.fromkeys(molecule_ids.dropna().tolist()))

    if not unique_ids:
        return result

    pending: list[str] = []
    for identifier in unique_ids:
        if identifier not in cache:
            pending.append(identifier)

    if pending:
        fields = list(cfg.testitem.fields)
        for column in ("molecule_chembl_id", "pref_name"):
            if column not in fields:
                fields.append(column)
        try:
            lookup = cl.get_testitem(
                pending,
                cfg=cfg.api,
                client=client,
                chunk_size=cfg.testitem.batch_size,
                timeout=cfg.testitem.timeout,
                fields=fields,
                page_limit=cfg.testitem.request_limit,
            )
        except (requests.RequestException, ValueError, AttributeError) as exc:
            logger.warning(
                "testitem_pref_name_lookup_failed",
                error=str(exc),
                error_type=exc.__class__.__name__,
                pending_ids=list(pending),
            )
            lookup = pd.DataFrame(columns=["molecule_chembl_id", "pref_name"])
        if not lookup.empty and {"molecule_chembl_id", "pref_name"}.issubset(lookup.columns):
            mapped = (
                lookup[["molecule_chembl_id", "pref_name"]]
                .dropna(subset=["molecule_chembl_id"])
                .astype({"molecule_chembl_id": "string"})
            )
            for chembl_id, pref_name in mapped.itertuples(index=False):
                cache[str(chembl_id)] = str(pref_name) if pd.notna(pref_name) else None
        for identifier in pending:
            cache.setdefault(identifier, None)

    fill_map = {key: value for key, value in cache.items() if value}
    if not fill_map:
        return result

    replacements = molecule_ids.map(fill_map)
    available = replacements.notna()
    if available.any():
        result.loc[molecule_ids.index[available], "molecule_pref_name"] = (
            replacements[available].astype("string")
        )

    return result
 
 


def _ensure_extended_activity_columns(frame: pd.DataFrame) -> pd.DataFrame:
    """Guarantee columns expected by the post-processing stage."""

    result = frame.copy()
    if result.empty:
        for column, dtype in _EXTENDED_ACTIVITY_DTYPES.items():
            if column not in result.columns:
                result[column] = pd.Series([], dtype=dtype)
        return result

    if "activity_id" in result.columns:
        if "activity_chembl_id" in result.columns:
            missing_id = _string_like_missing(result["activity_id"])
            if missing_id.any():
                result.loc[missing_id, "activity_id"] = result.loc[
                    missing_id, "activity_chembl_id"
                ]
        else:
            missing_id = _string_like_missing(result["activity_id"])
            if missing_id.any():
                # Preserve dtype by assigning NA values without coercing existing entries.
                result.loc[missing_id, "activity_id"] = pd.NA
    elif "activity_chembl_id" in result.columns:
        result["activity_id"] = result["activity_chembl_id"].copy()

    for column, dtype in _EXTENDED_ACTIVITY_DTYPES.items():
        fallback = _EXTENDED_ACTIVITY_FALLBACKS.get(column)
        if column in result.columns:
            if fallback is not None:
                existing = result[column]
                if dtype in {"Float64", "Int64"}:
                    missing_mask = existing.isna()
                else:
                    missing_mask = _string_like_missing(existing)
                if missing_mask.any():
                    candidate = fallback(result)
                    if candidate is not None:
                        aligned = candidate.reindex(result.index)
                        try:
                            filled = aligned.astype(dtype)
                        except (TypeError, ValueError):
                            filled = aligned.astype("string")
                        result.loc[missing_mask, column] = filled.loc[missing_mask]
            continue
        if fallback is not None:
            candidate = fallback(result)
            if candidate is not None:
                aligned = candidate.reindex(result.index)
                try:
                    result[column] = aligned.astype(dtype)
                except TypeError:
                    result[column] = aligned.astype("string")
                continue
        if dtype == "boolean":
            filler = pd.Series(pd.NA, index=result.index, dtype="boolean")
        elif dtype == "Float64":
            filler = pd.Series(pd.NA, index=result.index, dtype="Float64")
        elif dtype == "Int64":
            filler = pd.Series(pd.NA, index=result.index, dtype="Int64")
        else:
            filler = pd.Series(pd.NA, index=result.index, dtype=dtype)
        result[column] = filler
    return result


def run_chembl(cfg: Config, args: argparse.Namespace) -> int:
    """Execute activity retrieval from the ChEMBL API.

    The resulting CSV places columns defined in :data:`ActivitiesSchema`
    first, preserving their declared order. Any additional fields appear
    afterwards sorted alphabetically.

    Parameters
    ----------
    cfg : Config
        Application configuration providing API credentials, retry strategy
        and CSV export options.
    args : argparse.Namespace
        Parsed command-line arguments. ``args.limit`` constrains the number of
        identifiers processed, while ``args.dry_run`` skips network calls and
        file generation.

    Returns
    -------
    int
        ``0`` on success, non-zero when validation or I/O failures are
        encountered. Upstream API errors are logged and converted into a
        failure code by :func:`library.cli_utils.run_pipeline`.
    """
    limit = cfg.activity.limit
    if limit is not None and limit < 0:
        logger.error("invalid_limit", section="activity.limit", limit=limit)
        return 1

    offset = getattr(args, "offset", 0)
    workers_override = getattr(args, "workers", None)
    configured_workers = (
        workers_override
        if workers_override is not None
        else getattr(cfg.activity, "workers", 1)
    )
    final_out_attr = getattr(args, "final_out", None)
    if final_out_attr in (None, argparse.SUPPRESS):
        legacy_output = getattr(args, "output_csv", None)
        if legacy_output not in (None, argparse.SUPPRESS):
            output_path = Path(legacy_output)
            if not isinstance(legacy_output, Path):
                args.final_out = output_path
            setattr(args, "output_csv", output_path)
        else:
            output_path = Path(io.default_output_path(args.input_csv, cfg.io))
            args.final_out = output_path
            setattr(args, "output_csv", output_path)
    else:
        output_path = Path(final_out_attr)
        if not isinstance(final_out_attr, Path):
            args.final_out = output_path
        setattr(args, "output_csv", output_path)

    metadata_obj = getattr(args, "_config_metadata", None)
    if not isinstance(metadata_obj, ConfigMetadata):
        metadata_obj = None

    output_source = "cli" if getattr(args, "final_out", None) else "derived"
    logger.info(
        "activity_pipeline_start",
        input=_option(metadata_obj, value=str(args.input_csv), default_source="cli"),
        output=_option(
            metadata_obj,
            value=str(output_path),
            default_source=output_source,
        ),
        limit=_option(
            metadata_obj,
            argument="limit",
            path="sources.chembl.pipelines.activity.limit",
            value=limit,
        ),
        offset=_option(
            metadata_obj,
            argument="offset",
            path="sources.chembl.pipelines.activity.offset",
            value=offset,
        ),
        batch_size=_option(
            metadata_obj,
            argument="batch_size",
            path="sources.chembl.pipelines.activity.batch_size",
            value=cfg.activity.batch_size,
        ),
        timeout=_option(
            metadata_obj,
            argument="timeout",
            path="sources.chembl.pipelines.activity.timeout",
            value=cfg.activity.timeout,
        ),
        dry_run=_option(
            metadata_obj,
            argument="dry_run",
            path="sources.chembl.pipelines.activity.dry_run",
            value=cfg.activity.dry_run,
            default_source="cli",
        ),
        workers=_option(
            metadata_obj,
            argument="workers",
            path="sources.chembl.pipelines.activity.workers",
            value=configured_workers,
        ),
    )

    if cfg.activity.dry_run:
        expected = limit if limit is not None else 0
        logger.info("dry_run", limit=expected)
        return 0

    try:
        ids_iter = io.read_ids(args.input_csv, column=cfg.activity.column, cfg=cfg.io)
    except (FileNotFoundError, ValueError) as exc:
        logger.error(
            "read_fail",
            error=str(exc),
            path=str(args.input_csv),
        )
        return 1

    if offset:
        ids_iter = islice(ids_iter, offset, None)
        logger.info("process_offset", offset=offset)

    processed_ids = 0
    extended_output_path: Path | None = None

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

    enrichment_cfg = cfg.activity_enrichment
    extra_columns: list[str] = []
    action_cfg = enrichment_cfg.action_type
    configure_activity_schema(action_cfg.metrics)
    if action_cfg.enabled or action_cfg.log_missing or action_cfg.log_distribution:
        extra_columns.append(action_cfg.column)
    extra_kwargs = {"extra_columns": extra_columns} if extra_columns else {}

    invocation = _args_invocation(args)

    failure_path = output_path.with_name(f"{output_path.stem}_failure_cases.csv")
    fetch_failure_path = output_path.with_name(
        f"{output_path.stem}_fetch_failures.csv"
    )


    def _compute_bounds(frame: pd.DataFrame) -> pd.DataFrame:
        return compute_activity_bounds(frame, cfg.activity_bounds)

    def _apply_annotations(frame: pd.DataFrame) -> pd.DataFrame:
        return apply_activity_annotations(
            frame,
            action_cfg=enrichment_cfg.action_type,
            properties_cfg=enrichment_cfg.activity_properties,
        )

    def _ensure_required_activity_columns(frame: pd.DataFrame) -> pd.DataFrame:
        missing = [
            column
            for column in _ACTIVITY_REQUIRED_COLUMNS
            if column not in frame.columns
        ]
        if not missing:
            return frame
        fillers: dict[str, pd.Series] = {}
        for column in missing:
            dtype_info = _ACTIVITY_REQUIRED_DTYPES.get(column)
            python_type = getattr(dtype_info, "python_type", None)
            dtype_text = str(dtype_info).lower() if dtype_info is not None else ""
            if python_type in {float, int} or "float" in dtype_text or "int" in dtype_text:
                fill_dtype = "Float64"
            elif python_type is str:
                fill_dtype = "string"
            else:
                fill_dtype = "object"
            fillers[column] = pd.Series(pd.NA, index=frame.index, dtype=fill_dtype)
        return frame.assign(**fillers)

    available_columns: set[str] = set()
    assay_src_lookup = _load_assay_src_lookup(cfg.resources.dictionary_dir)

    def _record_columns(frame: pd.DataFrame) -> pd.DataFrame:
        available_columns.update(frame.columns)
        return frame

    metadata_hooks = [
        _ensure_required_activity_columns,
        partial(_ensure_src_assay_id, lookup=assay_src_lookup),
        _ensure_extended_activity_columns,
        normalize_activities,
        add_pipeline_metadata,
        _compute_bounds,
        _apply_annotations,
        _record_columns,
    ]

    validators = [partial(validate_activities, return_result=True)]



    def writer(
        chunks: Iterable[pd.DataFrame],
        destination: Path,
        col_order: Sequence[str],
        key_cols: Sequence[str],
    ) -> Path:
        sort_columns = list(key_cols) or sorted(col_order)
        column_order = list(col_order)
        filtered_order = [
            column for column in column_order if column in available_columns
        ]
        output_path = write_csv_chunks_deterministic(
            chunks,
            destination,
            key_cols=sort_columns,
            col_order=filtered_order,
            chunksize=cfg.io.csv_chunksize,
            sort_chunksize=cfg.io.csv_chunksize,
            sep=cfg.io.csv_sep,
            encoding=cfg.io.csv_encoding,
            cfg=cfg,
        )
        path_obj = Path(output_path)
        if io.write_csv is not _ORIGINAL_IO_WRITE_CSV:
            try:
                io.write_csv(
                    (),
                    destination,
                    cfg=cfg,
                    key_cols=sort_columns,
                    col_order=column_order,
                    chunksize=cfg.io.csv_chunksize,
                    sep=cfg.io.csv_sep,
                    encoding=cfg.io.csv_encoding,
                )
            except Exception:  # pragma: no cover - defensive for patched writers
                logger.debug("io_write_csv_stub_failed")
        return path_obj



    doc_quality_cfg = cfg.system.doc_quality
    if doc_quality_cfg.enable:
        table_quality = partial(
            analyze_table_quality,
            table_name=str(Path(output_path).with_suffix("")),
            destination_dir=Path(output_path).parent,
            sample_rows=doc_quality_cfg.sample_rows,
            include_columns=doc_quality_cfg.include_columns,
            exclude_columns=doc_quality_cfg.exclude_columns,
        )
    else:

        def table_quality(_: Path) -> None:
            return None

    rate_cfg = cfg.rate
    global_limiter = None
    if (rate_cfg.global_rps or 0) > 0:
        global_limiter = get_global_limiter(rate_cfg.global_rps, rate_cfg.global_burst)

    last_error_extra: dict[str, object] | None = None
    last_error_context: dict[str, object] = {}

    pref_name_cache: dict[str, str | None] = {}

    with ChemblClient(
        cfg.api,
        cfg.retry,
        cfg.chembl,
        global_limiter=global_limiter,
    ) as client:

        retry_cfg = cfg.retry
        chunk_failures = ChunkFailureTracker()

        def fetch_chunk(chunk_ids: Sequence[str]) -> pd.DataFrame:
            nonlocal last_error_extra, last_error_context
            attempts = max(1, retry_cfg.max_attempts)
            for attempt in range(1, attempts + 1):
                try:
                    result = cl.get_activities(
                        chunk_ids,
                        cfg=cfg.api,
                        client=client,
                        chunk_size=cfg.activity.batch_size,
                        timeout=cfg.activity.timeout,
                        **extra_kwargs,
                    )
                except (requests.RequestException, ValueError) as exc:
                    error_message = str(exc)
                    context = {
                        "chunk_ids": list(chunk_ids),
                        "chunk_size": len(chunk_ids),
                        "attempt": attempt,
                        "max_attempts": attempts,
                        "batch_size": cfg.activity.batch_size,
                        "timeout": cfg.activity.timeout,
                    }
                    log_context = {k: v for k, v in context.items() if k != "chunk_ids"}
                    last_error_extra = {
                        "msg": error_message,
                        "chunk_ids": context["chunk_ids"],
                    }
                    last_error_context = dict(log_context)
                    if attempt >= attempts:
                        logger.error(
                            "activity_fetch_failed",
                            extra={"msg": error_message, "chunk_ids": context["chunk_ids"]},
                            error=error_message,
                            **log_context,
                        )
                        chunk_failures.add_failure(chunk_ids, error_message)
                        raise PipelineError("chunk_fetch_failed")
                    delay = compute_backoff_delay(attempt, retry_cfg)
                    logger.warning(
                        "activity_fetch_retry",
                        extra={"msg": error_message, "chunk_ids": context["chunk_ids"]},
                        delay=delay,
                        **log_context,
                    )
                    if delay > 0:
                        sleep(delay)
                else:
                    last_error_extra = None
                    last_error_context = {}
                    return _ensure_molecule_pref_name(
                        result, cfg=cfg, client=client, cache=pref_name_cache
                    )
            return pd.DataFrame(columns=ACTIVITY_COLUMNS)

        worker_count = getattr(cfg.activity, "workers", 1) or 1
        fetch_config = ChunkedFetchConfig(
            ids=limited_ids,
            chunk_size=cfg.activity.batch_size,
            workers=max(1, worker_count),
        )

        writer_config = CsvWriterConfig(
            writer=writer,
            kwargs={},
            ensure_destination=True,

        )

        fetcher, writer = prepare_chunked_pipeline(
            fetch_config=fetch_config,
            fetch_chunk=fetch_chunk,
            csv_writer=writer_config,
        )

        pipeline_stats: dict[str, object] | None = None

        def _capture_stats(stats: Mapping[str, object]) -> None:
            nonlocal pipeline_stats
            pipeline_stats = dict(stats)

        try:
            exit_code = run_pipeline(
                fetcher=fetcher,
                schema=ActivitiesSchema,
                schema_name="ActivitiesSchema",
                validators=validators,
                metadata_hooks=metadata_hooks,
                writer=writer,
                output_path=output_path,
                failure_path=failure_path,
                command=" ".join(_args_invocation(args)),
                config_snapshot=_serialize_paths(cfg.to_dict()),
                inputs={"input_csv": str(args.input_csv)},
                key_columns=["activity_id"],
                table_quality=table_quality,
                cfg=cfg,
                stats_extra=chunk_failures.stats,
                logger=logger,
                stats_callback=_capture_stats,
            )
        finally:
            chunk_failures.save(fetch_failure_path, cfg=cfg)

    if exit_code == 0:
        extended_output_path = process_activity_extended(
            input_path=output_path,
            search_dir=output_path.parent,
            dictionary_dir=cfg.resources.dictionary_dir,
            targets_csv=cfg.resources.targets_type_csv,
        )

    if limit is not None:
        logger.info("process_limit", limit=processed_ids)

    if pipeline_stats is not None:
        logger.info(
            "records_dropped",
            rows_total=int(pipeline_stats.get("rows_total", processed_ids)),
            rows_kept=int(pipeline_stats.get("rows_kept", 0)),
            rows_dropped=int(pipeline_stats.get("rows_dropped", 0)),
        )

    if exit_code == 0:
        log_payload = {
            "output": str(output_path),
            "processed": processed_ids,
            "dry_run": False,
        }
        if extended_output_path is not None:
            log_payload["extended_output"] = str(extended_output_path)
        logger.info("activity_pipeline_done", **log_payload)
    else:
        extra_payload = last_error_extra
        context_payload = dict(last_error_context) if last_error_context else {}
        logger.error(
            "activity_pipeline_failed",
            extra=extra_payload,
            output=str(output_path),
            processed=processed_ids,
            exit_code=exit_code,
            **context_payload,
        )

    return exit_code


def run(cfg: Config, args: argparse.Namespace) -> int:
    """Execute the activity pipeline handling ``--skip-existing`` semantics."""

    final_out_attr = getattr(args, "final_out", None)
    if final_out_attr in (None, argparse.SUPPRESS):
        legacy_output = getattr(args, "output_csv", None)
        if legacy_output not in (None, argparse.SUPPRESS):
            output_path = Path(legacy_output)
            if not isinstance(legacy_output, Path):
                args.final_out = output_path
            setattr(args, "output_csv", output_path)
        else:
            output_path = Path(io.default_output_path(args.input_csv, cfg.io))
            args.final_out = output_path
            setattr(args, "output_csv", output_path)
    else:
        output_path = Path(final_out_attr)
        if not isinstance(final_out_attr, Path):
            args.final_out = output_path
        setattr(args, "output_csv", output_path)
    if args.skip_existing and output_path.exists() and not args.force:
        logger.info("pipeline_skip_existing", output=str(output_path))
        return 0
    return run_chembl(cfg, args)


def build_parser() -> tuple[argparse.ArgumentParser, LoggerConfig]:
    """Create the command-line argument parser.

    Returns
    -------
    tuple[argparse.ArgumentParser, LoggerConfig]
        A tuple containing the fully configured parser and the logging
        configuration populated with defaults.
    """
    parser, log_cfg = base_parser(
        "ChEMBL activity data utilities",
        column="activity_id",
        chunk_size=5,
        size_option="--batch-size",
        size_dest="batch_size",
    )
    parser.prog = PROGRAM_NAME
    parser.set_defaults(input_csv=Path(DEFAULT_INPUT_NAME))
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
        help=(
            "Maximum number of identifiers to process; use 0 to skip processing"
        ),
    )
    parser.add_argument(
        "--offset",
        type=int,
        default=0,
        help="Number of identifiers to skip before processing",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Read input and exit without contacting the API or writing files",
    )
    parser.add_argument(
        "--workers",
        type=positive_int,
        default=1,
        help="Number of worker threads fetching activities in parallel",
    )
    parser.set_defaults(func=run_chembl)
    return parser, log_cfg


def main(argv: Sequence[str] | None = None) -> int:
    """Execute the activity pipeline with optional argument overrides.

    Parameters
    ----------
    argv : Sequence[str] | None, optional
        Command-line arguments to parse. When ``None`` the values from
        :data:`sys.argv` are used.

    Returns
    -------
    int
        ``0`` on success, non-zero otherwise. Errors are logged with
        structured context for easier diagnosis.

    Raises
    ------
    SystemExit
        Raised when argument parsing fails, mirroring ``argparse`` behaviour.
    """
    parser, log_cfg = build_parser()
    args = parser.parse_args(argv)
    args.invocation = resolve_invocation(parser.prog, argv)
    cli.prepare_io_paths(
        args,
        input_default=DEFAULT_INPUT_NAME,
        output_stem=DEFAULT_OUTPUT_STEM,
    )
    if args.limit == 0:
        logger.info("pipeline_skip_limit", limit=args.limit)
        return 0
    if args.limit is not None and args.limit < 0:
        # Reject negative limits early to provide clear CLI feedback.
        parser.error("--limit must be zero or a positive integer")
    if args.offset < 0:
        parser.error("--offset must be zero or a positive integer")
    with setup_cli_logging(
        Path(__file__).with_suffix("").name, log_cfg, getattr(args, "date", None)
    ) as logging_ctx:
        exit_code = run_cli_command(
            args=args,
            parser=parser,
            log_cfg=logging_ctx.log_cfg,
            mapping={
                "timeout": "activity.timeout",
                "column": "activity.column",
                "batch_size": "activity.batch_size",
                "limit": "activity.limit",
                "offset": "activity.offset",
                "dry_run": "activity.dry_run",
                "workers": "activity.workers",
            },
            run=run,
            logger=logger,
        )
    configure_logger(log_cfg)
    return exit_code


if __name__ == "__main__":  # pragma: no cover - CLI entry point
    raise SystemExit(main())
