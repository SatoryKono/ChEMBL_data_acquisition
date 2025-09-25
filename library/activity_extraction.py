"""High level helpers for downloading ChEMBL activity data.

The :func:`extract_activities` function orchestrates the full pipeline used in
other scripts of the project.  It reads activity identifiers from a CSV file,
queries the public ChEMBL API in batches and persists the validated result in a
new CSV file alongside metadata artefacts.
"""

from __future__ import annotations

import logging
from collections.abc import Iterable, Sequence
from itertools import islice
from pathlib import Path
from typing import Final

import requests
from pandera.errors import SchemaErrors

from library import chembl_library as cl
from library import io
from library.chembl_client import ChemblClient
from library.config import Config, IoCfg, _serialize_paths
from library.metadata import Stats, file_sha256, write_meta_yaml
from library.sidecar import SidecarErrors
from library.table_quality import analyze_table_quality
from library.validation import validate_activities
from schemas import ActivitiesSchema, normalize_activities

from .pipeline_metadata import add_pipeline_metadata

LOG: Final = logging.getLogger(__name__)


def _prepare_identifiers(
    *,
    input_csv: Path,
    column: str,
    cfg: IoCfg,
    limit: int | None,
    sep: str | None,
    encoding: str | None,
    na_markers: Sequence[str] | None,
) -> Iterable[str] | list[str]:
    """Read identifiers from ``input_csv`` applying ``limit`` when provided."""

    ids_iter: Iterable[str] = io.read_ids(
        input_csv,
        column=column,
        cfg=cfg,
        sep=sep,
        encoding=encoding,
        na_markers=na_markers,
    )
    if limit is None:
        return ids_iter
    limited = list(islice(ids_iter, limit))
    LOG.info("Limiting processing to the first %d identifiers", len(limited))
    return limited


def extract_activities(
    input_csv: Path | str,
    output_csv: Path | str | None,
    *,
    column: str = "activity_id",
    chunk_size: int = 5,
    timeout: float = 30.0,
    limit: int | None = None,
    sep: str | None = None,
    encoding: str | None = None,
    na_markers: Sequence[str] | None = None,
    cfg: Config | None = None,
    command: str | None = None,
) -> int:
    """Retrieve activity data from ChEMBL and write the validated CSV.

    Parameters
    ----------
    input_csv:
        CSV file containing a column with activity identifiers.
    output_csv:
        Destination for the resulting CSV. If ``None`` the default output path
        derived from :class:`~library.config.IoCfg` is used.
    column:
        Name of the column holding ``activity_id`` values.
    chunk_size:
        Maximum number of identifiers per HTTP request.
    timeout:
        Request timeout in seconds for API calls.
    limit:
        Optional upper bound on the number of identifiers read from the input
        CSV. ``None`` processes the entire file.
    sep:
        Optional field delimiter overriding the configuration defaults.
    encoding:
        Optional file encoding overriding the configuration defaults.
    na_markers:
        Optional list of markers interpreted as missing values when reading the
        input CSV.
    cfg:
        Optional :class:`~library.config.Config` instance providing runtime
        defaults. When ``None`` a new configuration with default values is
        created.
    command:
        String stored in the metadata sidecar to describe the invocation. When
        omitted it defaults to ``"extract_activities"``.

    Returns
    -------
    int
        ``0`` on success, ``1`` when a recoverable error occurred during
        retrieval, validation or file writing.
    """

    if chunk_size <= 0:
        LOG.error("chunk_size must be positive, got %d", chunk_size)
        return 1
    if timeout < 0:
        LOG.error("timeout must be non-negative, got %s", timeout)
        return 1
    if limit is not None and limit < 0:
        LOG.error("limit must be non-negative, got %d", limit)
        return 1

    cfg = cfg or Config()
    io_cfg = cfg.io
    input_path = Path(input_csv)
    ids: Iterable[str] | list[str]
    try:
        ids = _prepare_identifiers(
            input_csv=input_path,
            column=column,
            cfg=io_cfg,
            limit=limit,
            sep=sep,
            encoding=encoding,
            na_markers=na_markers,
        )
    except (FileNotFoundError, ValueError) as exc:
        LOG.error("%s", exc)
        return 1

    output_path = (
        Path(output_csv)
        if output_csv is not None
        else io.default_output_path(input_path, io_cfg)
    )

    with ChemblClient(cfg.api, cfg.retry, cfg.chembl) as client:
        try:
            df = cl.get_activities(
                ids,
                cfg=cfg.api,
                client=client,
                chunk_size=chunk_size,
                timeout=timeout,
            )
        except (requests.RequestException, ValueError) as exc:
            LOG.error("failed to retrieve activities: %s", exc)
            return 1

    df = normalize_activities(df)
    df = add_pipeline_metadata(df)

    required = {name for name, col in ActivitiesSchema.columns.items() if col.required}
    missing_required = required - set(df.columns)
    exit_code = 0
    if missing_required:
        LOG.warning(
            "Skipping validation; missing columns: %s", sorted(missing_required)
        )
    else:
        try:
            validation_result = validate_activities(df, return_result=True)
        except SchemaErrors as exc:
            sidecar = SidecarErrors()
            for row in exc.failure_cases.to_dict("records"):
                sidecar.add_error(row)
            failure = Path(output_path).with_name(
                f"{Path(output_path).stem}_failure_cases.csv"
            )
            sidecar.save(failure, cfg=cfg)
            LOG.error("validation failed; see %s", failure)
            df = getattr(exc, "validated_data", df)
            exit_code = 1
        else:
            df = validation_result.data
            if not validation_result.failure_cases.empty:
                sidecar = SidecarErrors()
                for row in validation_result.failure_cases.to_dict("records"):
                    sidecar.add_error(row)
                failure = Path(output_path).with_name(
                    f"{Path(output_path).stem}_failure_cases.csv"
                )
                sidecar.save(failure, cfg=cfg)
                LOG.error("validation failed; see %s", failure)
                exit_code = 1

    schema_cols = list(ActivitiesSchema.columns)
    head = [c for c in schema_cols if c in df.columns]
    tail = sorted(c for c in df.columns if c not in schema_cols)
    col_order = head + tail

    key_cols = [col for col in ["activity_id"] if col in df.columns]

    try:
        csv_path = io.write_csv(
            df,
            output_path,
            cfg=cfg,
            key_cols=key_cols or None,
            col_order=col_order,
        )
    except OSError as exc:
        LOG.error("failed to write output CSV: %s", exc)
        return 1

    stats: Stats = {
        "rows_total": len(df),
        "rows_kept": len(df),
        "rows_dropped": 0,
        "output_sha256": file_sha256(csv_path),
    }

    write_meta_yaml(
        csv_path=csv_path,
        command=command or "extract_activities",
        config_subset=_serialize_paths(cfg.to_dict()),
        inputs={"input_csv": str(input_csv)},
        stats=stats,
        schema="ActivitiesSchema",
    )

    try:
        analyze_table_quality(df, table_name=str(Path(csv_path).with_suffix("")))
    except ValueError as exc:
        LOG.error("failed to generate quality report: %s", exc)
        return 1

    return exit_code
