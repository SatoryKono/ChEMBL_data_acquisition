"""I/O helpers for CSV-based data pipelines.

This module centralises reading and writing of CSV files to ensure
consistent handling of delimiters, encodings and error reporting across
command line utilities. Writing functions automatically create parent
directories and generate YAML sidecar files describing the runtime
environment and configuration.

Examples
--------
>>> from pathlib import Path
>>> import pandas as pd
>>> from library.config import Config
>>> from library.io import write_csv, read_csv
>>> cfg = Config()  # doctest: +SKIP
>>> df = pd.DataFrame({"id": [1]})
>>> out = write_csv(df, Path("same_document/example.csv"), cfg=cfg)  # doctest: +SKIP
>>> read_csv(out, cfg=cfg).to_dict("list")  # doctest: +SKIP
{'id': [1]}
"""

from __future__ import annotations

import csv
import locale
import sys
from collections.abc import Hashable, Iterable, Iterator, Mapping, Sequence
from datetime import datetime
from importlib import import_module
from pathlib import Path
from typing import TYPE_CHECKING, Any

import pandas as pd
import yaml

from . import validation
from .config import Config, IoCfg, _serialize_paths
from .git_utils import _git_sha
from .log import logger


class _EncodingDecodeError(Exception):
    """Wrap :class:`UnicodeDecodeError` with the attempted encoding."""

    def __init__(self, encoding: str, error: UnicodeDecodeError) -> None:
        super().__init__(str(error))
        self.encoding = encoding
        self.error = error


def _stream_ids_with_encoding(
    path: Path,
    *,
    column: str,
    sep: str,
    encoding: str,
    marker_set: set[str],
) -> Iterator[str]:
    """Yield identifier values using a specific ``encoding``."""

    try:
        with path.open("r", encoding=encoding, newline="") as fh:
            reader = csv.DictReader(fh, delimiter=sep)
            if reader.fieldnames is None or column not in reader.fieldnames:
                raise ValueError(
                    f"column '{column}' not found in {path}; available columns: {reader.fieldnames}"
                )
            for row in reader:
                value = (row.get(column) or "").strip()
                if value and value not in marker_set:
                    yield value
    except UnicodeDecodeError as exc:  # pragma: no cover - exercised via fallback tests
        raise _EncodingDecodeError(encoding, exc) from exc

if TYPE_CHECKING:  # pragma: no cover - only for type checking
    import pandera as pa
else:  # pragma: no cover - exercised in tests via monkeypatch
    try:
        pa = import_module("pandera")
    except (ImportError, TypeError):
        pa = None


def read_ids(
    path: str | Path,
    *,
    column: str,
    cfg: IoCfg,
    sep: str | None = None,
    encoding: str | None = None,
    na_markers: Sequence[str] | None = None,
) -> Iterator[str]:
    """Yield identifier values from ``column`` in ``path``.

    Parameters
    ----------
    path:
        Location of the CSV file.
    column:
        Name of the column that contains identifiers.
    cfg:
        I/O configuration providing default CSV parameters.
    sep:
        Field delimiter used in the CSV file. Defaults to ``cfg.csv_sep``.
    encoding:
        Character encoding of the CSV file. Defaults to ``cfg.csv_encoding``.
    na_markers:
        Strings indicating missing values. Defaults to ``cfg.na_markers``.

    Yields
    ------
    str
        Identifier values in the order they appear. Empty strings and values
        present in ``na_markers`` are discarded.

    Raises
    ------
    FileNotFoundError
        If ``path`` does not exist.
    ValueError
        If the CSV file is malformed or ``column`` is missing.

    """
    sep = sep or cfg.csv_sep
    encoding = encoding or cfg.csv_encoding
    marker_set = set(na_markers or cfg.na_markers or ())

    def _append_candidate(values: Sequence[str] | str | None, seen: set[str], out: list[str]) -> None:
        if values is None:
            return
        if isinstance(values, str):
            iterable: Sequence[str] = (values,)
        else:
            iterable = values
        for value in iterable:
            if not value:
                continue
            key = value.lower()
            if key in seen:
                continue
            out.append(value)
            seen.add(key)

    candidates: list[str] = []
    seen_candidates: set[str] = set()
    _append_candidate((encoding,) if encoding else None, seen_candidates, candidates)
    fallbacks = getattr(cfg, "csv_fallback_encodings", ()) or ()
    _append_candidate(fallbacks, seen_candidates, candidates)

    locale_encoding = locale.getpreferredencoding(False)
    _append_candidate(locale_encoding, seen_candidates, candidates)

    if not candidates:
        candidates.append("utf-8")
        seen_candidates.add("utf-8")

    path_obj = Path(path)

    def _iter_candidates() -> Iterator[str]:
        last_error: UnicodeDecodeError | None = None
        for candidate in candidates:
            try:
                yield from _stream_ids_with_encoding(
                    path_obj,
                    column=column,
                    sep=sep,
                    encoding=candidate,
                    marker_set=marker_set,
                )
                return
            except _EncodingDecodeError as exc:
                last_error = exc.error
                logger.warning(
                    "csv_decode_failed",
                    path=str(path_obj),
                    encoding=exc.encoding,
                    error=str(exc.error),
                )
                continue
            except LookupError as exc:
                logger.warning(
                    "csv_encoding_lookup_failed",
                    path=str(path_obj),
                    encoding=candidate,
                    error=str(exc),
                )
                continue
        attempted = ", ".join(candidates)
        message = (
            f"failed to decode CSV {path_obj} with encodings: {attempted}. "
            "Update 'io.csv_encoding' or 'io.csv_fallback_encodings' in the configuration."
        )
        raise ValueError(message) from last_error

    try:
        yield from _iter_candidates()
    except FileNotFoundError:
        raise
    except csv.Error as exc:
        raise ValueError(f"malformed CSV in file: {path}: {exc}") from exc


def read_csv(
    path: str | Path,
    *,
    cfg: IoCfg,
    sep: str | None = None,
    encoding: str | None = None,
    required_columns: Iterable[str] | None = None,
    dtype: Mapping[Hashable, Any] | type | None = None,
    na_values: Sequence[str] | str | None = None,
    parse_dates: Sequence[str] | None = None,
    schema: pa.DataFrameSchema | type[pa.DataFrameModel] | None = None,
) -> pd.DataFrame:
    """Load a CSV file into a :class:`pandas.DataFrame` with optional schema validation.

    Parameters
    ----------
    path:
        Location of the CSV file.
    cfg:
        I/O configuration providing default CSV parameters.
    sep:
        Field delimiter used in the CSV file. Defaults to ``cfg.csv_sep``.
    encoding:
        Character encoding of the CSV file. Defaults to ``cfg.csv_encoding``.
    required_columns:
        Optional list of column names that must be present in the loaded
        DataFrame. A :class:`ValueError` is raised if any are missing.
    dtype:
        Data type specification forwarded to :func:`pandas.read_csv`.
        Can be a single type applied to all columns or a mapping of
        column names to types.
    na_values:
        Additional strings to recognize as NA/NaN. Passed to
        :func:`pandas.read_csv`.
    parse_dates:
        Column names to parse as dates using :func:`pandas.read_csv`.
    schema:
        Optional :class:`pa.DataFrameSchema` or
        :class:`pa.DataFrameModel` used for advanced validation and
        dtype coercion. Requires :mod:`pandera` when provided.

    Returns
    -------
    pandas.DataFrame
        DataFrame containing the CSV contents.

    """
    sep = sep or cfg.csv_sep
    encoding = encoding or cfg.csv_encoding
    df = pd.read_csv(
        path,
        sep=sep,
        encoding=encoding,
        dtype=dtype,
        na_values=na_values,
        parse_dates=list(parse_dates) if parse_dates is not None else None,
    )
    if schema is not None:
        if pa is None:
            raise RuntimeError(
                "pandera is required for schema validation; install pandera to use the 'schema' argument"
            )
        if isinstance(schema, pa.DataFrameSchema):
            df = schema.validate(df)
        else:
            df = schema.to_schema().validate(df)
    elif required_columns is not None:
        validation.validate_columns(df, required_columns)
    return df


def write_csv(
    df: pd.DataFrame,
    path: str | Path,
    *,
    cfg: Config,
    sep: str | None = None,
    encoding: str | None = None,
    key_cols: Iterable[str] | None = None,
    col_order: Iterable[str] | None = None,
    chunksize: int | None = None,
) -> Path:
    """Write ``df`` to ``path`` as CSV and store metadata.

    This is a thin wrapper around :func:`library.csv_utils.write_csv_deterministic`
    which handles deterministic sorting and metadata sidecar creation.
    Parent directories of ``path`` are created automatically.

    Parameters
    ----------
    df:
        DataFrame to serialise.
    path:
        Destination file path. Nested paths are allowed and will be
        created as needed.
    cfg:
        Full configuration used for metadata sidecars.
    sep:
        Field delimiter used in the CSV file. Defaults to ``cfg.io.csv_sep``.
    encoding:
        Character encoding of the CSV file. Defaults to ``cfg.io.csv_encoding``.
    key_cols:
        Columns used to determine row order. Defaults to all columns if
        omitted.
    col_order:
        Optional preferred column ordering. Columns not listed here are
        appended in lexicographical order.
    chunksize:
        Optional number of rows per chunk to stream when writing the CSV.
        Passed through to :meth:`pandas.DataFrame.to_csv` and
        :func:`library.csv_utils.write_csv_deterministic`.

    Examples
    --------
    >>> from pathlib import Path
    >>> df = pd.DataFrame({"id": [1]})
    >>> cfg = Config()  # doctest: +SKIP
    >>> write_csv(df, Path("independent/example.csv"), cfg=cfg)  # doctest: +SKIP
    PosixPath('independent/example.csv')

    Returns
    -------
    pathlib.Path
        Path to the written CSV file.
    """
    sep = sep or cfg.io.csv_sep
    encoding = encoding or cfg.io.csv_encoding
    key_cols_list = list(key_cols) if key_cols is not None else sorted(df.columns)
    missing_keys = [col for col in key_cols_list if col not in df.columns]
    if missing_keys:
        logger.error(
            "missing_key_columns",
            requested=key_cols_list,
            missing=missing_keys,
        )
        raise ValueError(f"Missing key columns: {missing_keys}")
    col_order_list = list(col_order) if col_order is not None else None
    from .csv_utils import write_csv_deterministic

    return write_csv_deterministic(
        df,
        path,
        key_cols=key_cols_list,
        col_order=col_order_list,
        chunksize=chunksize,
        sep=sep,
        encoding=encoding,
        cfg=cfg,
    )


def default_output_path(input_path: str | Path, cfg: IoCfg) -> Path:
    """Return the default output path for ``input_path``.

    Parameters
    ----------
    input_path:
        Source file path used only for deriving the stem name.
    cfg:
        I/O configuration providing the base directory for generated files.

    Returns
    -------
    pathlib.Path
        Path pointing inside ``cfg.output_dir`` following the
        ``output_<stem>_YYYYMMDD.csv`` naming scheme.
    """
    inp = Path(input_path)
    date_str = datetime.now().strftime("%Y%m%d")
    return Path(cfg.output_dir) / f"output_{inp.stem}_{date_str}.csv"


def write_meta_yaml(
    path: Path | str,
    cfg: Config | None = None,
    *,
    columns: Sequence[str] | None = None,
    dtypes: Mapping[str, str] | None = None,
) -> Path:
    """Write metadata for ``path`` to ``<path>.meta.yaml``.

    Parameters
    ----------
    path:
        Target file for which the sidecar should be created.
    cfg:
        Optional configuration instance. When provided all :class:`~pathlib.Path`
        objects inside the configuration are serialised as plain strings to make
        the resulting YAML portable. When ``None`` an empty mapping is written
        under the ``config`` key.
    columns:
        Optional sequence of column names present in the exported table.
    dtypes:
        Optional mapping from column names to their data types. When omitted
        and ``columns`` are provided each column is assumed to contain
        ``"string"`` values.

    Returns
    -------
    pathlib.Path
        Path to the generated ``.meta.yaml`` file.
    """

    if dtypes is None and columns is not None:
        dtypes = {col: "string" for col in columns}

    meta = {
        "generated_at": datetime.now().isoformat(),
        "git_sha": _git_sha(),
        "command": " ".join(sys.argv),
        "columns": list(columns or []),
        "dtypes": dict(dtypes or {}),
        "config": _serialize_paths(cfg.to_dict()) if cfg is not None else {},
    }
    meta_path = Path(f"{path}.meta.yaml")
    with meta_path.open("w", encoding="utf8") as fh:
        yaml.safe_dump(meta, fh, sort_keys=False)
    return meta_path
