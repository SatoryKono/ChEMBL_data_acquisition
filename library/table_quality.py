"""Tabular data profiling utilities.

This module provides the :func:`analyze_table_quality` function which
profiles every column in a given :class:`pandas.DataFrame` and computes
pairwise correlations for numeric columns.  It is intended for offline
use without external dependencies beyond :mod:`pandas` and :mod:`numpy`.
"""

from __future__ import annotations

import re
import warnings
from collections import Counter
from collections.abc import Iterable, Sequence, Sized
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from pandas.errors import DtypeWarning

from .common.log import logger

# Precompiled regular expressions for pattern coverage
_DOI_RE = re.compile(r"^10\.\d{4,9}/\S+$")
_ISSN_RE = re.compile(r"^\d{4}-\d{3}[\dXx]$")
_URL_RE = re.compile(r"^(?:https?|ftp)://|^www\.")
_EMAIL_RE = re.compile(r"^[^@\s]+@[^@\s]+\.[^@\s]+$")

BOOL_LIKE = {"true", "false", "yes", "no", "y", "n", "1", "0", "t", "f"}


def _load_table(table: pd.DataFrame | str | Path) -> pd.DataFrame:
    """Load a table from ``table``.

    Parameters
    ----------
    table:
        Either an existing :class:`pandas.DataFrame` or path to a CSV file.

    Returns
    -------
    pandas.DataFrame

    """
    # Importing pandas locally guards against environments where the global
    # ``pd`` alias might be missing.  This avoids ``NameError`` exceptions when
    # :func:`analyze_table_quality` is executed in constrained runtimes.
    import pandas as pd

    if isinstance(table, pd.DataFrame):
        return table.copy()

    path = Path(table)
    encodings = ["utf-8-sig", "utf-8", "cp1251", "latin-1"]
    for enc in encodings:
        try:
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore", category=DtypeWarning)
                return pd.read_csv(path, encoding=enc, low_memory=False)
        except UnicodeDecodeError:
            logger.debug("failed to decode %s with %s", path, enc)
            continue
    raise UnicodeDecodeError("utf-8", b"", 0, 1, "Unable to decode CSV")


def _apply_sampling_and_filters(
    frame: pd.DataFrame,
    *,
    table_name: str,
    sample_rows: int | None,
    include_columns: tuple[str, ...] | None,
    exclude_columns: tuple[str, ...] | None,
    include_warning_logged: list[bool],
    exclude_warning_logged: list[bool],
    no_columns_logged: list[bool],
) -> tuple[pd.DataFrame, int | None]:
    """Apply row sampling and column filters to ``frame``."""

    remaining: int | None = sample_rows
    filtered = frame
    if remaining is not None:
        if remaining <= 0:
            return frame.iloc[0:0], 0
        filtered = filtered.head(remaining)

    if include_columns:
        missing = sorted(set(include_columns) - set(frame.columns))
        if missing and not include_warning_logged[0]:
            logger.warning(
                "include_columns_missing",
                columns=missing,
                table_name=table_name,
            )
            include_warning_logged[0] = True
        include_set = set(include_columns)
        filtered = filtered.loc[
            :, [col for col in filtered.columns if col in include_set]
        ]

    if exclude_columns:
        missing = sorted(set(exclude_columns) - set(frame.columns))
        if missing and not exclude_warning_logged[0]:
            logger.warning(
                "exclude_columns_missing",
                columns=missing,
                table_name=table_name,
            )
            exclude_warning_logged[0] = True
        exclude_set = set(exclude_columns)
        filtered = filtered.loc[
            :, [col for col in filtered.columns if col not in exclude_set]
        ]

    if (
        (include_columns or exclude_columns)
        and filtered.shape[1] == 0
        and not no_columns_logged[0]
    ):
        logger.warning("no_columns_after_filter", table_name=table_name)
        no_columns_logged[0] = True

    if remaining is not None:
        remaining = max(remaining - len(filtered), 0)

    return filtered, remaining


def _is_isbn(value: str) -> bool:
    """Return ``True`` if ``value`` is a valid ISBN10/13."""
    digits = re.sub(r"[-\s]", "", value)
    if re.fullmatch(r"\d{9}[\dXx]", digits):
        total = sum(
            (10 - i) * (10 if ch.upper() == "X" else int(ch))
            for i, ch in enumerate(digits)
        )
        return total % 11 == 0
    if re.fullmatch(r"\d{13}", digits):
        total = sum(
            (1 if i % 2 == 0 else 3) * int(ch) for i, ch in enumerate(digits[:-1])
        )
        check = (10 - total % 10) % 10
        return check == int(digits[-1])
    return False


def _non_empty_mask(series: pd.Series) -> pd.Series:
    """Boolean mask of values that are not empty by content."""

    def is_non_empty(val: Any) -> bool:
        if isinstance(val, str):
            return bool(val.strip())
        if isinstance(val, Sized) and not isinstance(val, bytes | bytearray):
            # For sequences and other sized containers, consider length
            return len(val) > 0
        try:
            return not pd.isna(val)
        except (TypeError, ValueError):
            # Objects that cannot be evaluated by pandas are treated as non-empty
            return True

    return series.map(is_non_empty)


def _string_values(series: pd.Series, mask: pd.Series) -> pd.Series:
    """Return non-empty string values from ``series`` according to ``mask``."""
    return (
        series[mask & series.map(lambda x: isinstance(x, str))].astype(str).str.strip()
    )


def _pattern_cov(strings: pd.Series, pattern: re.Pattern[str]) -> float:
    """Fraction of ``strings`` matching ``pattern``."""
    if strings.empty:
        return 0.0
    return float(strings.str.match(pattern).mean())


def _isbn_cov(strings: pd.Series) -> float:
    if strings.empty:
        return 0.0
    return float(strings.map(_is_isbn).mean())


def _bool_like_cov(values: pd.Series) -> float:
    if values.empty:
        return 0.0
    return float(values.astype(str).str.strip().str.lower().isin(BOOL_LIKE).mean())


def _numeric_stats(series: pd.Series) -> tuple[pd.Series, float, dict[str, float]]:
    """Convert ``series`` to numeric values and compute summary statistics."""
    numeric = pd.to_numeric(series, errors="coerce").astype(float)

    coverage = float(numeric.notna().mean())
    stats = {
        "numeric_min": float(numeric.min()) if coverage else np.nan,
        "numeric_p50": float(numeric.quantile(0.5)) if coverage else np.nan,
        "numeric_p95": float(numeric.quantile(0.95)) if coverage else np.nan,
        "numeric_max": float(numeric.max()) if coverage else np.nan,
        "numeric_mean": float(numeric.mean()) if coverage else np.nan,
        "numeric_std": float(numeric.std(ddof=0)) if coverage else np.nan,
    }

    return numeric, coverage, stats


def _parse_dates(series: pd.Series) -> tuple[float, dict[str, pd.Timestamp | float]]:
    def normalise(val: object) -> object:
        if isinstance(val, str) and re.fullmatch(r"\d{4}", val.strip()):
            return f"{val.strip()}-07-01"
        return val

    normalised = series.map(normalise)

    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message="Parsed string.*included an un-recognized timezone",
            category=FutureWarning,
        )
        dates = pd.to_datetime(normalised, errors="coerce", utc=True, format="mixed")

    coverage = float(dates.notna().mean())
    if coverage:
        dt = dates.dropna().dt.tz_convert(None)
        stats: dict[str, pd.Timestamp | float] = {
            "date_min": dt.min(),
            "date_p50": dt.quantile(0.5),
            "date_max": dt.max(),
        }
    else:
        stats = {"date_min": np.nan, "date_p50": np.nan, "date_max": np.nan}
    return coverage, stats


def _text_length_stats(series: pd.Series) -> dict[str, float]:
    values = series.dropna().map(lambda x: str(x).strip())
    if values.empty:
        return {
            "text_len_min": np.nan,
            "text_len_p50": np.nan,
            "text_len_p95": np.nan,
            "text_len_max": np.nan,
        }
    lengths = values.map(len)
    return {
        "text_len_min": float(lengths.min()),
        "text_len_p50": float(lengths.quantile(0.5)),
        "text_len_p95": float(lengths.quantile(0.95)),
        "text_len_max": float(lengths.max()),
    }


def _top_values(series: pd.Series) -> str:
    counts = series.dropna().map(lambda x: str(x).strip()).value_counts().head(3)
    parts = [f"{str(val)[:60]} ({cnt})" for val, cnt in counts.items()]
    return "; ".join(parts)


class _ColumnAccumulator:
    """Incrementally gather profiling metrics for a single column."""

    __slots__ = (
        "name",
        "total",
        "non_null",
        "non_empty",
        "string_count",
        "pattern_counts",
        "bool_like_matches",
        "unique_hashes",
        "top_counter",
        "numeric_values",
        "numeric_full",
        "numeric_valid",
        "date_values",
        "date_valid",
        "text_lengths",
        "numeric_cov",
    )

    def __init__(self, name: str, prefilled: int = 0) -> None:
        self.name = name
        self.total = 0
        self.non_null = 0
        self.non_empty = 0
        self.string_count = 0
        self.pattern_counts = {
            "doi": 0,
            "issn": 0,
            "isbn": 0,
            "url": 0,
            "email": 0,
        }
        self.bool_like_matches = 0
        self.unique_hashes: set[int] = set()
        self.top_counter: Counter[str] = Counter()
        self.numeric_values: list[float] = []
        self.numeric_full: list[float] = []
        self.numeric_valid = 0
        self.date_values: list[pd.Timestamp] = []
        self.date_valid = 0
        self.text_lengths: list[int] = []
        self.numeric_cov = 0.0
        if prefilled:
            self.pad(prefilled)

    def pad(self, count: int) -> None:
        """Register ``count`` missing values for future chunks."""

        if count <= 0:
            return
        self.total += count
        self.numeric_full.extend([np.nan] * count)

    def process(self, series: pd.Series) -> None:
        """Update statistics based on ``series`` values."""

        length = len(series)
        self.total += length
        if length == 0:
            return

        self.non_null += int(series.notna().sum())

        mask = _non_empty_mask(series)
        non_empty = int(mask.sum())
        self.non_empty += non_empty

        strings = _string_values(series, mask)
        string_len = len(strings)
        if string_len:
            self.string_count += string_len
            self.pattern_counts["doi"] += int(strings.str.match(_DOI_RE).sum())
            self.pattern_counts["issn"] += int(strings.str.match(_ISSN_RE).sum())
            self.pattern_counts["url"] += int(strings.str.match(_URL_RE).sum())
            self.pattern_counts["email"] += int(strings.str.match(_EMAIL_RE).sum())
            isbn_matches = strings.map(_is_isbn)
            self.pattern_counts["isbn"] += int(isbn_matches.sum())
            lower = strings.str.lower()
            bool_mask = lower.isin(BOOL_LIKE)
            self.bool_like_matches += int(bool_mask.sum())

        non_na = series.dropna()
        if not non_na.empty:
            trimmed = non_na.map(lambda x: str(x).strip())
            self.top_counter.update(trimmed)
            self.text_lengths.extend(trimmed.map(len).tolist())
            try:
                hashes = pd.util.hash_pandas_object(non_na, index=False)
            except TypeError:
                hashes = pd.util.hash_pandas_object(
                    non_na.map(lambda x: str(x)), index=False
                )
            self.unique_hashes.update(hashes.astype(np.uint64).tolist())

        numeric = pd.to_numeric(series, errors="coerce").astype(float)
        self.numeric_full.extend(numeric.tolist())
        numeric_valid = numeric.dropna()
        valid_len = len(numeric_valid)
        if valid_len:
            self.numeric_valid += valid_len
            self.numeric_values.extend(numeric_valid.tolist())

        def _normalise(val: object) -> object:
            if isinstance(val, str) and re.fullmatch(r"\d{4}", val.strip()):
                return f"{val.strip()}-07-01"
            return val

        normalised = series.map(_normalise)
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore",
                message="Parsed string.*included an un-recognized timezone",
                category=FutureWarning,
            )
            dates = pd.to_datetime(
                normalised, errors="coerce", utc=True, format="mixed"
            )

        valid_dates = dates.dropna()
        if not valid_dates.empty:
            naive = valid_dates.dt.tz_convert(None)
            self.date_values.extend(naive.tolist())
            self.date_valid += len(naive)

    def finalize(self) -> dict[str, object]:
        """Return a summary row mirroring ``analyze_table_quality`` output."""

        empty_pct = float(1 - self.non_empty / self.total) if self.total else 0.0
        unique_cnt = len(self.unique_hashes)
        unique_pct = float(unique_cnt / self.non_empty) if self.non_empty else np.nan

        doi_cov = (
            float(self.pattern_counts["doi"] / self.string_count)
            if self.string_count
            else 0.0
        )
        issn_cov = (
            float(self.pattern_counts["issn"] / self.string_count)
            if self.string_count
            else 0.0
        )
        isbn_cov = (
            float(self.pattern_counts["isbn"] / self.string_count)
            if self.string_count
            else 0.0
        )
        url_cov = (
            float(self.pattern_counts["url"] / self.string_count)
            if self.string_count
            else 0.0
        )
        email_cov = (
            float(self.pattern_counts["email"] / self.string_count)
            if self.string_count
            else 0.0
        )
        bool_like_cov = (
            float(self.bool_like_matches / self.string_count)
            if self.string_count
            else 0.0
        )

        numeric_cov = float(self.numeric_valid / self.total) if self.total else 0.0
        self.numeric_cov = numeric_cov
        if self.numeric_values:
            numeric_arr = np.array(self.numeric_values, dtype=float)
            numeric_min = float(np.min(numeric_arr))
            numeric_p50 = float(np.quantile(numeric_arr, 0.5))
            numeric_p95 = float(np.quantile(numeric_arr, 0.95))
            numeric_max = float(np.max(numeric_arr))
            numeric_mean = float(np.mean(numeric_arr))
            numeric_std = float(np.std(numeric_arr, ddof=0))
        else:
            numeric_min = np.nan
            numeric_p50 = np.nan
            numeric_p95 = np.nan
            numeric_max = np.nan
            numeric_mean = np.nan
            numeric_std = np.nan

        date_cov = float(self.date_valid / self.total) if self.total else 0.0
        if self.date_values:
            date_series = pd.Series(self.date_values)
            date_min = date_series.min()
            date_p50 = date_series.quantile(0.5)
            date_max = date_series.max()
        else:
            date_min = np.nan
            date_p50 = np.nan
            date_max = np.nan

        if self.text_lengths:
            lengths = np.array(self.text_lengths, dtype=float)
            text_len_min = float(lengths.min())
            text_len_p50 = float(np.quantile(lengths, 0.5))
            text_len_p95 = float(np.quantile(lengths, 0.95))
            text_len_max = float(lengths.max())
        else:
            text_len_min = np.nan
            text_len_p50 = np.nan
            text_len_p95 = np.nan
            text_len_max = np.nan

        roles: list[str] = []
        if url_cov > 0:
            roles.append("url")
        if doi_cov > 0:
            roles.append("doi")
        if issn_cov > 0:
            roles.append("issn")
        if isbn_cov > 0:
            roles.append("isbn")
        if email_cov > 0:
            roles.append("email")
        if bool_like_cov >= 0.8:
            roles.append("boolean")
        if date_cov >= 0.8:
            roles.append("date")
        if numeric_cov >= 0.8:
            roles.append("numeric")
        if self.non_empty and unique_cnt / self.non_empty >= 0.98:
            roles.append("identifier-like")
        elif unique_cnt <= min(100, 0.05 * self.non_empty):
            roles.append("categorical")
        else:
            roles.append("free-text")

        parts = [
            f"{value[:60]} ({count})"
            for value, count in self.top_counter.most_common(3)
        ]
        top_values = "; ".join(parts)

        return {
            "column": self.name,
            "non_null": self.non_null,
            "non_empty": self.non_empty,
            "empty_pct": empty_pct,
            "unique_cnt": unique_cnt,
            "unique_pct_of_non_empty": unique_pct,
            "pattern_cov_doi": doi_cov,
            "pattern_cov_issn": issn_cov,
            "pattern_cov_isbn": isbn_cov,
            "pattern_cov_url": url_cov,
            "pattern_cov_email": email_cov,
            "bool_like_cov": bool_like_cov,
            "numeric_cov": numeric_cov,
            "numeric_min": numeric_min,
            "numeric_p50": numeric_p50,
            "numeric_p95": numeric_p95,
            "numeric_max": numeric_max,
            "numeric_mean": numeric_mean,
            "numeric_std": numeric_std,
            "date_cov": date_cov,
            "date_min": date_min,
            "date_p50": date_p50,
            "date_max": date_max,
            "text_len_min": text_len_min,
            "text_len_p50": text_len_p50,
            "text_len_p95": text_len_p95,
            "text_len_max": text_len_max,
            "guessed_roles": "|".join(roles),
            "top_values": top_values,
        }


class TableQualityProfiler:
    """Accumulate quality metrics across streamed ``DataFrame`` chunks."""

    def __init__(self) -> None:
        self._columns: list[str] = []
        self._accumulators: dict[str, _ColumnAccumulator] = {}
        self._rows_processed = 0

    def consume(self, frame: pd.DataFrame) -> None:
        if not isinstance(frame, pd.DataFrame):
            raise TypeError("TableQualityProfiler.consume expects a pandas DataFrame")

        new_columns = [col for col in frame.columns if col not in self._accumulators]
        for column in new_columns:
            accumulator = _ColumnAccumulator(column, prefilled=self._rows_processed)
            self._columns.append(column)
            self._accumulators[column] = accumulator

        missing_columns = [col for col in self._columns if col not in frame.columns]
        for column in missing_columns:
            self._accumulators[column].pad(len(frame))

        processed: set[str] = set()
        for column, series in frame.items():
            if column in processed:
                continue
            processed.add(column)
            self._accumulators[column].process(series)

        self._rows_processed += len(frame)

    def build(
        self,
        table_name: str,
        destination_dir: Path | str | None = None,
    ) -> tuple[pd.DataFrame, pd.DataFrame]:
        rows: list[dict[str, object]] = []
        numeric_candidates: dict[str, pd.Series] = {}
        for column in self._columns:
            accumulator = self._accumulators[column]
            row = accumulator.finalize()
            rows.append(row)
            if accumulator.numeric_cov >= 0.8:
                numeric_candidates[column] = pd.Series(
                    accumulator.numeric_full, dtype=float
                )

        column_order = [
            "column",
            "non_null",
            "non_empty",
            "empty_pct",
            "unique_cnt",
            "unique_pct_of_non_empty",
            "pattern_cov_doi",
            "pattern_cov_issn",
            "pattern_cov_isbn",
            "pattern_cov_url",
            "pattern_cov_email",
            "bool_like_cov",
            "numeric_cov",
            "numeric_min",
            "numeric_p50",
            "numeric_p95",
            "numeric_max",
            "numeric_mean",
            "numeric_std",
            "date_cov",
            "date_min",
            "date_p50",
            "date_max",
            "text_len_min",
            "text_len_p50",
            "text_len_p95",
            "text_len_max",
            "guessed_roles",
            "top_values",
        ]
        quality_report = pd.DataFrame(rows, columns=column_order)
        destination = (
            Path(destination_dir) if destination_dir is not None else Path(".")
        )
        if destination_dir is not None:
            destination.mkdir(parents=True, exist_ok=True)

        table_label = Path(str(table_name)).name
        if not table_label:
            raise ValueError(
                "table_name must contain at least one non-separator character"
            )

        quality_path = destination / f"{table_label}_quality_report_table.csv"
        quality_report.to_csv(quality_path, index=False, encoding="utf-8-sig")

        if numeric_candidates:
            corr_report = pd.DataFrame(numeric_candidates).corr(method="pearson")
        else:
            corr_report = pd.DataFrame()

        corr_path = destination / f"{table_label}_data_correlation_report_table.csv"
        corr_report.reset_index().to_csv(corr_path, index=False, encoding="utf-8-sig")

        return quality_report, corr_report


def analyze_table_quality(
    table: pd.DataFrame | str | Path | Iterable[pd.DataFrame] | TableQualityProfiler,
    table_name: str,
    *,
    destination_dir: Path | str | None = None,
    sample_rows: int | None = None,
    include_columns: Sequence[str] | None = None,
    exclude_columns: Sequence[str] | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Profile ``table`` and compute correlations for numeric columns.

    Parameters
    ----------
    table:
        :class:`pandas.DataFrame` or path to a CSV file.
    table_name:
        Base name used for output files.
    destination_dir:
        Directory where output files are stored. When ``None`` the current
        working directory is used.

    Returns
    -------
    tuple[pandas.DataFrame, pandas.DataFrame]
        Quality report and correlation matrix.

    """
    # Import pandas locally to ensure the alias is available even if the module
    # level import was stripped by the execution environment.
    import pandas as pd

    include_tuple: tuple[str, ...] | None = (
        tuple(include_columns) if include_columns is not None else None
    )
    exclude_tuple: tuple[str, ...] | None = (
        tuple(exclude_columns) if exclude_columns is not None else None
    )

    include_warning_logged = [False]
    exclude_warning_logged = [False]
    no_columns_logged = [False]

    profiler: TableQualityProfiler
    if isinstance(table, TableQualityProfiler):
        profiler = table
    else:
        profiler = TableQualityProfiler()
        if isinstance(table, pd.DataFrame):
            filtered, _ = _apply_sampling_and_filters(
                table,
                table_name=table_name,
                sample_rows=sample_rows,
                include_columns=include_tuple,
                exclude_columns=exclude_tuple,
                include_warning_logged=include_warning_logged,
                exclude_warning_logged=exclude_warning_logged,
                no_columns_logged=no_columns_logged,
            )
            profiler.consume(filtered)
        elif isinstance(table, str | Path):
            df = _load_table(table)
            filtered, _ = _apply_sampling_and_filters(
                df,
                table_name=table_name,
                sample_rows=sample_rows,
                include_columns=include_tuple,
                exclude_columns=exclude_tuple,
                include_warning_logged=include_warning_logged,
                exclude_warning_logged=exclude_warning_logged,
                no_columns_logged=no_columns_logged,
            )
            profiler.consume(filtered)
        elif isinstance(table, Iterable):
            remaining = sample_rows
            for frame in table:
                if not isinstance(frame, pd.DataFrame):
                    raise TypeError(
                        "Streaming quality analysis requires pandas DataFrame chunks"
                    )
                filtered, remaining = _apply_sampling_and_filters(
                    frame,
                    table_name=table_name,
                    sample_rows=remaining,
                    include_columns=include_tuple,
                    exclude_columns=exclude_tuple,
                    include_warning_logged=include_warning_logged,
                    exclude_warning_logged=exclude_warning_logged,
                    no_columns_logged=no_columns_logged,
                )
                profiler.consume(filtered)
                if remaining is not None and remaining <= 0:
                    break
        else:
            raise TypeError(
                "Unsupported table type. Expected DataFrame, path or chunk iterable."
            )

    return profiler.build(table_name, destination_dir=destination_dir)


if __name__ == "__main__":  # pragma: no cover - illustrative usage
    demo = pd.DataFrame(
        {
            "id": [1, 2, 3, 4],
            "url": ["https://example.com", "http://test.com", None, ""],
            "value": [10, "20", "30", "not"],
            "year": ["2020", "1999", None, "no"],
        }
    )
    analyze_table_quality(demo, table_name="demo")
