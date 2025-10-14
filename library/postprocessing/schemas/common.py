from __future__ import annotations

import logging
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from types import MappingProxyType
from typing import Any

import pandas as pd
from pandera.errors import SchemaError, SchemaErrors

from library._compat.pandera import pa
from library.postprocessing.common.logging import get_logger
from library.postprocessing.common.types import SchemaValidationError


@dataclass(frozen=True, slots=True)
class TableSchema:
    """Immutable schema declaration backed by :mod:`pandera`."""

    name: str | None = None
    required_columns: Sequence[str] = field(default_factory=tuple)
    optional_columns: Sequence[str] = field(default_factory=tuple)
    dtypes: Mapping[str, Any] = field(default_factory=dict)
    nullable_columns: Sequence[str] | None = None
    column_order: Sequence[str] | None = None
    sort_by: Sequence[str] | None = None
    coerce: bool = True
    _nullable_set: frozenset[str] = field(init=False, repr=False)
    _pandera_schema: pa.DataFrameSchema = field(init=False, repr=False)

    def __post_init__(self) -> None:
        required = tuple(self.required_columns)
        optional = tuple(self.optional_columns)
        nullable = set(optional)
        if self.nullable_columns is not None:
            nullable.update(self.nullable_columns)
        object.__setattr__(self, "required_columns", required)
        object.__setattr__(self, "optional_columns", optional)
        object.__setattr__(self, "nullable_columns", tuple(sorted(nullable)))
        dtypes = MappingProxyType(dict(self.dtypes))
        object.__setattr__(self, "dtypes", dtypes)
        column_order = tuple(self.column_order) if self.column_order else None
        sort_by = tuple(self.sort_by) if self.sort_by else None
        object.__setattr__(self, "column_order", column_order)
        object.__setattr__(self, "sort_by", sort_by)
        object.__setattr__(self, "_nullable_set", frozenset(nullable))
        object.__setattr__(self, "_pandera_schema", self._build_pandera_schema())

    def all_columns(self) -> Sequence[str]:
        """Return the ordered union of required and optional columns."""

        ordered = list(self.required_columns) + list(self.optional_columns)
        if self.column_order:
            seen: set[str] = set()
            custom_order: list[str] = []
            for column in self.column_order:
                if column not in seen and column in ordered:
                    custom_order.append(column)
                    seen.add(column)
            remaining = [
                column for column in ordered if column not in seen
            ]
            return tuple(custom_order + remaining)
        return tuple(ordered)

    def pandera_schema(self) -> pa.DataFrameSchema:
        """Return the underlying :class:`pandera.DataFrameSchema`."""

        return self._pandera_schema

    def validate(
        self,
        df: pd.DataFrame,
        *,
        context: str,
        logger: logging.Logger | None = None,
    ) -> pd.DataFrame:
        """Validate ``df`` and raise :class:`SchemaValidationError` on failure."""

        return validate_with_pandera(df, self, context=context, logger=logger)

    def _build_pandera_schema(self) -> pa.DataFrameSchema:
        columns: dict[str, pa.Column] = {}
        for column in (*self.required_columns, *self.optional_columns):
            dtype = self.dtypes.get(column, object)
            columns[column] = pa.Column(
                dtype,
                required=column in self.required_columns,
                nullable=column in self._nullable_set,
                coerce=self.coerce,
            )
        return pa.DataFrameSchema(columns, coerce=self.coerce)


def _log_schema_failures(
    *,
    schema: TableSchema,
    context: str,
    error: SchemaErrors,
    logger: logging.Logger,
) -> None:
    failure_cases = error.failure_cases
    unique_columns = sorted(
        {str(column) for column in failure_cases.get("column", pd.Series()).dropna()}
    )
    logger.error(
        "Schema validation failed | schema=%s | context=%s | errors=%d | columns=%s",
        schema.name or "DataFrameSchema",
        context,
        len(failure_cases.index),
        ",".join(unique_columns) or "<unknown>",
    )
    preview = failure_cases.head(10)
    for row in preview.itertuples(index=False):
        logger.error(
            "Validation failure | column=%s | check=%s | failure_case=%s | index=%s",
            getattr(row, "column", "<unknown>"),
            getattr(row, "check", "<unknown>"),
            getattr(row, "failure_case", "<unknown>"),
            getattr(row, "index", "<unknown>"),
        )


def validate_with_pandera(
    df: pd.DataFrame,
    schema: TableSchema,
    *,
    context: str,
    logger: logging.Logger | None = None,
) -> pd.DataFrame:
    """Validate ``df`` against ``schema`` using Pandera with detailed logging."""

    log = logger or get_logger(f"chembl.postprocess.schema.{schema.name or 'anonymous'}")
    pandera_schema = schema.pandera_schema()
    try:
        validated = pandera_schema.validate(df, lazy=True)
    except SchemaErrors as exc:
        _log_schema_failures(schema=schema, context=context, error=exc, logger=log)
        message = f"{len(exc.failure_cases.index)} schema violations"
        raise SchemaValidationError(context, message, cause=exc) from exc
    except SchemaError as exc:
        log.error(
            "Schema validation error | schema=%s | context=%s | error=%s",
            schema.name or "DataFrameSchema",
            context,
            exc,
        )
        raise SchemaValidationError(context, str(exc), cause=exc) from exc
    return validated


__all__ = ["TableSchema", "validate_with_pandera"]
