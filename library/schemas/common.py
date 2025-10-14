"""Factory helpers for reusable Pandera column declarations."""

from __future__ import annotations

from typing import Any

import pandas as pd

from library._compat.pandera import pa

__all__ = [
    "boolean_column",
    "float_column",
    "int_column",
    "object_column",
    "string_column",
]


def string_column(
    *, required: bool = False, nullable: bool = True, coerce: bool = True
) -> pa.Column:
    """Return a :class:`pandera.Column` configured for string data."""

    dtype = pd.StringDtype()
    return pa.Column(dtype, required=required, nullable=nullable, coerce=coerce)


def int_column(
    *, required: bool = False, nullable: bool = True, coerce: bool = True
) -> pa.Column:
    """Return a nullable integer column accepting standard pandas integers."""

    return pa.Column("Int64", required=required, nullable=nullable, coerce=coerce)


def float_column(
    *, required: bool = False, nullable: bool = True, coerce: bool = True
) -> pa.Column:
    """Return a numeric column with optional coercion to :class:`float`."""

    return pa.Column(float, required=required, nullable=nullable, coerce=coerce)


def boolean_column(
    *, required: bool = False, nullable: bool = True, coerce: bool = True
) -> pa.Column:
    """Return a boolean column backed by :class:`pandas.BooleanDtype`."""

    dtype = pd.BooleanDtype()
    return pa.Column(dtype, required=required, nullable=nullable, coerce=coerce)


def object_column(
    dtype: Any | None = object,
    *,
    required: bool = False,
    nullable: bool = True,
    coerce: bool = True,
) -> pa.Column:
    """Return a column accepting arbitrary Python objects."""

    return pa.Column(dtype, required=required, nullable=nullable, coerce=coerce)
