"""Shared type hints and dataclasses for postprocessing steps."""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field
from types import MappingProxyType
from typing import Any, Protocol

import pandas as pd


class StepFn(Protocol):
    """Protocol describing a pure DataFrame transformation."""

    def __call__(self, df: pd.DataFrame, **kwargs: Any) -> pd.DataFrame:
        """Return a new :class:`pandas.DataFrame` derived from ``df``."""


@dataclass(frozen=True)
class StepDefinition:
    """Metadata describing a transformation step."""

    name: str
    func: StepFn
    description: str | None = None
    params: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:  # pragma: no cover - dataclass validation
        if not callable(self.func):
            raise TypeError("StepDefinition.func must be callable")
        params = dict(self.params)
        object.__setattr__(self, "params", MappingProxyType(params))


class StepError(RuntimeError):
    """Base error raised when a step fails."""

    def __init__(
        self, step_name: str, message: str, *, cause: BaseException | None = None
    ):
        self.step_name = step_name
        self.cause = cause
        error_message = f"Step '{step_name}' failed: {message}"
        super().__init__(error_message)


class SchemaValidationError(StepError):
    """Raised when a DataFrame does not comply with the expected schema."""


class ImportResolutionError(RuntimeError):
    """Raised when a dotted import path cannot be resolved."""


StepIterable = Iterable[StepDefinition]


__all__ = [
    "ImportResolutionError",
    "SchemaValidationError",
    "StepDefinition",
    "StepError",
    "StepFn",
    "StepIterable",
]
