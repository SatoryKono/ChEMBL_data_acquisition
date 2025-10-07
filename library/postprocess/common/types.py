"""Shared type hints and dataclasses for postprocessing steps."""
from __future__ import annotations

from dataclasses import dataclass
from types import MappingProxyType
from typing import Any, Iterable, Mapping, Optional, Protocol, Tuple, Union

import pandas as pd


class StepFn(Protocol):
    """Protocol describing a pure DataFrame transformation."""

    def __call__(self, df: pd.DataFrame, **params: Any) -> pd.DataFrame:  # pragma: no cover - structural
        """Return a new :class:`pandas.DataFrame` derived from ``df``."""


@dataclass(frozen=True)
class StepDefinition:
    """Metadata describing a transformation step."""

    name: str
    func: StepFn
    description: Optional[str] = None
    parameters: Mapping[str, Any] | None = None

    def __post_init__(self) -> None:  # pragma: no cover - dataclass validation
        if not callable(self.func):
            raise TypeError("StepDefinition.func must be callable")
        if self.parameters is None:
            object.__setattr__(self, "parameters", MappingProxyType({}))
        elif isinstance(self.parameters, Mapping):
            object.__setattr__(
                self,
                "parameters",
                MappingProxyType(dict(self.parameters)),
            )
        else:
            raise TypeError("StepDefinition.parameters must be a mapping if provided")


class StepError(RuntimeError):
    """Base error raised when a step fails."""

    def __init__(self, step_name: str, message: str, *, cause: Optional[BaseException] = None):
        self.step_name = step_name
        self.cause = cause
        error_message = f"Step '{step_name}' failed: {message}"
        super().__init__(error_message)


class SchemaValidationError(StepError):
    """Raised when a DataFrame does not comply with the expected schema."""


class ImportResolutionError(RuntimeError):
    """Raised when a dotted import path cannot be resolved."""


StepLike = Union[
    StepDefinition,
    StepFn,
    Tuple[StepFn],
    Tuple[StepFn, Mapping[str, Any] | None],
    Tuple[str, StepFn],
    Tuple[str, StepFn, Mapping[str, Any] | None],
]
StepIterable = Iterable[StepLike]


__all__ = [
    "ImportResolutionError",
    "SchemaValidationError",
    "StepDefinition",
    "StepError",
    "StepFn",
    "StepIterable",
]
