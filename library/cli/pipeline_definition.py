"""Data structures describing reusable pipeline configuration."""

from __future__ import annotations

from collections.abc import Callable, Iterable, Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Protocol, TypeVar

import pandas as pd

from ..common.metadata import Stats

SchemaT = TypeVar("SchemaT")


class ValidationResult(Protocol):
    """Protocol describing the return type of validator callables."""

    data: pd.DataFrame
    failure_cases: pd.DataFrame


class Validator(Protocol):
    """Callable interface for data frame validators."""

    def __call__(
        self, df: pd.DataFrame
    ) -> ValidationResult:  # pragma: no cover - Protocol
        ...


MetadataHook = Callable[[pd.DataFrame], pd.DataFrame]
Writer = Callable[
    [Iterable[pd.DataFrame], Path, Sequence[str] | None, Sequence[str]], Path
]
TableQualityHook = Callable[[Path], None]
Fetcher = Callable[[], Iterable[pd.DataFrame] | pd.DataFrame]


@dataclass(slots=True)
class PipelineDefinition:
    """Configuration bundle describing how :func:`run_pipeline` should behave."""

    schema: SchemaT | None
    schema_name: str
    writer: Writer
    validators: Sequence[Validator] = field(default_factory=tuple)
    metadata_hooks: Sequence[MetadataHook] = field(default_factory=tuple)
    command: str | None = None
    invocation: Sequence[str] | None = None
    config_snapshot: Mapping[str, object] = field(default_factory=dict)
    inputs: Mapping[str, object] = field(default_factory=dict)
    key_columns: Sequence[str] = field(default_factory=tuple)
    table_quality: TableQualityHook | None = None
    stats_extra: Mapping[str, object] | Callable[[], Mapping[str, object]] | None = None
    #: Additional statistics merged into the metadata output produced by
    #: :func:`library.cli_utils.run_pipeline`. The mapping must only contain
    #: JSON-serialisable values and may be provided directly or via a factory.
    stats_callback: Callable[[Stats], None] | None = None
    strict_mode: bool = False
    dictionary_resources: Sequence[str] | None = None

    def __post_init__(self) -> None:  # pragma: no cover - simple normalisation
        self.validators = tuple(self.validators)
        self.metadata_hooks = tuple(self.metadata_hooks)
        if self.invocation is not None:
            self.invocation = tuple(self.invocation)
        self.key_columns = tuple(self.key_columns)
        if self.dictionary_resources is not None:
            self.dictionary_resources = tuple(self.dictionary_resources)

    @classmethod
    def from_legacy_kwargs(
        cls,
        *,
        schema: SchemaT | None,
        schema_name: str,
        writer: Writer,
        validators: Sequence[Validator] | None = None,
        metadata_hooks: Sequence[MetadataHook] | None = None,
        command: str | None = None,
        invocation: Sequence[str] | None = None,
        config_snapshot: Mapping[str, object] | None = None,
        inputs: Mapping[str, object] | None = None,
        key_columns: Sequence[str] | None = None,
        table_quality: TableQualityHook | None = None,
        stats_extra: (
            Mapping[str, object] | Callable[[], Mapping[str, object]] | None
        ) = None,
        stats_callback: Callable[[Stats], None] | None = None,
        strict_mode: bool = False,
        dictionary_resources: Sequence[str] | None = None,
    ) -> PipelineDefinition:
        """Construct a :class:`PipelineDefinition` using legacy keyword arguments."""

        if table_quality is None:
            raise TypeError(
                "'table_quality' must be provided when using legacy arguments"
            )

        return cls(
            schema=schema,
            schema_name=schema_name,
            writer=writer,
            validators=validators or (),
            metadata_hooks=metadata_hooks or (),
            command=command,
            invocation=invocation,
            config_snapshot=dict(config_snapshot or {}),
            inputs=dict(inputs or {}),
            key_columns=tuple(key_columns or ()),
            table_quality=table_quality,
            stats_extra=stats_extra,
            stats_callback=stats_callback,
            strict_mode=strict_mode,
            dictionary_resources=dictionary_resources,
        )


LEGACY_DEFINITION_KEYS = {
    "schema",
    "schema_name",
    "writer",
    "validators",
    "metadata_hooks",
    "command",
    "invocation",
    "config_snapshot",
    "inputs",
    "key_columns",
    "table_quality",
    "stats_extra",
    "stats_callback",
    "strict_mode",
    "dictionary_resources",
}


def normalise_definition(
    definition: PipelineDefinition | None,
    extra_kwargs: Mapping[str, object],
) -> PipelineDefinition:
    """Return a :class:`PipelineDefinition` using legacy keyword fallbacks."""

    if definition is not None:
        if extra_kwargs:
            unexpected = ", ".join(sorted(extra_kwargs))
            raise TypeError(
                "run_pipeline received unexpected legacy keyword arguments: "
                f"{unexpected}"
            )
        return definition

    remaining = dict(extra_kwargs)
    legacy_args: dict[str, object] = {}
    for key in LEGACY_DEFINITION_KEYS:
        if key in remaining:
            legacy_args[key] = remaining.pop(key)

    if remaining:
        unexpected = ", ".join(sorted(remaining))
        raise TypeError(
            "run_pipeline received unexpected keyword arguments: " f"{unexpected}"
        )

    try:
        return PipelineDefinition.from_legacy_kwargs(**legacy_args)
    except TypeError as exc:  # pragma: no cover - mirrors legacy errors
        raise TypeError(
            "run_pipeline requires a PipelineDefinition or the legacy keyword "
            "arguments (schema, schema_name, writer, table_quality, ...)."
        ) from exc
