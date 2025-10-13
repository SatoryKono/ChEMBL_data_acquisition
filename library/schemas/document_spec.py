"""Utilities for loading the declarative document schema."""

from __future__ import annotations

from collections import OrderedDict
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from types import MappingProxyType
from typing import Any

from library._compat.pandera import pa
import yaml

from config.paths import SCHEMA_DIR


@dataclass(frozen=True)
class ColumnSpec:
    """Declarative attributes describing a single column."""

    name: str
    dtype: str
    required: bool
    nullable: bool
    coerce: bool
    description: str | None

    def to_pandera_column(self) -> pa.Column:
        """Return a :class:`pandera.Column` for the specification."""

        dtype = _resolve_dtype(self.dtype)
        kwargs: dict[str, Any] = {"required": self.required, "nullable": self.nullable}
        if self.coerce:
            kwargs["coerce"] = True
        return pa.Column(dtype, **kwargs)


@dataclass(frozen=True)
class ColumnGroup:
    """Logical grouping of related columns."""

    name: str
    title: str | None
    description: str | None
    columns: tuple[ColumnSpec, ...]


@dataclass(frozen=True)
class DocumentDeclaration:
    """Fully materialised declaration of the document schema."""

    groups: tuple[ColumnGroup, ...]
    column_specs: Mapping[str, ColumnSpec]
    ordered_columns: tuple[str, ...]
    schema: pa.DataFrameSchema
    export_columns: tuple[str, ...]


def load_document_declaration(path: str | Path | None = None) -> DocumentDeclaration:
    """Load the document schema declaration from ``path``."""

    schema_path = _resolve_schema_path(path)
    with schema_path.open(encoding="utf-8") as handle:
        raw_data = yaml.safe_load(handle) or {}

    groups_data = raw_data.get("groups", [])
    if not isinstance(groups_data, Sequence):
        raise TypeError("document schema 'groups' must be a sequence")

    groups: list[ColumnGroup] = []
    ordered_specs: list[ColumnSpec] = []
    specs_by_name: dict[str, ColumnSpec] = {}

    for group_data in groups_data:
        if not isinstance(group_data, Mapping):
            raise TypeError("each group entry must be a mapping")
        group_name = str(group_data.get("name"))
        if not group_name:
            raise ValueError("group entries must define a non-empty 'name'")

        column_entries = group_data.get("columns", [])
        if not isinstance(column_entries, Sequence):
            raise TypeError(f"group '{group_name}' columns must be a sequence")

        group_columns: list[ColumnSpec] = []
        for column_data in column_entries:
            if not isinstance(column_data, Mapping):
                raise TypeError(
                    f"column declaration in group '{group_name}' must be a mapping"
                )
            column_name = str(column_data.get("name"))
            if not column_name:
                raise ValueError(f"column in group '{group_name}' is missing a 'name'")
            if column_name in specs_by_name:
                raise ValueError(f"duplicate column declaration for '{column_name}'")

            spec = ColumnSpec(
                name=column_name,
                dtype=str(column_data.get("dtype", "object")),
                required=bool(column_data.get("required", False)),
                nullable=bool(column_data.get("nullable", True)),
                coerce=bool(column_data.get("coerce", False)),
                description=column_data.get("description"),
            )
            specs_by_name[column_name] = spec
            ordered_specs.append(spec)
            group_columns.append(spec)

        groups.append(
            ColumnGroup(
                name=group_name,
                title=group_data.get("title"),
                description=group_data.get("description"),
                columns=tuple(group_columns),
            )
        )

    schema_columns = OrderedDict(
        (spec.name, spec.to_pandera_column()) for spec in ordered_specs
    )
    schema = pa.DataFrameSchema(schema_columns, ordered=True)

    export_columns = _parse_export_columns(raw_data.get("export", {}))

    return DocumentDeclaration(
        groups=tuple(groups),
        column_specs=MappingProxyType(dict(specs_by_name)),
        ordered_columns=tuple(spec.name for spec in ordered_specs),
        schema=schema,
        export_columns=export_columns,
    )


def _resolve_schema_path(path: str | Path | None) -> Path:
    if path is None:
        return SCHEMA_DIR / "document.yaml"
    resolved = Path(path)
    if resolved.is_dir():
        resolved = resolved / "document.yaml"
    return resolved


def _resolve_dtype(value: str) -> Any:
    key = value.lower().strip()
    if key in {"string", "str"}:
        return str
    if key in {"int", "int64"}:
        return int
    if key in {"float", "float64"}:
        return float
    if key in {"object", "any"}:
        return object
    raise ValueError(f"unsupported column dtype '{value}' in document schema")


def _parse_export_columns(export_data: Mapping[str, Any]) -> tuple[str, ...]:
    if not export_data:
        return tuple()
    columns = export_data.get("columns", [])
    if not isinstance(columns, Sequence):
        raise TypeError("export.columns must be a sequence")
    resolved: list[str] = []
    for entry in columns:
        if isinstance(entry, Mapping):
            name = entry.get("name")
        else:
            name = entry
        if not name:
            raise ValueError("export column entries must define a name")
        resolved.append(str(name))
    return tuple(resolved)


_DECLARATION = load_document_declaration()

DOCUMENT_DECLARATION: DocumentDeclaration = _DECLARATION
"""Materialised declaration of the document schema."""

DOCUMENT_SCHEMA_COLUMNS: list[str] = list(_DECLARATION.ordered_columns)
"""Canonical schema column ordering for document metadata."""

DOCUMENT_COLUMN_GROUPS: Mapping[str, tuple[str, ...]] = MappingProxyType(
    {
        group.name: tuple(spec.name for spec in group.columns)
        for group in _DECLARATION.groups
    }
)
"""Mapping of group name to the ordered column names it contains."""

DOCUMENT_COLUMN_SPECS: Mapping[str, ColumnSpec] = _DECLARATION.column_specs
"""Mapping of column names to their declarative specification."""

DOCUMENT_EXPORT_COLUMNS: tuple[str, ...] = _DECLARATION.export_columns
"""Default export projection used by the document CLI."""

__all__ = [
    "ColumnGroup",
    "ColumnSpec",
    "DocumentDeclaration",
    "DOCUMENT_DECLARATION",
    "DOCUMENT_SCHEMA_COLUMNS",
    "DOCUMENT_COLUMN_GROUPS",
    "DOCUMENT_COLUMN_SPECS",
    "DOCUMENT_EXPORT_COLUMNS",
    "load_document_declaration",
]
