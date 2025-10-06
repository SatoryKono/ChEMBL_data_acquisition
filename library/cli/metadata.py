"""Helper utilities shared across CLI metadata consumers."""

from __future__ import annotations

from typing import Any

from library.config import ConfigMetadata

DEFAULT_OPTION_SOURCE = "unknown"
"""Default provenance value used when metadata is unavailable."""

OPTION_UNSET = object()
"""Sentinel value identifying parameters without an explicit override."""


def prepare_option(
    metadata: ConfigMetadata | None,
    *,
    argument: str | None = None,
    path: str | None = None,
    value: Any = OPTION_UNSET,
    default_source: str = DEFAULT_OPTION_SOURCE,
    default_detail: str | None = None,
) -> dict[str, Any]:
    """Return metadata entry for an option, handling missing configuration."""

    if metadata is not None:
        if value is OPTION_UNSET:
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

    actual = None if value is OPTION_UNSET else value
    entry: dict[str, Any] = {"value": actual, "source": default_source}
    if default_detail is not None:
        entry["detail"] = default_detail
    return entry


__all__ = ["DEFAULT_OPTION_SOURCE", "OPTION_UNSET", "prepare_option"]
