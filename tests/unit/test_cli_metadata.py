"""Unit tests for :mod:`library.cli.metadata`."""

from __future__ import annotations

import pytest

from library.cli.metadata import DEFAULT_OPTION_SOURCE, prepare_option
from library.config import ConfigMetadata, ConfigSource


@pytest.mark.unit
def test_prepare_option__metadata_absent_uses_defaults() -> None:
    """Return placeholder metadata when the config snapshot is unavailable."""

    result = prepare_option(None)

    assert result == {"value": None, "source": DEFAULT_OPTION_SOURCE}


@pytest.mark.unit
def test_prepare_option__metadata_absent_applies_value_and_detail() -> None:
    """Include override information when no metadata object is provided."""

    result = prepare_option(
        None,
        value="example",
        default_source="cli",
        default_detail="runtime",
    )

    assert result == {"value": "example", "source": "cli", "detail": "runtime"}


@pytest.mark.unit
def test_prepare_option__metadata_present_respects_snapshot() -> None:
    """Reuse existing metadata entries and override the value when requested."""

    metadata = ConfigMetadata(
        snapshot={
            "assay": {
                "limit": {"value": 50, "source": "config", "detail": "config.yaml"}
            }
        },
        sources={("assay", "limit"): ConfigSource("config", "config.yaml")},
        cli_paths={"limit": ("assay", "limit")},
    )

    default_result = prepare_option(metadata, argument="limit")
    override_result = prepare_option(metadata, argument="limit", value=100)

    assert default_result == {"value": 50, "source": "config", "detail": "config.yaml"}
    assert override_result == {
        "value": 100,
        "source": "config",
        "detail": "config.yaml",
    }
