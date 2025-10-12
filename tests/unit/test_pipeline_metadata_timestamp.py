"""Tests for :func:`library.pipelines.common.metadata.get_timestamp_utc`."""

from __future__ import annotations

from datetime import UTC
from datetime import datetime as dt
from types import SimpleNamespace

import pytest

from library.common.run_context import RunContext, set_current
from library.pipelines.common import metadata as pipeline_metadata


@pytest.mark.unit
def test_get_timestamp_utc__uses_current_time_when_generated_at_blank(monkeypatch):
    """Fallback to the current UTC timestamp when the context has no value."""

    pipeline_metadata.get_timestamp_utc.cache_clear()
    set_current(RunContext(run_id="test", generated_at=""))

    fixed_now = dt(2024, 1, 2, 3, 4, 5, tzinfo=UTC)
    monkeypatch.setattr(
        pipeline_metadata,
        "datetime",
        SimpleNamespace(now=lambda tz=None: fixed_now),
    )

    try:
        timestamp = pipeline_metadata.get_timestamp_utc()
        assert timestamp == fixed_now.isoformat()
    finally:
        set_current(None)
        pipeline_metadata.get_timestamp_utc.cache_clear()
