from __future__ import annotations

import pandas as pd
import pytest

from library import chembl_library as cl
from library.chembl_client import ChemblClient
from library.dtype_inspector import inspect_dtypes, logger


def _make_df(name: str) -> pd.DataFrame:
    """Return a single-row data frame with deterministic dtypes."""

    return pd.DataFrame({f"{name}_id": pd.Series([1], dtype="int64")})


def test_inspect_dtypes_uses_mocked_calls(monkeypatch: pytest.MonkeyPatch) -> None:
    """Inspector returns dtypes and avoids network access."""

    monkeypatch.setattr(cl, "get_activities", lambda *a, **k: _make_df("activity"))
    monkeypatch.setattr(cl, "get_assays", lambda *a, **k: _make_df("assay"))
    monkeypatch.setattr(cl, "get_documents", lambda *a, **k: _make_df("document"))
    monkeypatch.setattr(cl, "get_targets", lambda *a, **k: _make_df("target"))
    monkeypatch.setattr(cl, "get_testitem", lambda *a, **k: _make_df("testitem"))

    def fail_request_json(
        self, url: str, *, cfg, timeout=None
    ) -> None:  # pragma: no cover - safety
        raise AssertionError("network access not allowed")

    monkeypatch.setattr(ChemblClient, "request_json", fail_request_json)

    log_calls: list[tuple[str, dict[str, object]]] = []
    monkeypatch.setattr(
        logger,
        "info",
        lambda event, **kw: log_calls.append((event, kw)),
    )

    result = inspect_dtypes()

    assert result["activities"]["activity_id"] == "int64"
    assert any(
        ev == "dtypes" and data.get("dataset") == "activities" for ev, data in log_calls
    )
