from __future__ import annotations

import pandas as pd
import pytest

from library.config import ActivityBoundsCfg
from schemas import normalize_activities
from scripts.get_activity_data import compute_activity_bounds, logger as activity_logger


def _compute(data: list[dict[str, object]], cfg: ActivityBoundsCfg | None = None) -> pd.DataFrame:
    frame = pd.DataFrame(data)
    normalised = normalize_activities(frame)
    return compute_activity_bounds(normalised, cfg or ActivityBoundsCfg())


def test_bounds_equal_relation() -> None:
    result = _compute([
        {"standard_value": 5.0, "standard_relation": "="},
    ])
    row = result.iloc[0]
    assert row["lower_value"] == pytest.approx(5.0)
    assert row["upper_value"] == pytest.approx(5.0)


def test_bounds_greater_relation() -> None:
    result = _compute([
        {"standard_value": 3.0, "standard_relation": ">"},
    ])
    row = result.iloc[0]
    assert row["lower_value"] == pytest.approx(3.0)
    assert pd.isna(row["upper_value"])


def test_bounds_less_relation() -> None:
    result = _compute([
        {"standard_value": 4.0, "standard_relation": "<"},
    ])
    row = result.iloc[0]
    assert pd.isna(row["lower_value"])
    assert row["upper_value"] == pytest.approx(4.0)


def test_bounds_range_pair() -> None:
    result = _compute([
        {"standard_value": 1.0, "standard_upper_value": 2.0},
    ])
    row = result.iloc[0]
    assert row["lower_value"] == pytest.approx(1.0)
    assert row["upper_value"] == pytest.approx(2.0)


def test_bounds_unknown_relation_logs(monkeypatch: pytest.MonkeyPatch) -> None:
    events: list[str] = []

    def fake_warning(event: str, **_: object) -> None:
        events.append(event)

    monkeypatch.setattr(activity_logger, "warning", fake_warning)
    result = _compute([
        {"standard_value": 7.0, "standard_relation": "??"},
    ])
    assert "activity_bounds_unknown_relation" in events
    assert result["lower_value"].isna().all()
    assert result["upper_value"].isna().all()
