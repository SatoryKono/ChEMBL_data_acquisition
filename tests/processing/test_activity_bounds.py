from __future__ import annotations

import pandas as pd
import pytest

from library.config import ActivityBoundsCfg
from library.log import logger as activity_logger
from library.processing.activity import compute_activity_bounds
from library.schemas import normalize_activities


def _compute_frame(
    data: list[dict[str, object]], cfg: ActivityBoundsCfg | None = None
) -> pd.DataFrame:
    frame = pd.DataFrame(data)
    normalised = normalize_activities(frame)
    return compute_activity_bounds(normalised, cfg or ActivityBoundsCfg())


def test_compute_activity_bounds_unit_normalization() -> None:
    df = _compute_frame(
        [
            {
                "standard_value": 1.0,
                "standard_relation": "=",
                "value": 1000.0,
                "units": "pM",
            },
            {
                "standard_value": 1.0,
                "standard_relation": "=",
                "value": 0.001,
                "units": "uM",
            },
        ]
    )
    assert df["lower_value"].tolist() == pytest.approx([1.0, 1.0])
    assert df["upper_value"].tolist() == pytest.approx([1.0, 1.0])


def test_compute_activity_bounds_swaps_conflicts(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    events: list[str] = []

    def fake_warning(event: str, **_: object) -> None:
        events.append(event)

    monkeypatch.setattr(activity_logger, "warning", fake_warning)
    df = _compute_frame(
        [
            {"standard_lower_value": 5.0, "standard_upper_value": 2.0},
        ]
    )
    assert df.loc[0, "lower_value"] == pytest.approx(2.0)
    assert df.loc[0, "upper_value"] == pytest.approx(5.0)
    assert "activity_bounds_conflict_swapped" in events


def test_compute_activity_bounds_round_and_clamp() -> None:
    cfg = ActivityBoundsCfg()
    cfg.rounding_digits = 3
    df = _compute_frame(
        [
            {
                "standard_lower_value": -0.00049,
                "standard_upper_value": 0.123456,
                "standard_type": "IC50",
                "standard_units": "nM",
            }
        ],
        cfg=cfg,
    )
    assert df.loc[0, "lower_value"] == 0.0
    assert df.loc[0, "upper_value"] == pytest.approx(0.123)


def test_compute_activity_bounds_no_clamp_when_disabled() -> None:
    cfg = ActivityBoundsCfg(clamp_nonnegative=False)
    df = _compute_frame(
        [
            {
                "standard_lower_value": -1.0,
                "standard_upper_value": 1.0,
                "standard_type": "IC50",
                "standard_units": "nM",
            }
        ],
        cfg=cfg,
    )
    assert df.loc[0, "lower_value"] == pytest.approx(-1.0)


def test_compute_activity_bounds_uncertainty_toggle() -> None:
    row = {
        "standard_value": 8.0,
        "standard_text_value": "8 ± 0.25",
        "standard_type": "Ki",
        "standard_units": "nM",
    }
    disabled = _compute_frame([row])
    assert disabled["lower_value"].isna().all()
    assert disabled["upper_value"].isna().all()

    cfg = ActivityBoundsCfg(enable_from_uncertainty=True)
    enabled = _compute_frame([row], cfg=cfg)
    assert enabled.loc[0, "lower_value"] == pytest.approx(7.75)
    assert enabled.loc[0, "upper_value"] == pytest.approx(8.25)


def test_compute_activity_bounds_golden_samples() -> None:
    input_df = pd.read_csv("tests/data/activity_bounds_input.csv")
    expected = pd.read_csv("tests/data/activity_bounds_expected.csv")
    result = compute_activity_bounds(
        normalize_activities(input_df), ActivityBoundsCfg()
    )
    pd.testing.assert_frame_equal(
        result[["activity_id", "lower_value", "upper_value"]],
        expected,
        check_like=True,
    )


def test_compute_activity_bounds_preserves_other_columns() -> None:
    df = _compute_frame(
        [
            {"standard_value": 5.0, "standard_relation": "=", "extra": "keep"},
        ]
    )
    assert df.loc[0, "extra"] == "keep"


def test_compute_activity_bounds_disable_relation() -> None:
    cfg = ActivityBoundsCfg(enable_from_relation=False)
    df = _compute_frame(
        [
            {"standard_value": 4.2, "standard_relation": "="},
        ],
        cfg=cfg,
    )
    assert df["lower_value"].isna().all()
    assert df["upper_value"].isna().all()
