import pandas as pd

from library.config import ActivityBoundsCfg
from library.processing.activity.bounds import compute_activity_bounds


def test_compute_activity_bounds__preserves_standard_ranges() -> None:
    rows: list[dict[str, object]] = []
    for idx in range(20):
        rows.append(
            {
                "activity_id": f"ACT{idx}",
                "standard_lower_value": float(idx),
                "standard_upper_value": float(idx) + 1.0,
                "standard_value": float(idx) + 0.5,
                "standard_type": "IC50",
                "standard_units": "nM",
            }
        )

    rows[0]["standard_lower_value"] = None
    rows[0]["standard_upper_value"] = None

    frame = pd.DataFrame(rows)
    cfg = ActivityBoundsCfg()

    result = compute_activity_bounds(frame, cfg)

    lower_ratio = result["lower_value"].notna().mean()
    upper_ratio = result["upper_value"].notna().mean()

    assert lower_ratio >= 0.95
    assert upper_ratio >= 0.95

    populated = result.dropna(subset=["standard_lower_value", "standard_upper_value"])
    pd.testing.assert_series_equal(
        populated["lower_value"].reset_index(drop=True),
        pd.to_numeric(populated["standard_lower_value"], errors="coerce").reset_index(
            drop=True
        ),
        check_names=False,
        check_dtype=False,
    )
    pd.testing.assert_series_equal(
        populated["upper_value"].reset_index(drop=True),
        pd.to_numeric(populated["standard_upper_value"], errors="coerce").reset_index(
            drop=True
        ),
        check_names=False,
        check_dtype=False,
    )
