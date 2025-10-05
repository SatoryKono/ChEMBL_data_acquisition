import pandas as pd

from library.postprocessing.activity_extended import _augment_activity_frame


def test_augment_activity_frame__fills_missing_log_value():
    frame = pd.DataFrame(
        {
            "log_value": pd.Series(["not-a-number", pd.NA], dtype="string"),
            "standard_value": [1.0, 1.0],
            "standard_units": ["nM", "nM"],
        }
    )

    result, _ = _augment_activity_frame(frame)

    expected = pd.Series([9.0, 9.0], dtype="Float64", name="log_value")

    assert str(result["log_value"].dtype) == "Float64"
    pd.testing.assert_series_equal(result["log_value"], expected)
