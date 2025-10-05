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


def test_augment_activity_frame__creates_log_value_from_standard_value():
    frame = pd.DataFrame(
        {
            "standard_value": [1.0, 10.0, pd.NA],
            "standard_units": ["nM", "µM", "nM"],
        }
    )

    result, filled = _augment_activity_frame(frame)

    expected = pd.Series([9.0, 5.0, pd.NA], dtype="Float64", name="log_value")

    assert str(result["log_value"].dtype) == "Float64"
    pd.testing.assert_series_equal(result["log_value"], expected)
    assert "log_value" in filled
