from __future__ import annotations

import pandas as pd
import pytest

from scripts import get_activity_data


@pytest.mark.unit
def test_filter_activity_output_columns__drops_expected_columns(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    frame = pd.DataFrame(
        [
            {
                "activity_id": "CHEMBL1",
                "approx_cited_activity": True,
                "value": 1.23,
                "standard_upper_value": 2.34,
                "shuffled_cit": False,
            }
        ]
    )

    messages: list[str] = []

    def _capture_info(message: str, *args: object, **kwargs: object) -> None:
        if args:
            try:
                formatted = message % args
            except TypeError:
                formatted = message
        else:
            formatted = message
        messages.append(formatted)

    monkeypatch.setattr(get_activity_data.logger, "info", _capture_info)

    result = get_activity_data._filter_activity_output_columns(frame)

    assert list(result.columns) == ["activity_id", "value"]
    assert result.iloc[0].to_dict() == {"activity_id": "CHEMBL1", "value": 1.23}
    assert (
        "Dropped columns from output.activity_*: approx_cited_activity, standard_upper_value, shuffled_cit"
        in messages
    )
