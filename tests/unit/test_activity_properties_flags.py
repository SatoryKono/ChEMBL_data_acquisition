import pandas as pd
import pytest

from library.postprocessing.activity_extended import _extract_activity_properties_flags


@pytest.mark.unit
@pytest.mark.pipeline_scenario("transformation_rules")
def test_extract_activity_properties_flags__parses_standard_payload() -> None:
    frame = pd.DataFrame(
        {
            "activity_properties": [
                '{"effect_features": {"allosteric": true, "positive": true, "negative": false}}',
                '{"effect_features": {"allosteric": false, "positive": false, "negative": true}}',
                '{"effect_features": {"allosteric": false, "positive": false, "negative": false}}',
            ]
        }
    )

    result = _extract_activity_properties_flags(frame)

    assert result["allosteric"].tolist() == [True, False, False]
    assert result["pam"].tolist() == [True, False, False]
    assert result["nam"].tolist() == [False, True, False]
    for column in ("allosteric", "pam", "nam"):
        assert result[column].dtype == pd.BooleanDtype()


@pytest.mark.unit
def test_extract_activity_properties_flags__handles_split_and_replace_payload() -> None:
    raw = '"{""effect_features"": {""allosteric"": ""true"", ""positive"": ""1"", ""negative"": ""0""}}"'
    frame = pd.DataFrame({"activity_properties": [raw]})

    result = _extract_activity_properties_flags(frame)

    assert bool(result.loc[0, "allosteric"]) is True
    assert bool(result.loc[0, "pam"]) is True
    assert bool(result.loc[0, "nam"]) is False
    assert result["allosteric"].dtype == pd.BooleanDtype()


@pytest.mark.unit
def test_extract_activity_properties_flags__malformed_payload_results_in_na() -> None:
    frame = pd.DataFrame({"activity_properties": ["not-json"]})

    result = _extract_activity_properties_flags(frame)

    assert result.loc[0, "allosteric"] is pd.NA
    assert result.loc[0, "pam"] is pd.NA
    assert result.loc[0, "nam"] is pd.NA
    assert result["nam"].dtype == pd.BooleanDtype()
