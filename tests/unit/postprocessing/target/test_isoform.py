import pandas as pd
import pytest

from library.postprocessing.target import isoform


@pytest.mark.unit
def test_transform__falls_back_to_uniprot_columns():
    frame = pd.DataFrame(
        {
            "isoform_synonyms": ["Alpha|Beta"],
            "isoform_names": ["Alpha|Beta"],
            "isoform_ids": ["ID_A|ID_B"],
            "uniProtkbId": ["Q11111"],
        }
    )

    transform = isoform._transform(frame)

    result = transform.result
    assert result["uniprot_id_primary"].unique().tolist() == ["Q11111"]
    assert result["target_chembl_id"].unique().tolist() == ["Q11111"]
    assert result.isna().sum().sum() == 0


@pytest.mark.unit
def test_transform__empty_input_returns_empty_frame():
    frame = pd.DataFrame(columns=["isoform_synonyms", "isoform_names", "isoform_ids"])

    transform = isoform._transform(frame)

    result = transform.result
    assert result.empty
    assert list(result.columns) == list(isoform._OUTPUT_COLUMNS)


@pytest.mark.unit
def test_resolve_source_columns__includes_fallback_aliases_in_error():
    frame = pd.DataFrame(columns=["unrelated_column"])

    with pytest.raises(KeyError) as excinfo:
        isoform._resolve_source_columns(frame)

    message = str(excinfo.value)
    assert "uniProtkbId" in message
    assert "target_chembl_id" in message
