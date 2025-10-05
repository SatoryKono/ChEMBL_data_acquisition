import pandas as pd
import pytest

from library.integration.iuphar_library import IUPHARData


@pytest.mark.unit
def test_map_uniprot_file__self_referential_family(tmp_path, caplog):
    target_df = pd.DataFrame(
        [
            {
                "target_id": "T1",
                "uniprot_id": "P12345",
                "target_name": "Cycle target",
                "family_id": "F1",
                "type": "Enzyme.Transferase",
                "synonyms": "",
                "gene_name": "",
                "hgnc_name": "",
                "hgnc_id": "",
            }
        ]
    )
    family_df = pd.DataFrame(
        [
            {
                "family_id": "F1",
                "parent_family_id": "F1",
                "family_name": "Family One",
                "type": "Enzyme.Transferase",
            }
        ]
    )

    data = IUPHARData(target_df=target_df, family_df=family_df)

    input_path = tmp_path / "input.csv"
    output_path = tmp_path / "output.csv"
    pd.DataFrame(
        [
            {
                "uniprot_id": "P12345",
                "GuidetoPHARMACOLOGY": "T1",
            }
        ]
    ).to_csv(input_path, index=False)

    caplog.set_level("WARNING")

    result = data.map_uniprot_file(input_path, output_path)

    assert result.loc[0, "IUPHAR_chain"].split(">") == ["F1"]
    assert result.loc[0, "full_id_path"] == "T1#F1"
    assert result.loc[0, "full_name_path"] == "Cycle target#Family One"
    assert any("family_chain_cycle" in record.message for record in caplog.records)
