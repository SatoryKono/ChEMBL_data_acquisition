import pandas as pd

from library.pipelines.testitem import pubchem


def test_normalise_pubchem_column_dtype__string_assignment_to_numeric_column():
    result = pd.DataFrame(
        {"pubchem_iupac_name": pd.Series([pd.NA, pd.NA], dtype="Int64")}
    )

    replacement = pd.Series(
        ["2-(1-benzofuran-2-yl)-4,5-dihydro-1H-imidazole", None],
        index=[0, 1],
        dtype="string",
    )

    pubchem._normalise_pubchem_column_dtype(result, "pubchem_iupac_name")

    non_null = replacement.dropna()
    result.loc[non_null.index, "pubchem_iupac_name"] = non_null

    assert result.loc[0, "pubchem_iupac_name"] == "2-(1-benzofuran-2-yl)-4,5-dihydro-1H-imidazole"
    assert result["pubchem_iupac_name"].dtype == "string"
