import pandas as pd
import pytest
from pandera.errors import SchemaError

from library.schemas.testitems import TestitemsSchema


@pytest.mark.integration
def test_testitems_schema__valid_frame() -> None:
    frame = pd.DataFrame(
        {
            "molecule_chembl_id": ["CHEMBL1"],
            "parent_molecule_chembl_id": ["CHEMBL10"],
            "natural_product": pd.Series([pd.NA], dtype="boolean"),
        }
    )

    validated = TestitemsSchema.validate(frame)
    assert list(validated["molecule_chembl_id"]) == ["CHEMBL1"]


@pytest.mark.integration
def test_testitems_schema__missing_required_column_raises() -> None:
    frame = pd.DataFrame({"parent_molecule_chembl_id": ["CHEMBL10"]})

    with pytest.raises(SchemaError):
        TestitemsSchema.validate(frame)
