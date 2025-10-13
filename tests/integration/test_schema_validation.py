import pandas as pd
import pytest
from pandera.errors import SchemaError

from library.schemas import AssaysSchema
from library.schemas.testitems import TestitemsSchema


@pytest.mark.integration
def test_testitems_schema__valid_frame() -> None:
    frame = pd.DataFrame(
        {
            "raw.index": [0],
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


@pytest.mark.integration
def test_assays_schema__accepts_group_and_strain_columns() -> None:
    frame = pd.DataFrame(
        {
            "assay_chembl_id": ["CHEMBL1"],
            "assay_group": ["Group"],
            "assay_strain": ["Strain"],
        }
    )

    validated = AssaysSchema.validate(frame)

    assert list(validated["assay_group"]) == ["Group"]
    assert list(validated["assay_strain"]) == ["Strain"]


@pytest.mark.integration
def test_assays_schema__missing_optional_group_and_strain() -> None:
    frame = pd.DataFrame({"assay_chembl_id": ["CHEMBL1"]})

    validated = AssaysSchema.validate(frame)

    assert "assay_group" not in validated.columns
    assert "assay_strain" not in validated.columns

    schema_columns = list(AssaysSchema.columns)
    group_index = schema_columns.index("assay_group")
    strain_index = schema_columns.index("assay_strain")
    assert group_index < strain_index
