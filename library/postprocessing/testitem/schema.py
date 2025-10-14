"""Schema specification for test item postprocessing outputs."""

from __future__ import annotations

from collections.abc import Sequence

import pandas as pd

from library.postprocessing.common import DataFrameSchema, validate_schema
from library.schemas.testitem import TestitemSchema

BOOLEAN_COLUMNS: Sequence[str] = (
    "natural_product",
    "prodrug",
    "polymer_flag",
    "oral",
    "parenteral",
    "topical",
    pandera_schema=TestitemSchema,
)


TESTITEM_SCHEMA = DataFrameSchema(
    required_columns=("molecule_chembl_id",),
    optional_columns=(
        "parent_molecule_chembl_id",
        "salt_chembl_id",
        "natural_product",
        "prodrug",
        "polymer_flag",
        "black_box_warning",
        "first_approval",
        "max_phase",
        "canonical_smiles",
        "standard_inchi",
        "standard_inchi_key",
        "molecule_type",
        "oral",
        "parenteral",
        "pref_name",
        "pubchem_canonical_smiles",
        "pubchem_cid",
        "pubchem_inchi",
        "pubchem_inchikey",
        "pubchem_isomeric_smiles",
        "pubchem_iupac_name",
        "pubchem_molecular_formula",
        "structure_type",
        "topical",
        "pipeline_version",
        "timestamp_utc",
    ),
    dtypes={
        "molecule_chembl_id": "string",
        "parent_molecule_chembl_id": "string",
        "salt_chembl_id": "string",
        "natural_product": pd.BooleanDtype(),
        "prodrug": pd.BooleanDtype(),
        "polymer_flag": pd.BooleanDtype(),
        "canonical_smiles": "string",
        "standard_inchi": "string",
        "standard_inchi_key": "string",
        "molecule_type": "string",
        "oral": pd.BooleanDtype(),
        "parenteral": pd.BooleanDtype(),
        "pref_name": "string",
        "pubchem_canonical_smiles": "string",
        "pubchem_cid": object,
        "pubchem_inchi": "string",
        "pubchem_inchikey": "string",
        "pubchem_isomeric_smiles": "string",
        "pubchem_iupac_name": "string",
        "pubchem_molecular_formula": "string",
        "structure_type": "string",
        "topical": pd.BooleanDtype(),
        "pipeline_version": "string",
        "timestamp_utc": "string",
    },
    sort_by=("molecule_chembl_id",),
    column_order=(
        "molecule_chembl_id",
        "parent_molecule_chembl_id",
        "salt_chembl_id",
        "pref_name",
        "molecule_type",
        "structure_type",
        "max_phase",
        "natural_product",
        "prodrug",
        "polymer_flag",
        "oral",
        "parenteral",
        "topical",
        "black_box_warning",
        "first_approval",
        "canonical_smiles",
        "standard_inchi",
        "standard_inchi_key",
        "pubchem_cid",
        "pubchem_canonical_smiles",
        "pubchem_isomeric_smiles",
        "pubchem_inchi",
        "pubchem_inchikey",
        "pubchem_iupac_name",
        "pubchem_molecular_formula",
        "timestamp_utc",
        "pipeline_version",
    ),
)


def validate_testitems(df: pd.DataFrame, *, context: str = "testitems") -> pd.DataFrame:
    """Validate ``df`` using :data:`TESTITEM_SCHEMA`."""

    return validate_schema(df, TESTITEM_SCHEMA, context=context)


__all__ = ["BOOLEAN_COLUMNS", "TESTITEM_SCHEMA", "validate_testitems"]
