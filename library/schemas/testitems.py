"""Schema definitions for test item data."""

from __future__ import annotations

from typing import cast

import pandas as pd
from pandera.dtypes import DataType

from library._compat.pandera import Column, DataFrameSchema

PA_ANY = cast(DataType, None)

TestitemsSchema: DataFrameSchema = DataFrameSchema(
    {
        "raw.index": Column(int, required=True, nullable=False, coerce=True),
        "molecule_chembl_id": Column(str, required=True, nullable=True),
        "parent_molecule_chembl_id": Column(str, required=False, nullable=True),
        "salt_chembl_id": Column(str, required=False, nullable=True),
        "natural_product": Column(pd.BooleanDtype(), required=False, nullable=True),
        "prodrug": Column(pd.BooleanDtype(), required=False, nullable=True),
        "polymer_flag": Column(pd.BooleanDtype(), required=False, nullable=True),
        "black_box_warning": Column(PA_ANY, required=False, nullable=True),
        "first_approval": Column(PA_ANY, required=False, nullable=True),
        "max_phase": Column(str, required=False, nullable=True),
        "canonical_smiles": Column(str, required=False, nullable=True),
        "standard_inchi": Column(str, required=False, nullable=True),
        "standard_inchi_key": Column(str, required=False, nullable=True),
        "molecule_type": Column(str, required=False, nullable=True),
        "oral": Column(PA_ANY, required=False, nullable=True),
        "parenteral": Column(PA_ANY, required=False, nullable=True),
        "pref_name": Column(str, required=False, nullable=True),
        "pubchem_canonical_smiles": Column(str, required=False, nullable=True),
        # ``pubchem_cid`` may appear as either a string or an integer depending on
        # the data source; use a generic ``object`` dtype to accept mixed representations.
        "pubchem_cid": Column(object, required=False, nullable=True),
        "pubchem_inchi": Column(str, required=False, nullable=True),
        "pubchem_inchikey": Column(str, required=False, nullable=True),
        "pubchem_isomeric_smiles": Column(str, required=False, nullable=True),
        "pubchem_iupac_name": Column(str, required=False, nullable=True),
        "pubchem_molecular_formula": Column(str, required=False, nullable=True),
        "structure_type": Column(str, required=False, nullable=True),
        "topical": Column(PA_ANY, required=False, nullable=True),
        "pipeline_version": Column(str, required=False, nullable=True),
        "timestamp_utc": Column(str, required=False, nullable=True),
    }
)

"""DataFrameSchema: Validation schema for test items."""
