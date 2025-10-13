"""Schema definitions for test item data."""

from __future__ import annotations

from typing import cast

import pandas as pd
from library._compat.pandera import pa
from pandera.dtypes import DataType

PA_ANY = cast(DataType, None)

TestitemsSchema: pa.DataFrameSchema = pa.DataFrameSchema(
    {
        "raw.index": pa.Column(int, required=True, nullable=False, coerce=True),
        "molecule_chembl_id": pa.Column(str, required=True, nullable=True),
        "parent_molecule_chembl_id": pa.Column(str, required=False, nullable=True),
        "salt_chembl_id": pa.Column(str, required=False, nullable=True),
        "natural_product": pa.Column(pd.BooleanDtype(), required=False, nullable=True),
        "prodrug": pa.Column(pd.BooleanDtype(), required=False, nullable=True),
        "polymer_flag": pa.Column(pd.BooleanDtype(), required=False, nullable=True),
        "black_box_warning": pa.Column(PA_ANY, required=False, nullable=True),
        "first_approval": pa.Column(PA_ANY, required=False, nullable=True),
        "max_phase": pa.Column(str, required=False, nullable=True),
        "canonical_smiles": pa.Column(str, required=False, nullable=True),
        "standard_inchi": pa.Column(str, required=False, nullable=True),
        "standard_inchi_key": pa.Column(str, required=False, nullable=True),
        "molecule_type": pa.Column(str, required=False, nullable=True),
        "oral": pa.Column(PA_ANY, required=False, nullable=True),
        "parenteral": pa.Column(PA_ANY, required=False, nullable=True),
        "pref_name": pa.Column(str, required=False, nullable=True),
        "pubchem_canonical_smiles": pa.Column(str, required=False, nullable=True),
        # ``pubchem_cid`` may appear as either a string or an integer depending on
        # the data source; use a generic ``object`` dtype to accept mixed representations.
        "pubchem_cid": pa.Column(object, required=False, nullable=True),
        "pubchem_inchi": pa.Column(str, required=False, nullable=True),
        "pubchem_inchikey": pa.Column(str, required=False, nullable=True),
        "pubchem_isomeric_smiles": pa.Column(str, required=False, nullable=True),
        "pubchem_iupac_name": pa.Column(str, required=False, nullable=True),
        "pubchem_molecular_formula": pa.Column(str, required=False, nullable=True),
        "structure_type": pa.Column(str, required=False, nullable=True),
        "topical": pa.Column(PA_ANY, required=False, nullable=True),
        "pipeline_version": pa.Column(str, required=False, nullable=True),
        "timestamp_utc": pa.Column(str, required=False, nullable=True),
    }
)

"""pa.DataFrameSchema: Validation schema for test items."""
