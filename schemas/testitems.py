"""Schema definitions for test item data."""

from __future__ import annotations

import pandera.pandas as pa

TestitemsSchema: pa.DataFrameSchema = pa.DataFrameSchema(
    {
        "molecule_chembl_id": pa.Column(str, required=True),
        "black_box_warning": pa.Column(str, required=False),
        "first_approval": pa.Column(str, required=False),
        "max_phase": pa.Column(str, required=False),
        "molecule_structures.canonical_smiles": pa.Column(str, required=False),
        "molecule_structures.standard_inchi": pa.Column(str, required=False),
        "molecule_structures.standard_inchi_key": pa.Column(str, required=False),
        "molecule_type": pa.Column(str, required=False),
        "oral": pa.Column(str, required=False),
        "parenteral": pa.Column(str, required=False),
        "pref_name": pa.Column(str, required=False),
        "pubchem_canonical_smiles": pa.Column(str, required=False),
        "pubchem_cid": pa.Column(str, required=False),
        "pubchem_inchi": pa.Column(str, required=False),
        "pubchem_inchikey": pa.Column(str, required=False),
        "pubchem_isomeric_smiles": pa.Column(str, required=False),
        "pubchem_iupac_name": pa.Column(str, required=False),
        "pubchem_molecular_formula": pa.Column(str, required=False),
        "structure_type": pa.Column(str, required=False),
        "topical": pa.Column(str, required=False),

    }
)

"""pa.DataFrameSchema: Validation schema for test items."""
