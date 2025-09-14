"""Schema definitions for test item data."""

from __future__ import annotations

import pandera.pandas as pa

TestitemsSchema: pa.DataFrameSchema = pa.DataFrameSchema(
    {
        "molecule_chembl_id": pa.Column(str, required=True),
        "black_box_warning": pa.Column(str, required=False, nullable=True),
        "first_approval": pa.Column(str, required=False, nullable=True),
        "max_phase": pa.Column(str, required=False, nullable=True),
        "molecule_structures.canonical_smiles": pa.Column(str, required=False, nullable=True),
        "molecule_structures.standard_inchi": pa.Column(str, required=False, nullable=True),
        "molecule_structures.standard_inchi_key": pa.Column(str, required=False, nullable=True),
        "molecule_type": pa.Column(str, required=False, nullable=True),
        "oral": pa.Column(str, required=False, nullable=True),
        "parenteral": pa.Column(str, required=False, nullable=True),
        "pref_name": pa.Column(str, required=False, nullable=True),
        "pubchem_canonical_smiles": pa.Column(str, required=False, nullable=True),
        "pubchem_cid": pa.Column(str, required=False, nullable=True),
        "pubchem_inchi": pa.Column(str, required=False, nullable=True),
        "pubchem_inchikey": pa.Column(str, required=False, nullable=True),
        "pubchem_isomeric_smiles": pa.Column(str, required=False, nullable=True),
        "pubchem_iupac_name": pa.Column(str, required=False, nullable=True),
        "pubchem_molecular_formula": pa.Column(str, required=False, nullable=True),
        "structure_type": pa.Column(str, required=False, nullable=True),
        "topical": pa.Column(str, required=False, nullable=True),

    }
)

"""pa.DataFrameSchema: Validation schema for test items."""
