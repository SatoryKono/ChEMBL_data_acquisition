"""Pandera schema describing test item postprocessing outputs."""

from __future__ import annotations

from library._compat.pandera import pa

from .common import boolean_column, object_column, string_column

__all__ = ["TestitemSchema"]


TestitemSchema: pa.DataFrameSchema = pa.DataFrameSchema(
    {
        "molecule_chembl_id": string_column(required=True, nullable=False),
        "parent_molecule_chembl_id": string_column(required=False, nullable=True),
        "salt_chembl_id": string_column(required=False, nullable=True),
        "natural_product": boolean_column(required=False, nullable=True),
        "prodrug": boolean_column(required=False, nullable=True),
        "polymer_flag": boolean_column(required=False, nullable=True),
        "black_box_warning": string_column(required=False, nullable=True),
        "first_approval": string_column(required=False, nullable=True),
        "max_phase": string_column(required=False, nullable=True),
        "canonical_smiles": string_column(required=False, nullable=True),
        "standard_inchi": string_column(required=False, nullable=True),
        "standard_inchi_key": string_column(required=False, nullable=True),
        "molecule_type": string_column(required=False, nullable=True),
        "oral": boolean_column(required=False, nullable=True),
        "parenteral": boolean_column(required=False, nullable=True),
        "pref_name": string_column(required=False, nullable=True),
        "pubchem_canonical_smiles": string_column(required=False, nullable=True),
        "pubchem_cid": object_column(object, required=False, nullable=True),
        "pubchem_inchi": string_column(required=False, nullable=True),
        "pubchem_inchikey": string_column(required=False, nullable=True),
        "pubchem_isomeric_smiles": string_column(required=False, nullable=True),
        "pubchem_iupac_name": string_column(required=False, nullable=True),
        "pubchem_molecular_formula": string_column(required=False, nullable=True),
        "structure_type": string_column(required=False, nullable=True),
        "topical": boolean_column(required=False, nullable=True),
        "pipeline_version": string_column(required=False, nullable=True),
        "timestamp_utc": string_column(required=False, nullable=True),
    }
)
