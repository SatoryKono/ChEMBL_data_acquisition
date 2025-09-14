"""Schema definitions for assay data.

This module defines :data:`AssaysSchema`, a
:class:`pandera.DataFrameSchema` describing the expected structure of the
``assay.csv`` input table.  The schema mirrors the columns distributed with
``data/input/assay.csv`` and assigns concrete dtypes.  Counts such as
``acts_per_assay_step5``, temporal fields (``month``, ``year`` and ``version``)
use :class:`int`, boolean indicators like ``cited_assay_corr`` are :class:`bool`,
and the remaining descriptive fields are :class:`str`.

Column overview
---------------
* ``assay_chembl_id`` (:class:`str`): Primary ChEMBL assay identifier.
* ``ASSAY_ID`` (:class:`str`): Secondary assay identifier from the source file.
* ``Target TYPE`` (:class:`str`): Category of the biological target.
* ``accession`` (:class:`str`): UniProt accession of the target protein.
* ``acts_per_assay_step5`` (:class:`int`): Number of activities per assay step 5.
* ``assay_cell_type`` (:class:`str`): Cell type used in the assay.
* ``assay_subcellular_fraction`` (:class:`str`): Sub-cellular fraction tested.
* ``assay_tissue`` (:class:`str`): Tissue or organ where the assay is performed.
* ``bao_format`` (:class:`str`): BAO format code.
* ``cited_assay_corr`` (:class:`bool`): Whether the assay is cited as correlated.
* ``description`` (:class:`str`): Textual description of the assay.
* ``document_chembl_id`` (:class:`str`): Identifier of the source document.
* ``error_assay_corr`` (:class:`bool`): Error flag for the correlation citation.
* ``higly_correlated_cit`` (:class:`bool`): Flag for highly correlated citations.
* ``isoform`` (:class:`str`): Protein isoform number or identifier.
* ``month`` (:class:`int`): Month in which the assay was reported.
* ``mutation`` (:class:`str`): Target mutation details.
* ``shuffled_cit`` (:class:`bool`): Indicator for shuffled citation.
* ``shuffled_target_assay`` (:class:`bool`): Indicator for shuffled target/assay pair.
* ``substrate_name`` (:class:`str`): Name of the substrate used in the assay.
* ``target_chembl_id`` (:class:`str`): ChEMBL identifier of the target.
* ``target_name`` (:class:`str`): Human-readable target name.
* ``version`` (:class:`int`): Internal version number of the assay.
* ``year`` (:class:`int`): Year of publication.
"""

from __future__ import annotations

import pandera.pandas as pa

AssaysSchema: pa.DataFrameSchema = pa.DataFrameSchema(
    {
        "assay_chembl_id": pa.Column(str, required=True),
        "ASSAY_ID": pa.Column(str, required=False),
        "Target TYPE": pa.Column(str, required=False),
        "accession": pa.Column(str, required=False),
        "acts_per_assay_step5": pa.Column(int, required=False),
        "assay_cell_type": pa.Column(str, required=False),
        "assay_subcellular_fraction": pa.Column(str, required=False),
        "assay_tissue": pa.Column(str, required=False),
        "bao_format": pa.Column(str, required=False),
        "cited_assay_corr": pa.Column(bool, required=False),
        "description": pa.Column(str, required=False),
        "document_chembl_id": pa.Column(str, required=False),
        "error_assay_corr": pa.Column(bool, required=False),
        "higly_correlated_cit": pa.Column(bool, required=False),
        "isoform": pa.Column(str, required=False),
        "month": pa.Column(int, required=False),
        "mutation": pa.Column(str, required=False),
        "shuffled_cit": pa.Column(bool, required=False),
        "shuffled_target_assay": pa.Column(bool, required=False),
        "substrate_name": pa.Column(str, required=False),
        "target_chembl_id": pa.Column(str, required=False),
        "target_name": pa.Column(str, required=False),
        "version": pa.Column(int, required=False),
        "year": pa.Column(int, required=False),
    }
)

"""pa.DataFrameSchema: Validation schema for assays."""
