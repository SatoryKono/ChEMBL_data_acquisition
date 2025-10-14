"""Schema definitions for assay data.

This module defines :data:`AssaysSchema`, a
:class:`pandera.DataFrameSchema` describing the expected structure of the
``assay.csv`` input table.  The schema mirrors the columns distributed with
``data/input/assay.csv`` and allows flexible dtypes.  Counts such as
 ``acts_per_assay_step-5`` and temporal fields (``month``, ``year`` and
 ``version``) previously used strict :class:`int` dtypes while boolean flags
like ``cited_assay_corr`` were :class:`bool`.  These fields now omit explicit
dtype enforcement so that CSVs representing numbers or booleans as strings
remain valid.  All columns are nullable to accommodate missing values.

Column overview
---------------
* ``assay_chembl_id`` (:class:`str`): Primary ChEMBL assay identifier.
* ``accession`` (:class:`str`): UniProt accession of the target protein.
* ``assay_cell_type`` (:class:`str`): Cell type used in the assay.
* ``assay_subcellular_fraction`` (:class:`str`): Sub-cellular fraction tested.
* ``assay_tissue`` (:class:`str`): Tissue or organ where the assay is performed.
* ``bao_format`` (:class:`str`): BAO format code.
* ``description`` (:class:`str`): Textual description of the assay.
* ``document_chembl_id`` (:class:`str`): Identifier of the source document.
* ``error_assay_corr`` (:class:`object`): Error flag for the correlation citation.
* ``isoform`` (:class:`str`): Protein isoform number or identifier.
* ``month`` (:class:`object`): Month in which the assay was reported.
* ``mutation`` (:class:`str`): Target mutation details.
* ``substrate_name`` (:class:`str`): Name of the substrate used in the assay.
* ``target_chembl_id`` (:class:`str`): ChEMBL identifier of the target.
* ``version`` (:class:`object`): Internal version number of the assay.
* ``year`` (:class:`object`): Year of publication.
"""

from __future__ import annotations

from typing import Any, Final, cast

from library._compat.pandera import Column, DataFrameSchema

# ``None`` disables dtype enforcement while still allowing schema validation.
FLEXIBLE_DTYPE: Final[Any] = cast(Any, None)

AssaysSchema: DataFrameSchema = DataFrameSchema(
    {
        "assay_chembl_id": Column(str, required=True, nullable=True),
        "accession": Column(str, required=False, nullable=True),
        "assay_cell_type": Column(str, required=False, nullable=True),
        "assay_subcellular_fraction": Column(str, required=False, nullable=True),
        "assay_group": Column(str, required=False, nullable=True),
        "assay_tissue": Column(str, required=False, nullable=True),
        "assay_strain": Column(str, required=False, nullable=True),
        "bao_format": Column(str, required=False, nullable=True),
        "description": Column(str, required=False, nullable=True),
        "document_chembl_id": Column(str, required=False, nullable=True),
        "isoform": Column(FLEXIBLE_DTYPE, required=False, nullable=True),
        "mutation": Column(str, required=False, nullable=True),
        "target_chembl_id": Column(str, required=False, nullable=True),
        "year": Column(FLEXIBLE_DTYPE, required=False, nullable=True),
        "pipeline_version": Column(str, required=False, nullable=True),
        "timestamp_utc": Column(str, required=False, nullable=True),
    }
)

"""DataFrameSchema: Validation schema for assays."""
