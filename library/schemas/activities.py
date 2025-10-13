"""Schema definitions for activity data.

This module provides the :data:`ActivitiesSchema` object describing
expected structure of activity dataframes.
"""

from __future__ import annotations

from collections.abc import Mapping
from typing import cast

from pandera.dtypes import DataType

from library._compat.pandera import pa

PA_ANY = cast(DataType, None)

# Definition of the schema describing the activities table.
_ALLOWED_ACTION_TYPES = {
    "PAM",
    "NAM",
    "activation",
    "inhibition",
    "binding",
    "triaged",
    "unknown",
}


def _derive_standard_type_values(metrics: Mapping[str, str] | None) -> tuple[str, ...]:
    """Derive permissible ``standard_type`` values from *metrics* keys."""

    if not metrics:
        return ()

    values: set[str] = set()
    for key in metrics.keys():
        token = str(key).strip()
        if not token:
            continue
        variants = {
            token,
            token.lower(),
            token.upper(),
            token.capitalize(),
        }
        if token[0].isalpha():
            variants.add(token[0].upper() + token[1:])
        values.update(filter(None, variants))
    return tuple(sorted(values))


_STANDARD_TYPE_CHECKS: list[pa.Check] = []


ActivitiesSchema: pa.DataFrameSchema = pa.DataFrameSchema(
    {
        "activity_id": pa.Column(PA_ANY, required=True, nullable=True),
        "molecule_chembl_id": pa.Column(str, required=True, nullable=True),
        "assay_chembl_id": pa.Column(str, required=True, nullable=True),
        "activity_comment": pa.Column(str, required=False, nullable=True),
        "assay_description": pa.Column(str, required=False, nullable=True),
        "assay_variant_accession": pa.Column(str, required=False, nullable=True),
        "assay_variant_mutation": pa.Column(str, required=False, nullable=True),
        "bao_format": pa.Column(str, required=False, nullable=True),
        "bao_label": pa.Column(str, required=False, nullable=True),
        "data_validity_comment": pa.Column(str, required=False, nullable=True),
        "data_validity_description": pa.Column(str, required=False, nullable=True),
        "document_chembl_id": pa.Column(str, required=False, nullable=True),
        "pchembl_value": pa.Column(object, required=False, nullable=True),
        "potential_duplicate": pa.Column(PA_ANY, required=False, nullable=True),
        "qudt_units": pa.Column(str, required=False, nullable=True),
        "record_id": pa.Column(PA_ANY, required=False, nullable=True),
        "relation": pa.Column(object, required=False, nullable=True),
        "src_assay_id": pa.Column(PA_ANY, required=False, nullable=True),
        "src_id": pa.Column(PA_ANY, required=False, nullable=True),
        "standard_relation": pa.Column(str, required=False, nullable=True),
        "standard_units": pa.Column(str, required=False, nullable=True),
        "type": pa.Column(str, required=False, nullable=True),
        "units": pa.Column(str, required=False, nullable=True),
        "value": pa.Column(object, required=False, nullable=True),
        "standard_type": pa.Column(
            str,
            checks=_STANDARD_TYPE_CHECKS,
            required=False,
            nullable=True,
        ),
        "standard_value": pa.Column(
            float, pa.Check.ge(0), required=True, nullable=True, coerce=True
        ),
        "lower_value": pa.Column(float, required=False, nullable=True, coerce=True),
        "upper_value": pa.Column(float, required=False, nullable=True, coerce=True),
        "activity_properties": pa.Column(str, required=False, nullable=True),
        "action_type": pa.Column(
            str,
            pa.Check.isin(sorted(_ALLOWED_ACTION_TYPES)),
            required=False,
            nullable=True,
        ),
        "properties_hash": pa.Column(str, required=False, nullable=True),
        "pipeline_version": pa.Column(str, required=False, nullable=True),
        "timestamp_utc": pa.Column(str, required=False, nullable=True),
        #    "pA_value": pa.Column(float, pa.Check.in_range(-14, 14), required=False),
    }
)

"""pa.DataFrameSchema: Validation schema for activities."""


def configure_activity_schema(metrics: Mapping[str, str] | None) -> pa.DataFrameSchema:
    """Update :data:`ActivitiesSchema` based on *metrics* configuration."""

    values = _derive_standard_type_values(metrics)

    _STANDARD_TYPE_CHECKS.clear()
    if values:
        _STANDARD_TYPE_CHECKS.append(pa.Check.isin(sorted(values)))

    # ``DataFrameSchema`` stores a shallow copy of the provided list, therefore
    # we synchronise the column checks explicitly to reflect runtime changes.
    ActivitiesSchema.columns["standard_type"].checks = list(_STANDARD_TYPE_CHECKS)
    return ActivitiesSchema
