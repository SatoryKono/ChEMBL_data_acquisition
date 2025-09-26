"""Schema definitions for activity data.

This module provides the :data:`ActivitiesSchema` object describing
expected structure of activity dataframes.
"""

from __future__ import annotations

from pathlib import Path
from typing import cast

import pandera.pandas as pa
from pandera.dtypes import DataType
import yaml

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


def _load_standard_type_values() -> tuple[str, ...]:
    """Return configured standard type values if available."""

    config_path = Path(__file__).resolve().parents[1] / "config.yaml"
    try:
        with config_path.open("r", encoding="utf8") as handle:
            config = yaml.safe_load(handle) or {}
    except (FileNotFoundError, yaml.YAMLError, TypeError):
        return ()

    metrics = (
        config
        .get("activity_enrichment", {})
        .get("action_type", {})
        .get("metrics", {})
    )
    if not isinstance(metrics, dict):
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


_STANDARD_TYPE_VALUES = _load_standard_type_values()
_STANDARD_TYPE_CHECKS: list[pa.Check] = []
if _STANDARD_TYPE_VALUES:
    _STANDARD_TYPE_CHECKS.append(pa.Check.isin(sorted(_STANDARD_TYPE_VALUES)))


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
