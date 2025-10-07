"""Target postprocessing pipeline."""
from .schema import TARGET_SCHEMA, validate_targets
from .steps import (
    PIPELINE_STEPS,
    enrich_target_synonyms,
    finalize_target_records,
    normalize_target_fields,
    run_target_pipeline,
)

__all__ = [
    "PIPELINE_STEPS",
    "TARGET_SCHEMA",
    "enrich_target_synonyms",
    "finalize_target_records",
    "normalize_target_fields",
    "run_target_pipeline",
    "validate_targets",
]
