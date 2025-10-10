"""Entry points for the target post-processing helpers."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from library.pipelines.target.postprocessing import (
    finalise_targets,
    postprocess_targets,
)
from library.schemas.targets import CELLULARITY_COLUMN_NAME, TARGETS_COLUMN_ORDER

from . import cellularity, helpers, multifunctional

ORGANISM_COLUMNS: list[str] = [
    "target_chembl_id",
    "target_type",
    "unicellular_organism",
    "multifunctional_enzyme",
    "IUPHAR_class",
    "IUPHAR_subclass",
    "sortorder.target",
    "gene_index",
    "taxon_index",
]


def _build_organism_lookup(frame: pd.DataFrame) -> pd.DataFrame:
    """Return the helper organism lookup mirroring the legacy PQ workbook."""

    lookup = pd.DataFrame(index=frame.index.copy())
    lookup["target_chembl_id"] = frame["target_chembl_id"].astype("string")
    lookup["target_type"] = cellularity.normalise_series(frame[CELLULARITY_COLUMN_NAME])
    lookup["unicellular_organism"] = cellularity.unicellular_flag(lookup["target_type"])

    if "multifunctional_enzyme" in frame.columns:
        mf_series = multifunctional.normalise_multifunctional(
            frame["multifunctional_enzyme"]
        )
    else:
        mf_series = pd.Series(pd.NA, index=frame.index, dtype="boolean")
    lookup["multifunctional_enzyme"] = mf_series.fillna(False).astype("boolean")

    for column in (
        "IUPHAR_class",
        "IUPHAR_subclass",
        "sortorder.target",
        "gene_index",
        "taxon_index",
    ):
        if column in frame.columns:
            lookup[column] = (
                frame[column]
                .map(helpers.normalise_text)
                .replace("", "-")
                .astype("string")
            )
        else:
            lookup[column] = "-"

    lookup = lookup.reindex(columns=ORGANISM_COLUMNS)
    lookup = helpers.ensure_string_columns(
        lookup,
        [
            "target_chembl_id",
            "target_type",
            "IUPHAR_class",
            "IUPHAR_subclass",
            "sortorder.target",
            "gene_index",
            "taxon_index",
        ],
    )
    lookup = lookup.sort_values(by="target_chembl_id", kind="mergesort").reset_index(
        drop=True
    )
    return lookup


def _write_outputs(frame: pd.DataFrame, target_path: Path) -> str:
    output_path = target_path.with_name(f"postprocessed.{target_path.name}")
    helpers.write_csv(frame, output_path, columns=TARGETS_COLUMN_ORDER)

    organism_frame = _build_organism_lookup(frame)
    organism_path = target_path.with_name(f"organism.{target_path.name}")
    helpers.write_csv(organism_frame, organism_path, columns=ORGANISM_COLUMNS)

    print(f"Postprocessed target table saved to: {output_path}")
    return str(output_path)


def postprocess_target_table(input_path: str) -> str:
    """Post-process the merged target table and emit organism helpers."""

    source_path = Path(input_path)
    frame = helpers.read_csv_with_fallbacks(source_path)
    frame = helpers.ensure_string_columns(frame, frame.columns)

    processed = postprocess_targets(frame)
    processed = finalise_targets(processed)
    processed = cellularity.ensure_cellularity(processed)

    if CELLULARITY_COLUMN_NAME in processed.columns:
        processed[CELLULARITY_COLUMN_NAME] = cellularity.normalise_series(
            processed[CELLULARITY_COLUMN_NAME]
        )

    processed = helpers.fill_missing(processed, TARGETS_COLUMN_ORDER)

    return _write_outputs(processed, source_path)
