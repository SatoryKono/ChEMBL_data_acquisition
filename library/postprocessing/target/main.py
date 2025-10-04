"""Entry point for the organism-level target table post-processing."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd

from .cellularity import Cellularity
from .multifunctional import compute_multifunctional_table

__all__ = ["postprocess_target_table"]

_SOURCE_COLUMNS: tuple[str, ...] = (
    "target_chembl_id",
    "uniprot_id_primary",
    "organism",
    "taxon_id",
    "lineage_superkingdom",
    "lineage_phylum",
    "lineage_class",
)
_OUTPUT_COLUMNS: tuple[str, ...] = _SOURCE_COLUMNS + ("cellularity", "multifunctional_enzyme")


def _lower_text(value: Any) -> Any:
    if value is None:
        return None
    if isinstance(value, float) and pd.isna(value):
        return None
    text = str(value)
    return text.lower()


def _prepare_source_frame(frame: pd.DataFrame) -> pd.DataFrame:
    source_base = frame.loc[:, list(_SOURCE_COLUMNS)].copy()
    for column in ("lineage_superkingdom", "lineage_phylum", "lineage_class"):
        source_base[column] = source_base[column].map(_lower_text)
    return source_base


def _run_pipeline(frame: pd.DataFrame) -> pd.DataFrame:
    source_base = _prepare_source_frame(frame)
    with_cellularity = Cellularity.add_cellularity_smart(
        source_base,
        "taxon_id",
        "lineage_superkingdom",
        "lineage_class",
    )

    multifunctional_part = compute_multifunctional_table(frame)
    join_columns = ["target_chembl_id", "multifunctional_enzyme"]
    if not all(column in multifunctional_part.columns for column in join_columns):
        raise KeyError("multifunctional helper must provide required columns")

    merged = with_cellularity.merge(
        multifunctional_part.loc[:, join_columns],
        on="target_chembl_id",
        how="left",
        sort=False,
    )
    return merged.loc[:, list(_OUTPUT_COLUMNS)]


def postprocess_target_table(input_path: str) -> str:
    """Run the Power Query-equivalent organism post-processing pipeline."""

    source_path = Path(input_path)
    frame = pd.read_csv(source_path, dtype=object, keep_default_na=False)
    processed = _run_pipeline(frame)

    output_path = source_path.with_name(f"organism.{source_path.name}")
    processed.to_csv(output_path, index=False, encoding="utf-8", lineterminator="\n")

    print(f"Postprocessed target table saved to: {output_path.name}")
    return str(output_path)
