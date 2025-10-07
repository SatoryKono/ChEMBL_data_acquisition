"""Main entry point for the organism helper post-processing pipeline."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from ...common.log import logger
from ..helpers import normalise_export_basename
from .cellularity import FetchLineageCallable, add_cellularity_smart
from .multifunctional import compute_multifunctional


def _lowercase_column(series: pd.Series) -> pd.Series:
    """Return the lower-cased variant of ``series`` preserving nulls."""

    return series.astype("string").str.lower()


def postprocess_target_table(
    input_path: str | Path,
    *,
    email: str | None = None,
    fetcher: FetchLineageCallable | None = None,
) -> str:
    """Replicate the Power Query post-processing for organism helpers."""

    path = Path(input_path)
    source = pd.read_csv(path, dtype=str, keep_default_na=False)

    base_columns = [
        "target_chembl_id",
        "uniprot_id_primary",
        "organism",
        "taxon_id",
        "lineage_superkingdom",
        "lineage_phylum",
        "lineage_class",
    ]
    missing_columns = [
        column for column in base_columns if column not in source.columns
    ]
    if missing_columns:
        logger.warning(
            "target_postprocess_missing_columns",
            path=str(path),
            columns=missing_columns,
        )

    source_base = source.reindex(columns=base_columns, fill_value="").copy()

    for column in ("lineage_superkingdom", "lineage_phylum", "lineage_class"):
        source_base[column] = _lowercase_column(source_base[column])

    with_cellularity = add_cellularity_smart(
        source_base,
        "taxon_id",
        "lineage_superkingdom",
        "lineage_class",
        email=email,
        fetcher=fetcher,
    )

    multifunctional_part = compute_multifunctional(source)
    joined = with_cellularity.merge(
        multifunctional_part[["target_chembl_id", "multifunctional_enzyme"]],
        on="target_chembl_id",
        how="left",
        sort=False,
    )

    base = normalise_export_basename(path)
    output_path = path.with_name(f"organism.{base}").with_suffix(".csv")
    joined.to_csv(output_path, index=False, encoding="utf-8", lineterminator="\n")
    print(f"Postprocessed target table saved to: {output_path.name}")
    return str(output_path)


__all__ = ["postprocess_target_table"]
