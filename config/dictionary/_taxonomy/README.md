# Assay taxonomy lookup

The CSV aggregates organism-level metadata used to enrich assay exports with
``assay_group`` and ``assay_strain`` values. The table maps the ``tax_id`` values
that ChEMBL publishes for assays to the curated group/strain labels extracted
from the historical workbook.

| Column    | Description                                   |
|-----------|-----------------------------------------------|
| `tax_id`  | NCBI taxonomy identifier referenced by assays |
| `assay_group` | Normalised group label (e.g. `Human`).     |
| `assay_strain` | Strain label captured in the legacy ETL. |

The enrichment code expects the file at
``config/dictionary/_taxonomy/taxonomy.csv`` by default. Update this file when
refreshing dictionary exports and keep the README in sync with the source and
validation procedure.
