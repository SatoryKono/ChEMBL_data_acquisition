# Taxonomy dictionary

This directory ships the canonical lookup used by the assay enrichment
pipeline to translate raw `assay_tax_id` values into curated assay group and
strain labels.  The CSV is a distilled subset of the full reference taxonomy
produced from the October 2025 ChEMBL export.  It focuses on the most common
organisms observed in the public datasets and mirrors the structure consumed by
`library.postprocessing.assay_extended`.

## Files

- `taxonomy.csv` — mapping of `assay_tax_id` to curated `assay_group` and
  `assay_strain` values.  The file is encoded in UTF-8 with LF line endings and
  sorted by `assay_tax_id` for deterministic processing.

## Schema

| Column name      | Type   | Description                                                      |
|------------------|--------|------------------------------------------------------------------|
| `assay_tax_id`   | string | NCBI taxonomy identifier referenced by assay records.            |
| `assay_group`    | string | Normalised organism group presented in assay exports.            |
| `assay_strain`   | string | Curated strain or lineage label; empty when the source omits it. |

## Regeneration

Run `tools/build_dictionary_resources.py --resource taxonomy` after updating the
source taxonomy snapshot.  The script re-generates the CSV, sorts it, verifies
its checksum, and updates `config/dictionary/manifest.yaml` with the new hash.
