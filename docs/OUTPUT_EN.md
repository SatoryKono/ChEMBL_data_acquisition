# Output Directory Structure

## Base location

Generated datasets are written beneath `local.io.output_dir` (default: `data/output`). When `--output` is omitted the tools call
`library.io.default_output_path`, producing files named `output_<input-stem>_<YYYYMMDD>.csv` inside the configured output
directory. The writer automatically creates parent directories when `local.io.exist_ok` is `true`.

Example layout:

```
data/output/
└── ChEMBL/
    └── processed/
        ├── activity.csv
        ├── activity.csv.meta.yaml
        ├── activity_failure_cases.csv
        ├── activity_quality_report_table.csv
        ├── activity_data_correlation_report_table.csv
        └── ...
```

## Metadata sidecars

Each CSV export is accompanied by `<name>.csv.meta.yaml` created via `library.metadata.write_meta_yaml`. The metadata contains:

* `generated_at` — ISO 8601 timestamp in UTC.
* `git_sha` — commit hash of the repository at runtime.
* `python_version` and `platform` — runtime details.
* `command` — exact CLI invocation.
* `config` — relevant configuration values with secrets masked.
* `inputs` — description of the source files and parameters.
* `stats` — counters for `rows_total`, `rows_kept`, `rows_dropped` and the `output_sha256` digest.
* `schema` — name of the validation schema applied to the dataset.

If the sidecar already exists the new metadata is merged, preserving manually added annotations.

## Validation artefacts

* When Pandera validation discovers inconsistent rows the failing records are saved to `<stem>_failure_cases.csv` next to the main
  CSV.
* `library.table_quality.analyze_table_quality` generates `<stem>_quality_report_table.csv` and
  `<stem>_data_correlation_report_table.csv`. CLI utilities place these files alongside the dataset, while
  `scripts/get_input_initialisation.py` stores them under `<output>/data_validity_report/`.

All reports are written using UTF-8 encoding and share the same deterministic ordering rules as the main exports.

## Deterministic exports

`library.io.write_csv` forwards to `library.csv_utils.write_csv_deterministic`, which sorts columns and rows by key columns to
ensure identical results on repeated runs. The helper honours `cfg.io.csv_sep`, `cfg.io.csv_encoding` and optional
`key_cols`/`col_order` arguments supplied by the pipeline.

## Pipeline metadata columns

All entity exports append two bookkeeping columns produced by `library.pipeline_metadata.add_pipeline_metadata` before schema
validation and CSV writing: `pipeline_version` captures the installed package version (or the value discovered in
`pyproject.toml`), while `timestamp_utc` stores the run start time as an ISO 8601 string.【F:library/pipeline_metadata.py†L24-L84】
The columns are declared in the validation schemas for activities, documents, test items and other tables, ensuring the values
are present even when the payload is empty.【F:schemas/activities.py†L52-L55】【F:schemas/documents.py†L111-L112】【F:schemas/testitems.py†L41-L42】

`pipeline_version` remains constant within a single run and is safe for equality joins across different tables exported during
the same execution. `timestamp_utc` reflects the orchestrator’s clock; downstream consumers should treat it as metadata rather
than a surrogate for record-level timestamps.

## Document classification columns

`scripts/get_document_data.py` enriches the merged metadata with deterministic publication scores and labels emitted by `library.document_pipeline.merge_metadata`.【F:scripts/get_document_data.py†L607-L676】【F:library/document_pipeline.py†L160-L208】 The following fields appear in `document.csv` and the associated validation schema:

| Column | Description |
| --- | --- |
| `publication_types_normalised` | Semicolon-separated list of distinct publication type tokens collected from ChEMBL, PubMed, Semantic Scholar, OpenAlex and CrossRef payloads. The sequence is sorted to guarantee reproducible diffs. |
| `publication_type_score_review` | Integer weight derived from the weighted voting of review-specific terms. The score is non-negative and grows with corroborating evidence. |
| `publication_type_score_experimental` | Integer weight summarising experimental evidence terms; follows the same deterministic weighting scheme as the review score. |
| `publication_type_score_unknown` | Integer weight associated with ambiguous or explicitly unknown labels. |
| `publication_class` | Final class selected from `review`, `experimental` or `unknown` once the weighted tallies are compared against the configured thresholds.【F:library/document_type_classifier.py†L7-L74】 |

Scores default to `0` for rows without recognised tokens. The classifier prefers the highest score that also clears the minimum thresholds, falling back to `unknown` when the signal is inconclusive, so dashboards must not treat the label as an indicator of curation quality.

## Activity bounds (`lower_value`, `upper_value`)

`activity.csv` now exposes canonical value ranges via the `lower_value` and `upper_value` columns produced after normalisation in `scripts/get_activity_data.py`. The pipeline works exclusively with ChEMBL `standard_*` fields so that all limits remain in the canonical units already validated by the schema. Priority is applied row-wise in the following order:

1. Explicit bounds from `standard_lower_value` and `standard_upper_value` when provided by the API.
2. Paired values such as `standard_value` + `standard_upper_value`, using the minimum as the lower bound and the maximum as the upper bound if one side was previously missing.
3. Relation-driven inference when `activity_bounds.enable_from_relation` is `true`: `=`/`≈`/`~` set both bounds, `>=` fills only `lower_value`, `<=` fills only `upper_value`, while `between`/`range` expects a second canonical number. Unknown relation markers are left empty and logged for diagnosis. Missing canonical values despite the presence of a raw `value` trigger an `activity_bounds_missing_standard_value` warning so that data issues can be addressed upstream.
4. Optional parsing of `±` expressions from `standard_text_value` when `activity_bounds.enable_from_uncertainty` is enabled, guarded by the same canonical-unit requirement.

Derived bounds are rounded to `activity_bounds.rounding_digits` decimal places (default `3`) and clamped to zero for concentration-like metrics when `activity_bounds.clamp_nonnegative` is `true`, using heuristics based on `standard_type`/`standard_units`. All operations preserve existing columns and honour the deterministic column ordering enforced by the schema.【F:scripts/get_activity_data.py†L1-L234】【F:config.yaml†L108-L147】【F:library/config.py†L358-L420】【F:schemas/activities.py†L32-L64】

## Document quality JSON report

`scripts/get_document_data.py` writes an additional `<stem>.quality.json` file after the CSV export. The helper
`library.document_pipeline.build_quality_report` emits a structured payload with the total row count, DOI coverage ratio,
publication class distribution and per-source error counters; `save_quality_report` serialises the mapping as UTF-8 JSON with
stable formatting for diffs.【F:scripts/get_document_data.py†L636-L674】【F:library/document_pipeline.py†L300-L356】 Use the JSON
artefact to monitor metadata completeness without loading the full CSV, e.g. alert when DOI coverage drops below an agreed
threshold.

The JSON document contains:

| Key | Type | Description |
| --- | --- | --- |
| `rows_total` | Integer | Number of rows in the exported dataframe. |
| `doi_coverage` | Float | Fraction of documents with a populated DOI (0.0–1.0 range). |
| `publication_class_counts` | Object | Mapping of `review` / `experimental` / `unknown` to row counts, defaulting to `unknown` when the label is absent. |
| `error_counts` | Object | Dictionary with counts of failed enrichments per upstream (`pubmed`, `semantic_scholar`, `openalex`, `crossref`). |

All keys are always present; missing categories are represented by zero counters so that monitoring dashboards can rely on stable schemas across runs.

## Housekeeping recommendations

* Keep historical runs in dated subdirectories (`YYYYMMDD/`) to simplify comparisons.
* Archive or compress obsolete artefacts to reclaim disk space; metadata sidecars retain sufficient provenance information.
* Monitor free space before long-running extractions, especially when running multiple pipelines in parallel.

## Test item exports

`scripts/get_testitem_data.py` produces `testitem.csv` plus the standard `*.meta.yaml`, optional
`*_failure_cases.csv`, and quality reports following the deterministic ordering rules described above.【F:scripts/get_testitem_data.py†L151-L299】
Each row combines ChEMBL fields (`molecule_chembl_id`, structure descriptors, lifecycle flags), PubChem
augmentation, and pipeline metadata, allowing the dataset to serve as the canonical compound dimension.【F:scripts/get_testitem_data.py†L36-L193】【F:schemas/testitems.py†L12-L31】

Before distributing the export, join it with the parent molecule catalogue to expose
`parent_molecule_chembl_id` for roll-ups. The mapping is stored in the JSON file configured at
`sources.chembl.molecule_catalog.cache_path` and loaded via
`library.molecule_catalog.load_parent_catalog`, which refreshes the cache from the ChEMBL API when needed.【F:config.yaml†L25-L33】【F:library/molecule_catalog.py†L43-L136】

An additional enrichment stage reads `dictionary/molecule_hierarchy.csv` and
`dictionary/molecule_catalog.csv` to populate the salt identifier and the
`natural_product`, `prodrug`, `polymer_flag` booleans before validation. The
pipeline detects salts when `parent_molecule_chembl_id` differs from
`molecule_chembl_id`, fills `salt_chembl_id` with the child identifier, and
normalises catalogue flags to the pandas nullable boolean dtype. Missing child
flags fall back to the parent entry when available, while absent parent/child
records are surfaced via warning events for troubleshooting.【F:scripts/get_testitem_data.py†L205-L233】【F:library/testitem_enrichment.py†L17-L216】

## Release notes

- Target finalisation now sources organism types directly from the
  auto-generated `dictionary/_Target/targets_type.csv`. The legacy
  organism lookup CSV is no longer bundled. Override
  `local.resources.targets_type_csv` when custom classifications are
  required.【F:library/target_postprocessing.py†L485-L608】【F:scripts/get_target_data.py†L1008-L1031】
- The `get_target_data.py all` sub-command exposes `--targets-type-csv` instead
  of the legacy `--organism-csv` flag. Update automation scripts that passed the
  old argument to keep using custom lookup tables.【F:scripts/get_target_data.py†L312-L371】
