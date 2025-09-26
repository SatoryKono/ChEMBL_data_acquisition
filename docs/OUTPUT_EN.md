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

## Activity bounds (`lower_value`, `upper_value`)

`activity.csv` now exposes canonical value ranges via the `lower_value` and `upper_value` columns produced after normalisation in `scripts/get_activity_data.py`. The pipeline works exclusively with ChEMBL `standard_*` fields so that all limits remain in the canonical units already validated by the schema. Priority is applied row-wise in the following order:

1. Explicit bounds from `standard_lower_value` and `standard_upper_value` when provided by the API.
2. Paired values such as `standard_value` + `standard_upper_value`, using the minimum as the lower bound and the maximum as the upper bound if one side was previously missing.
3. Relation-driven inference when `activity_bounds.enable_from_relation` is `true`: `=`/`≈`/`~` set both bounds, `>=` fills only `lower_value`, `<=` fills only `upper_value`, while `between`/`range` expects a second canonical number. Unknown relation markers are left empty and logged for diagnosis. Missing canonical values despite the presence of a raw `value` trigger an `activity_bounds_missing_standard_value` warning so that data issues can be addressed upstream.
4. Optional parsing of `±` expressions from `standard_text_value` when `activity_bounds.enable_from_uncertainty` is enabled, guarded by the same canonical-unit requirement.

Derived bounds are rounded to `activity_bounds.rounding_digits` decimal places (default `3`) and clamped to zero for concentration-like metrics when `activity_bounds.clamp_nonnegative` is `true`, using heuristics based on `standard_type`/`standard_units`. All operations preserve existing columns and honour the deterministic column ordering enforced by the schema.【F:scripts/get_activity_data.py†L1-L234】【F:config.yaml†L108-L147】【F:library/config.py†L358-L420】【F:schemas/activities.py†L32-L64】

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
`parent_molecule_chembl_id` for roll-ups. The mapping resides in the SQLite cache at
`sources.chembl.molecule_catalog.sqlite_path`; `library.molecule_catalog.load_parent_catalog`
initialises the database (migrating any legacy JSON configured via
`cache_path`) and refreshes it from the ChEMBL API when necessary.【F:config.yaml†L25-L33】【F:library/molecule_catalog.py†L108-L189】
