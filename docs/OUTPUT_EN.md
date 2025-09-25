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

## Housekeeping recommendations

* Keep historical runs in dated subdirectories (`YYYYMMDD/`) to simplify comparisons.
* Archive or compress obsolete artefacts to reclaim disk space; metadata sidecars retain sufficient provenance information.
* Monitor free space before long-running extractions, especially when running multiple pipelines in parallel.
