# Usage Guide

This manual explains how to run the ChEMBL data acquisition pipelines and their
supporting utilities. Every section has a Russian counterpart in
[`docs/USAGE_RU.md`](./USAGE_RU.md).

## General CLI pattern

Each entity pipeline exposes an installed console script and can also be executed
via `python -m scripts.<name>` during development. Command line options fall into
three tiers:

1. **Shared options** provided by `library.cli.parser.add_common_arguments`:
   - `--input / --final-out` – input CSV and destination for the cleaned export.
     `--output` and `--out` remain as deprecated aliases that trigger a
     deprecation warning; prefer `--final-out`.
   - `--log-level` – logging verbosity (`DEBUG`, `INFO`, `WARNING`, `ERROR`).
   - `--sep`, `--encoding` – CSV delimiter and encoding (`utf-8-sig` by default).
   - `--base-path`, `--input-dir`, `--output-dir`, `--date` – shortcuts used by
     the orchestrator to build consistent folder layouts and default filenames.
   - `--force`, `--skip-existing` – overwrite or skip when the destination file
     is already present.
   - `--config` – alternative YAML configuration file (defaults to
     `config/config.yaml`).
   - `--print-config` – output the effective configuration (after env/CLI
     overrides) and exit.
2. **Shared pagination options** exposed by the individual pipelines: `--column`
   (identifier column name), `--batch-size` / `--chunk-size`, `--timeout`,
   `--limit`, `--offset`, `--workers` (for concurrent fetchers), and `--dry-run`
   (activity pipeline only).
3. **Pipeline-specific switches** that tailor behaviour (for example staging
   flags in the target pipeline or DOI fallbacks in the document pipeline).

All commands exit with a non-zero status on validation errors, IO issues, or
upstream API failures.

### Input templates

The repository ships minimal CSV templates under `data/input` (one file per
pipeline: `document.csv`, `target.csv`, `assay.csv`, `activity.csv`,
`testitem.csv`) and extended samples under `data/input/full`. Copy the relevant
file, fill in your identifiers, and point `--input` to the new location. If you
maintain your own seed lists, place them anywhere accessible to the CLI and
ensure the headers match [`docs/DATA_SCHEMA_EN.md`](./DATA_SCHEMA_EN.md).

## Orchestrator (`get-data`)

```
get-data --base-path /data/chembl \
    --input-dir input --output-dir exports \
    --config /data/chembl/config.yaml \
    --date 20250101 --limit 100 --log-level INFO
```

The orchestrator resolves shared directories, prepares arguments, and invokes the
pipelines in the following order: documents (`all` sub-command), targets (`all`),
assays, test items, activities. With the sample arguments above, inputs are
loaded from `/data/chembl/input`—copy the templates there or provide your own
directory. Each delegated CLI receives `--final-out` so the individual pipelines
write to the canonical destination without relying on deprecated aliases.
`--limit 0` skips execution, `--dry-run` prints scheduled steps without touching
the filesystem.

## Document pipeline (`get-document-data`)

Sub-commands:

| Mode | Description | Key options |
|------|-------------|-------------|
| `chembl` | Retrieve metadata from the ChEMBL API. | `--column`, `--chunk-size`, `--timeout`, `--limit`, `--offset`. |
| `pubmed` | Enrich with PubMed, Semantic Scholar, OpenAlex, and CrossRef. | `--column`, `--sleep`, `--workers`, `--batch-size`, `--limit`, `--offset`, `--openalex-rps`, `--crossref-rps`, `--fallback-doi-*`. |
| `all` | Run ChEMBL, merge external services, and export the consolidated table. | Same options as `pubmed` plus `--fallback-doi-*` for DOI overrides. |

Typical command:

```
get-document-data all \
    --input data/input/document.csv \
    --final-out output/documents_$(date +%Y%m%d).csv \
    --config config/config.yaml \
    --limit 500 --log-level INFO
```

The pipeline writes a deterministic CSV, `<name>.meta.yaml`,
`<name>_quality_report_table.csv`, `<name>_data_correlation_report_table.csv`,
and `<name>.quality.json` with DOI coverage statistics.

## Target pipeline (`get-target-data`)

Sub-commands accept aliases in both Latin and Cyrillic layouts (to ease Windows
keyboard switching): `chembl`, `uniprot`, `iuphar`, `all`.

> **UniProt identifier requirement**: the input CSV must contain a column with
> UniProt accessions. By default the pipeline expects `uniprot_id`. Alternative
> headers can be supplied through `target.all.uniprot_column`, or detected
> automatically when the column name includes keywords like `uniprot` or
> `accession`.

### Staging switches

Only the target pipeline currently honours the staging flags:

- `--raw-out` – destination for the combined dataset before final normalisation.
- `--raw-format` – `csv` (default) or `parquet` for the raw snapshot.
- `--id-cols` – identifier columns used for deterministic ordering when writing
  raw snapshots.
- `--no-reindex-raw` – preserve the column order emitted by upstream APIs.
- `--normalize-at-export / --no-normalize-at-export` – control whether the final
  CSV applies the normalisation layer (`--no-normalize-at-export` keeps the raw
  payload byte-for-byte; useful when inspecting discrepancies).

### Modes

| Mode | Description | Key options |
|------|-------------|-------------|
| `chembl` | Fetch targets from ChEMBL, normalise, validate, export. | `--column`, `--chunk-size`, `--timeout`, `--limit`, `--offset`. |
| `uniprot` | Resolve UniProt accessions via local caches or the REST API. | `--column`, `--limit`. |
| `iuphar` | Map UniProt IDs to IUPHAR families using bundled dictionaries. | `--limit`. |
| `all` | Run `chembl`, `uniprot`, and `iuphar` sequentially, merge the results, and export. | All options above plus `--id-cols` / staging switches. |

Example with raw snapshot:

```
get-target-data all \
    --input data/input/target.csv \
    --final-out output/targets_$(date +%Y%m%d).csv \
    --raw-out output/targets_raw_$(date +%Y%m%d).parquet \
    --raw-format parquet --id-cols target_chembl_id uniprot_id \
    --config config/config.yaml --log-level INFO
```

Cached replay of the production pipeline is available through
`library.utils.cli_tools.pipeline_targets_main` (see the helper table below).

## Assay pipeline (`get-assay-data`)

```
get-assay-data --input data/input/assay.csv \
    --final-out output/assays_$(date +%Y%m%d).csv \
    --batch-size 100 --timeout 60 --limit 200
```

The script downloads assay metadata, computes per-target counters, normalises and
validates the table, then emits the standard CSV + sidecar + quality reports.

## Activity pipeline (`get-activity-data`)

```
get-activity-data --input data/input/activity.csv \
    --final-out output/activities_$(date +%Y%m%d).csv \
    --batch-size 50 --workers 4 --timeout 60 --limit 500
```

`--dry-run` is unique to this pipeline: it validates CLI arguments and input
files, then exits without contacting remote services. Activity-specific
post-processing derives `lower_value`/`upper_value` using the rules described in
[`docs/OUTPUT_EN.md`](./OUTPUT_EN.md).

## Test item pipeline (`get-testitem-data`)

```
get-testitem-data --input data/input/testitem.csv \
    --final-out output/testitems_$(date +%Y%m%d).csv \
    --batch-size 1000 --timeout 60 --limit 400
```

After downloading molecules from ChEMBL, the pipeline enriches unique SMILES with
PubChem properties, merges the results, and exports the combined dataset.

## Helper utilities

Use these modules for diagnostics, QA, or offline workflows. Each exposes a
`main(argv: Sequence[str] | None = None) -> int` entry point.

| Module | Console command | Purpose |
|--------|-----------------|---------|
| `library.utils.cli_tools.check_determinism` | `check-determinism --log-level INFO` | Verify deterministic CSV writers by hashing sample outputs. |
| `library.utils.cli_tools.chunk_io_main` | `chunk-io --input data.csv --final-out copy.csv` | Re-serialise CSV files in deterministic chunks. |
| `library.utils.cli_tools.csv_utils_main` | `csv-utils --input data.csv --final-out clean.csv --sep ,` | Normalise delimiters, quoting, and ordering for arbitrary CSV files. |
| `library.utils.cli_tools.dtype_inspector_main` | `python -m library.utils.cli_tools.dtype_inspector_main` | Inspect pandas dtypes produced by the pipelines. |
| `library.utils.cli_tools.get_activities` | `get-activities --limit 10` | Emit synthetic activity rows to verify logging and CLI wiring. |
| `library.utils.cli_tools.get_document_type` | `get-document-type --input docs.csv` | Apply the bundled publication-type heuristics. |
| `library.utils.cli_tools.get_input_initialisation` | `get-input-initialisation --same-doc init.xlsx --all-doc pairs.xlsx` | Combine Excel workbooks into canonical entity/relationship tables. |
| `library.utils.cli_tools.mapper_main` | `mapper --input ids.csv --final-out mapped.csv --column target_chembl_id` | Interactive mapper for quick lookups. |
| `library.utils.cli_tools.mapper_batch_main` | `python -m library.utils.cli_tools.mapper_batch_main --input ids.csv` | Batch mapper suitable for scripts or QA jobs. |
| `library.utils.cli_tools.pipeline_targets_main` | `python -m library.utils.cli_tools.pipeline_targets_main --input targets.csv` | Replay the target pipeline using cached API responses; honours the staging flags documented above. |
| `library.utils.cli_tools.table_quality_main` | `table-quality --input data.csv --table-name chembl_targets --final-out reports/` | Produce table-quality summaries for arbitrary datasets. |

## Tips for large runs

- Prefer fully qualified output paths (`--final-out`, `--raw-out`) when running
  multiple exports on the same day; default filenames include `output.<stem>`
  plus the date prefix from `--date`.
- When automating, set `CHEMBL_DA_BASE_PATH` to control where user-specific
  directories are created (see `local.io` in `config/config.yaml`).
- The document pipeline can honour DOI overrides via `--fallback-doi-csv` with
  two extra flags: `--fallback-doi-pmid-column` and `--fallback-doi-value-column`.
- All pipelines respect the `CHEMBL_DA_LOG_LEVEL` environment variable. Combine
  it with `--log-level DEBUG` for verbose troubleshooting without editing YAML.
