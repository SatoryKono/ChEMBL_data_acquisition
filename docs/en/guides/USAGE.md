# Usage Guide

This manual explains how to run the ChEMBL data acquisition pipelines and their
supporting utilities. Every section has a Russian counterpart in
[`../../ru/guides/USAGE.md`](../../ru/guides/USAGE.md).

## General CLI pattern

Each entity pipeline exposes an installed console script and can also be executed
via `python -m scripts.<name>` during development. Command line options fall into
three tiers:

1. **Shared options** provided by `library.cli.parser.add_common_arguments`:
   - `--input / --final-out` – input CSV and destination for the cleaned export.
     A hidden `--out` compatibility alias still exists for automation but emits a
     warning; documentation and new jobs should rely on `--final-out`.
   - `--log-level` – logging verbosity (`DEBUG`, `INFO`, `WARNING`, `ERROR`).
   - `--verbose` – shortcut enabling `DEBUG` logging without overriding
     configuration files.
   - `--run-id` – explicit identifier stamped into logs and metadata sidecars;
     defaults to a deterministic value or the `CHEMBL_DA_RUN_ID` environment
     variable when set.
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
upstream API failures. Each run writes a text log to
`logs/<script>_<YYYYMMDD>.log` in the repository root. Setting
`CHEMBL_DA_BASE_PATH` relocates the folder to `<base>/logs`. All entries use the
`[timestamp] [LEVEL] [logger] message` format so that warnings and errors can be
audited after the fact.

### Input templates

The repository ships minimal CSV templates under `data/input` (one file per
pipeline: `document.csv`, `target.csv`, `assay.csv`, `activity.csv`,
`testitem.csv`) and extended samples under `data/input/full`. Copy the relevant
file, fill in your identifiers, and point `--input` to the new location. If you
maintain your own seed lists, place them anywhere accessible to the CLI and
ensure the headers match [`../DATA_SCHEMA.md`](../DATA_SCHEMA.md).

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

## Document pipeline (`python scripts/get_document_data.py`)

Run the document workflow via the single entry point `python scripts/get_document_data.py --mode <chembl|pubmed|all>`. The `--mode`
flag replaces the legacy positional sub-commands while keeping the common CLI
arguments consistent with the other pipelines.

### Quick reference

| `--mode` | Purpose | Default column | Namespaced flags |
|----------|---------|----------------|------------------|
| `chembl` | Retrieve document metadata from the ChEMBL API. | `document_chembl_id` | `--chembl-chunk-size`, `--chembl-timeout` (aliases: `--chunk-size`, `--timeout`). |
| `pubmed` | Enrich with PubMed, Semantic Scholar, OpenAlex, and CrossRef. | `PMID` | `--pubmed-sleep`, `--pubmed-workers`, `--pubmed-batch-size`, `--openalex-rps`, `--crossref-rps`. |
| `all` | Execute the ChEMBL and PubMed stages sequentially, then merge the payloads. | `document_chembl_id` | Accepts both namespaces plus fallback DOI options. |

### Help excerpt

```
$ python scripts/get_document_data.py --mode all --help
...
  --chembl-chunk-size CHEMBL_CHUNK_SIZE, --chunk-size CHEMBL_CHUNK_SIZE
                        Maximum identifiers per ChEMBL request
  --pubmed-sleep PUBMED_SLEEP, --sleep PUBMED_SLEEP
                        Seconds to sleep between PubMed requests
  --pubmed-workers PUBMED_WORKERS, --workers PUBMED_WORKERS
                        Number of concurrent PubMed requests
  --pubmed-batch-size PUBMED_BATCH_SIZE, --batch-size PUBMED_BATCH_SIZE
                        Maximum PMIDs per PubMed request
  --chembl-timeout CHEMBL_TIMEOUT, --timeout CHEMBL_TIMEOUT
                        Timeout in seconds for each ChEMBL HTTP request
  --openalex-rps OPENALEX_RPS
                        Requests per second limit for OpenAlex
  --crossref-rps CROSSREF_RPS
                        Requests per second limit for CrossRef

Fallback DOI overrides:
  --fallback-doi-enabled
                        Enable lookup of DOI overrides from a CSV file
  --fallback-doi-path FALLBACK_DOI_PATH
                        CSV file containing DOI overrides keyed by PMID
  --fallback-doi-col-pmid FALLBACK_DOI_COL_PMID
                        Column containing PubMed identifiers in the fallback CSV
  --fallback-doi-col-doi FALLBACK_DOI_COL_DOI
                        Column containing DOI values in the fallback CSV
  --fallback-doi-delimiter FALLBACK_DOI_DELIMITER
                        Delimiter used when reading the fallback CSV (default: io.csv_sep)
  --fallback-doi-encoding FALLBACK_DOI_ENCODING
                        Encoding used for the fallback CSV (default: io.csv_encoding)
  --fallback-doi-overwrite
                        Allow replacing existing DOIs with fallback values
```

### Fallback DOI flags

| Flag | Default | Description |
|------|---------|-------------|
| `--fallback-doi-enabled` | Disabled | Activate reading overrides from a CSV file. |
| `--fallback-doi-path` | _required when enabled_ | Path to the CSV containing PMID → DOI rows. |
| `--fallback-doi-col-pmid` | `PMID` | Column storing PubMed identifiers inside the fallback CSV. |
| `--fallback-doi-col-doi` | `DOI` | Column storing DOI values inside the fallback CSV. |
| `--fallback-doi-delimiter` | `local.io.csv_sep` (defaults to `,`) | CSV delimiter applied when parsing overrides. |
| `--fallback-doi-encoding` | `local.io.csv_encoding` (defaults to `utf-8-sig`) | Encoding used for the fallback CSV. |
| `--fallback-doi-overwrite` | Disabled | Permit replacing existing DOIs with fallback values. |

### Example invocations

```bash
# ChEMBL-only export
python scripts/get_document_data.py --mode chembl \
    --input data/input/document.csv \
    --final-out output/documents_chembl.csv \
    --config config/config.yaml

# PubMed enrichment with throttled partner APIs
python scripts/get_document_data.py --mode pubmed \
    --input data/input/document.csv \
    --final-out output/documents_pubmed.csv \
    --config config/config.yaml \
    --openalex-rps 3 --crossref-rps 3

# Full merge with manual DOI corrections
python scripts/get_document_data.py --mode all \
    --input data/input/document.csv \
    --final-out output/documents_full.csv \
    --config config/config.yaml \
    --fallback-doi-enabled \
    --fallback-doi-path data/input/document_fallback.csv \
    --fallback-doi-overwrite
```

The pipeline writes a deterministic CSV, `<name>.meta.yaml`, `<name>_quality_report_table.csv`, `<name>_data_correlation_report_table.csv`, and `<name>.quality.json` containing DOI coverage metrics.

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
    --batch-size 20 --workers 4 --timeout 90 --limit 500
```

`--dry-run` is unique to this pipeline: it validates CLI arguments and input
files, then exits without contacting remote services. Activity-specific
post-processing derives `lower_value`/`upper_value` using the rules described in
[`../OUTPUT.md`](../OUTPUT.md).

## Test item pipeline (`get-testitem-data`)

```
get-testitem-data --input data/input/testitem.csv \
    --final-out output/testitems_$(date +%Y%m%d).csv \
    --batch-size 250 --timeout 90 --limit 400
```

After downloading molecules from ChEMBL, the pipeline enriches unique SMILES with
PubChem properties, merges the results, and exports the combined dataset. By default
only the dataset plus QA/correlation CSVs are written; pass `--emit-legacy-artifacts`
to restore failure-case dumps and metadata sidecars when debugging a run. 【F:library/pipelines/testitem/cli.py†L864-L1186】【F:library/cli/commands/get_testitem_data.py†L564-L738】

## Helper utilities

Use these modules for diagnostics, QA, or offline workflows. Each exposes a
`main(argv: Sequence[str] | None = None) -> int` entry point.

| Module | Console command | Purpose |
|--------|-----------------|---------|
| `library.utils.cli_tools.check_determinism` | `check-determinism --log-level INFO` | Verify deterministic CSV writers by hashing sample outputs. |
| `library.utils.cli_tools.chunk_io_main` | `chunk-io --input data.csv --final-out copy.csv` | Re-serialise CSV files in deterministic chunks with Unix newlines. |
| `library.utils.cli_tools.csv_utils_main` | `csv-utils --input data.csv --final-out clean.csv --sep ,` | Normalise delimiters, quoting, and ordering for arbitrary CSV files. |
| `library.utils.cli_tools.dtype_inspector_main` | `python -m library.utils.cli_tools.dtype_inspector_main` | Inspect pandas dtypes produced by the pipelines. |
| `library.utils.cli_tools.get_activities` | `python scripts/get_activities.py --limit 10 --dry-run` | Emit synthetic activity rows to verify logging and CLI wiring; defaults to the `activity_id` column. |
| `library.utils.cli_tools.get_activities` | `get-activities --limit 10 --final-out output/activities_smoke.csv` | Emit synthetic activity rows and write deterministic CSV + `.meta.yaml` artefacts for smoke verification; defaults to the `activity_id` column. |
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
  it with `--verbose` (or `--log-level DEBUG`) for verbose troubleshooting
  without editing YAML.
