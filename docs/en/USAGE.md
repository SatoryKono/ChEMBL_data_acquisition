# Usage Guide

This manual explains how to run the ChEMBL data acquisition pipelines and their supporting utilities. Every section has a Russian counterpart in [`../../ru/USAGE.md`](../../ru/USAGE.md).

## General CLI Pattern

All command-line tools share a common set of arguments for consistency. They generally fall into three categories:

1.  **Shared I/O and Configuration Options:**
    *   `--input`: Path to the input CSV file containing identifiers.
    *   `--final-out`: Path to the final, cleaned output file. The legacy aliases `--output` and `--out` are still available but will be removed in a future version.
    *   `--config`: Path to a YAML configuration file. Defaults to the packaged `config/config.yaml`.
    *   `--log-level`: Sets the logging verbosity (e.g., `INFO`, `DEBUG`).
    *   `--print-config`: Prints the final configuration after all overrides and exits.

2.  **Shared Execution Control Options:**
    *   `--limit`: Limits the number of records to process. Setting `--limit 0` is a useful way to test configuration without running the full pipeline.
    *   `--batch-size` / `--chunk-size`: Number of records to process in a single API call or batch.
    *   `--workers`: Number of parallel workers for concurrent operations.

3.  **Pipeline-Specific Options:**
    *   These are arguments unique to a specific pipeline, such as `--raw-out` for the target pipeline or `--fallback-doi-csv` for the document pipeline.

---

## Main Pipelines

These are the core scripts for fetching and processing data for major ChEMBL entities.

### Orchestrator (`get-data`)

This is the main entry point for running all data acquisition pipelines sequentially. It ensures that a consistent configuration is passed to each step.

```bash
get-data --base-path /data/chembl \
    --input-dir seeds --output-dir exports \
    --config /data/chembl/config.yaml \
    --date 20250101 --limit 100 --log-level INFO
```

### Document Pipeline (`get-document-data`)

Fetches and enriches publication data from ChEMBL, PubMed, CrossRef, and other sources.

```bash
get-document-data all \
    --input seeds/document_ids.csv \
    --final-out output/documents.csv \
    --limit 500
```

### Target Pipeline (`get-target-data`)

Fetches and enriches target data from ChEMBL, UniProt, and IUPHAR. This pipeline has unique staging flags for more granular control over outputs.

**Staging Flags:**
*   `--raw-out`: Saves the intermediate, pre-normalized data to a separate file.
*   `--raw-format`: The format for the raw output (`csv` or `parquet`).
*   `--id-cols`: Identifier columns for deterministic sorting.

**Mode Options:** Every sub-command (`chembl`, `uniprot`, `iuphar`, `all`) accepts the same selector and execution flags:

| Option | Purpose | Default values |
| --- | --- | --- |
| `--column` | Input column used to select identifiers. | `chembl`: `target_chembl_id`; `uniprot`/`iuphar`: `uniprot_id`; `all`: `target_chembl_id`. |
| `--chunk-size` | Number of identifiers requested per API batch. | `chembl`/`all`: `5`; `uniprot`/`iuphar`: `100`. |
| `--timeout` | Network timeout in seconds. | `30.0` for all modes. |
| `--limit` | Maximum number of records to process; omit to process all rows. | `null` (no limit). |
| `--offset` | Number of rows to skip before processing. | `0` for all modes. |

The `all` orchestrator forwards these values to the ChEMBL stage and exposes per-pipeline overrides prefixed with the mode name: `--chembl-chunk-size`, `--uniprot-timeout`, `--iuphar-limit`, etc. The merge step still uses `--uniprot-column` to choose the UniProt identifier from the ChEMBL export.

```bash
get-target-data all \
    --input seeds/target_ids.csv \
    --final-out output/targets_final.csv \
    --raw-out output/targets_raw.parquet \
    --raw-format parquet
```

### Assay Pipeline (`get-assay-data`)

Downloads and processes assay metadata from ChEMBL.

```bash
get-assay-data --input seeds/assay_ids.csv \
    --final-out output/assays.csv \
    --limit 200
```

### Activity Pipeline (`get-activity-data`)

Extracts activity data and calculates normalized value bounds.

```bash
get-activity-data --input seeds/activity_ids.csv \
    --final-out output/activities.csv \
    --limit 500
```

### Test Item Pipeline (`get-testitem-data`)

Retrieves molecule data and enriches it with properties from PubChem.

```bash
get-testitem-data --input seeds/molecule_ids.csv \
    --final-out output/testitems.csv \
    --limit 400
```

---

## Helper Utilities

These are supporting tools for diagnostics, data manipulation, and quality control.

| Command | Purpose | Example |
|---|---|---|
| `check-determinism` | Compares two CSV files to ensure they are identical. | `check-determinism --input a.csv --previous b.csv` |
| `chunk-io` | Reads and re-writes a CSV file in deterministic chunks. | `chunk-io --input data.csv --final-out copy.csv` |
| `csv-utils` | Normalizes the formatting (delimiters, quoting) of a CSV file. | `csv-utils --input data.csv --final-out clean.csv` |
| `dtype-inspector` | Inspects and reports the pandas dtypes for data produced by pipelines. | `python -m library.utils.cli_tools.dtype_inspector_main` |
| `get-activities` | Generates synthetic activity data for testing purposes. | `get-activities --limit 10 --dry-run` |
| `get-document-type` | Applies heuristics to classify publication types. | `get-document-type --input docs.csv` |
| `get-input-initialisation` | Merges Excel workbooks into canonical input files. | `get-input-initialisation --same-doc init.xlsx` |
| `mapper` | Maps identifiers between ChEMBL and UniProt. | `mapper --input ids.csv --final-out mapped.csv` |
| `table-quality` | Generates a quality and correlation report for a CSV file. | `table-quality --input data.csv --table-name my_data` |
| `pipeline-targets-main` | Replays the target pipeline using cached data, avoiding network calls. | `python -m library.utils.cli_tools.pipeline_targets_main --input targets.csv` |