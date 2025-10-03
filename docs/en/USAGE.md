# Usage Guide

This manual explains how to run the ChEMBL data acquisition pipelines and their supporting utilities.

## General CLI Pattern

Each entity pipeline is exposed as a console script (e.g., `get-document-data`) and can also be run as a Python module (e.g., `python -m scripts.get_document_data`).

### Common Arguments

Most pipeline scripts accept a shared set of arguments for consistency:

*   `--input`: Path to the input CSV file containing identifiers.
*   `--final-out`: Destination path for the final, processed CSV export. The legacy aliases `--output` and `--out` are still available but will be removed in a future version.
*   `--config`: Path to a YAML configuration file. Defaults to the packaged `config/config.yaml`.
*   `--log-level`: Sets the logging verbosity (e.g., `DEBUG`, `INFO`, `ERROR`).
*   `--limit`: Restricts processing to the first N identifiers.
*   `--offset`: Skips the first N identifiers before starting to process.
*   `--force`: Overwrites the output file if it already exists.
*   `--skip-existing`: Skips the pipeline run if the output file already exists.
*   `--print-config`: Prints the final, effective configuration (after all overrides) and exits.

### Pipeline-Specific Arguments

Individual pipelines may have additional arguments for controlling their specific logic, such as `--batch-size`, `--workers`, or sub-commands for different data sources. These are detailed in the sections below.

## Orchestrator (`get-data`)

The `get-data` script is the main entry point for running all pipelines in sequence. It provides a convenient way to perform a full data refresh with a single command.

```bash
get-data --base-path ./data \
    --input-dir input --output-dir output \
    --config config/config.yaml \
    --date 20240101 --limit 100 --log-level INFO
```

The orchestrator invokes the pipelines in the following order: `documents`, `targets`, `assays`, `testitems`, and `activities`. It constructs the necessary `--input` and `--final-out` paths for each step based on the shared directory arguments.

## Document Pipeline (`get-document-data`)

This pipeline retrieves and enriches publication metadata. It has three sub-commands: `chembl`, `pubmed`, and `all`.

| Sub-command | Description | Key Arguments |
| :--- | :--- | :--- |
| `chembl` | Retrieves document metadata from the ChEMBL API. | `--chunk-size`, `--timeout`, `--offset` |
| `pubmed` | Enriches with data from PubMed, Semantic Scholar, OpenAlex, and CrossRef. | `--sleep`, `--workers`, `--batch-size`, `--offset`, `--fallback-doi-csv` |
| `all` | Runs the `chembl` and `pubmed` pipelines, merging the results. | All arguments from `pubmed` and `chembl`. |

**Example:**

```bash
get-document-data all \
    --input data/input/document.csv \
    --final-out output/documents.csv \
    --limit 100 --offset 10 \
    --fallback-doi-csv data/input/doi_overrides.csv
```

## Target Pipeline (`get-target-data`)

This pipeline collects and merges target information from ChEMBL, UniProt, and IUPHAR. It supports keyboard aliases for sub-commands (e.g., `uniprot` can be typed as `унипрот` on a Russian keyboard layout).

### Staging Arguments

The target pipeline supports staged exports to separate raw data from the final normalized output:

*   `--raw-out`: Destination for the raw, combined dataset before normalization.
*   `--raw-format`: Format for the raw data (`csv` or `parquet`).
*   `--id-cols`: Identifier columns for deterministic sorting of the raw output.
*   `--no-reindex-raw`: Disables alphabetical re-indexing of columns in the raw output.
*   `--normalize-at-export` / `--no-normalize-at-export`: Controls whether the final normalization step is applied.

### Sub-commands

| Sub-command | Description | Key Arguments |
| :--- | :--- | :--- |
| `chembl` | Fetches target data from the ChEMBL API. | `--chunk-size`, `--timeout`, `--offset` |
| `uniprot` | Resolves UniProt accessions from local caches or the REST API. | `--data-dir`, `--offset` |
| `iuphar` | Maps UniProt IDs to IUPHAR families using local dictionaries. | `--target-csv`, `--family-csv`, `--offset` |
| `all` | Runs all three pipelines and merges the results. | `--chembl-out`, `--uniprot-out`, `--iuphar-out`, staging arguments. |

**Example:**

```bash
get-target-data all \
    --input data/input/target.csv \
    --final-out output/targets_final.csv \
    --raw-out output/targets_raw.parquet \
    --raw-format parquet \
    --id-cols target_chembl_id uniprot_id
```

## Assay Pipeline (`get-assay-data`)

Downloads assay metadata from ChEMBL.

**Arguments:**
*   `--batch-size` (default: 10)
*   `--timeout` (default: 30.0)
*   `--offset`
*   `--column`

**Example:**

```bash
get-assay-data --input data/input/assay.csv \
    --final-out output/assays.csv \
    --batch-size 20 --limit 100
```

## Activity Pipeline (`get-activity-data`)

Extracts and enriches activity data from ChEMBL.

**Arguments:**
*   `--batch-size` (default: 5)
*   `--workers` (default: 1)
*   `--timeout` (default: 30.0)
*   `--offset`
*   `--column`
*   `--dry-run`: Validates inputs and exits without making network calls.

**Example:**

```bash
get-activity-data --input data/input/activity.csv \
    --final-out output/activities.csv \
    --workers 4 --limit 500
```

## Test Item Pipeline (`get-testitem-data`)

Downloads and enriches compound data from ChEMBL and PubChem.

**Arguments:**
*   `--batch-size` (default: 1000)
*   `--timeout` (default: 30.0)
*   `--offset`
*   `--column`

**Example:**

```bash
get-testitem-data --input data/input/testitem.csv \
    --final-out output/testitems.csv \
    --limit 400
```

## Helper Utilities

The project includes several command-line helper utilities for diagnostics, QA, and offline workflows.

| Command | Purpose |
| :--- | :--- |
| `check-determinism` | Verifies that the deterministic CSV writers produce stable output. |
| `chunk-io` | Re-serializes CSV files in deterministic chunks. |
| `csv-utils` | Normalizes delimiters, quoting, and ordering for arbitrary CSV files. |
| `get-activities`| Generates synthetic activity data for testing and validation. |
| `get-document-type`| Applies publication-type classification heuristics to a document CSV. |
| `get-input-initialisation`| Combines Excel workbooks into canonical entity tables. |
| `mapper`| Provides an interactive tool for mapping identifiers between ChEMBL and UniProt. |
| `table-quality`| Generates a column-level quality and correlation report for a CSV file. |