# Usage Guide

## Shared CLI options

All `scripts/get_*_data.py` commands share a common interface:

| Option | Description |
| --- | --- |
| `--config` | Path to the YAML configuration file (default: `config.yaml`). |
| `--print-config` | Print the effective configuration after overrides and exit. |
| `--log-level` | Logging level (`DEBUG`, `INFO`, `WARNING`, `ERROR`, `CRITICAL`). |
| `--input` | Input CSV with identifiers (default: `input.csv`). |
| `--output` | Destination CSV. When omitted, a file named `output_<input-stem>_<YYYYMMDD>.csv` is created inside `local.io.output_dir`. |
| `--sep` | CSV delimiter forwarded to `cfg.io.csv_sep`. |
| `--encoding` | File encoding forwarded to `cfg.io.csv_encoding`. |
| `--column` | Name of the identifier column. Defaults are populated from the configuration during start-up. |
| `--batch-size` / `--chunk-size` | Maximum number of identifiers per API request (option name depends on the pipeline). |

Each parser may add domain-specific switches such as `--timeout`, `--limit` or `--dry-run`. After parsing, `apply_config_overrides`
loads `config.yaml`, applies environment variables, merges CLI overrides back into the configuration, and updates missing CLI
arguments with the final values.

Before any network calls the utilities invoke `library.config.ensure_dirs`, ensuring that `local.io.output_dir` and
`local.io.cache_dir` exist (subject to `local.io.exist_ok`).

## Activity data (`get_activity_data.py`)

```bash
python scripts/get_activity_data.py \
  --input data/input-smoke/activity.csv \
  --column activity_chembl_id \
  --batch-size 25 \
  --timeout 45
```

* Reads the column configured at `sources.chembl.pipelines.activity.column` (`activity_chembl_id` by default).
* Writes the main CSV, `*.meta.yaml`, optional `*_failure_cases.csv` and quality reports.
* Supports `--limit` to restrict the number of identifiers and `--dry-run` to validate inputs without API calls.

## Assay descriptions (`get_assay_data.py`)

```bash
python scripts/get_assay_data.py \
  --input data/input-smoke/assay.csv \
  --column assay_chembl_id \
  --batch-size 25
```

Fetches assay metadata from ChEMBL using the configured identifier column.

## Document metadata (`get_document_data.py`)

```bash
python scripts/get_document_data.py \
  --input data/input-smoke/documents.csv \
  --column document_chembl_id \
  --sources.chembl.pipelines.document.pubmed.batch_size 20
```

The script merges ChEMBL and PubMed sources. CLI overrides accept dotted paths, enabling fine-tuned adjustments such as increasing
the PubMed batch size.

## Target aggregation (`get_target_data.py`)

```bash
python scripts/get_target_data.py \
  --input data/input-smoke/targets.csv \
  --column target_chembl_id
```

Combines ChEMBL, UniProt and IUPHAR sources according to `sources.chembl.pipelines.target.*`.

## Test item enrichment (`get_testitem_data.py`)

```bash
python scripts/get_testitem_data.py \
  --input data/input-smoke/testitem.csv \
  --column molecule_chembl_id
```

Downloads compound-centric annotations for the supplied identifiers.

## Input initialisation (`get_input_initialisation.py`)

```bash
python scripts/get_input_initialisation.py \
  --same-doc data/input/ChEMBL/ChEMBL_same_document_20_05.xlsx \
  --all-doc data/input/ChEMBL/ChEMBL_all_10_05_step5.xlsx \
  --out-dir data/output/ChEMBL/processed
```

* Builds pair tables (`pairs_same_document.csv`, `pairs_independent.csv`, `pairs_non_independent.csv`).
* Produces entity-specific slices (`activity_*`, `assay_*`, `document_*`, `target_*`, `testitem_*`, `system_*`).
* Creates a `data_validity_report/` folder with quality reports for each exported table.

## Table quality profiler (`table_quality_main.py`)

```bash
python table_quality_main.py \
  --input data/input-smoke/activity.csv \
  --table-name activity
```

Generates `<table-name>_quality_report_table.csv` and `<table-name>_data_correlation_report_table.csv` using the CSV defaults from
`local.io`.

## Configuration overrides from the CLI

Any CLI option can target nested configuration keys using dotted notation. Examples:

```bash
# Increase the global ChEMBL rate limit for a single run
python scripts/get_activity_data.py --sources.chembl.api.rps 10

# Change the CSV delimiter without editing config.yaml
python scripts/get_assay_data.py --sep ';'
```

## Environment variables

Environment variables follow the `CHEMBL_DA__SECTION__...__KEY` pattern. For frequently used paths, short aliases such as
`CHEMBL_DA_OUTDIR` (maps to `local.io.output_dir`) are available. See `docs/CONFIG_EN.md` for the full list.

## Running tests

Create a virtual environment and install the development extras first:

```bash
python -m venv .venv
source .venv/bin/activate  # Windows: .venv\\Scripts\\activate
pip install .[dev]
```

Then execute the unit and smoke suites with `pytest`:

```bash
pytest
pytest tests/smoke
```

To capture code coverage for the main packages (`library/`, `scripts/` and
`activity_extraction_main.py`) run:

```bash
pytest --cov=library --cov=scripts --cov=activity_extraction_main \
       --cov-report=term-missing --cov-report=xml
```

The command prints uncovered lines in the terminal and generates
`coverage.xml`, which can be consumed by CI tools or IDEs.

Code quality checks can be executed locally via:

```bash
ruff check
black --check .
mypy
```
