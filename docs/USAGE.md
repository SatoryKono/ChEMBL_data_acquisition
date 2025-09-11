# Usage Examples

This guide demonstrates how to run the command line tools on the bundled
"smoke" datasets. Each example writes output next to the input file with a
prefix of `output_` followed by the input stem and current date.

All scripts call :func:`library.config.ensure_dirs` after loading the
configuration so that the configured output and cache directories exist before
processing begins.

## Activity data

```bash
python get_activity_data.py \
    --input data/input-smoke/activity.csv \
    --column activity_id
```

## Assay descriptions

```bash
python get_assay_data.py \
    --input data/input-smoke/assay.csv \
    --column assay_chembl_id
```

## Document metadata

```bash
python get_document_data.py \
    --input data/input-smoke/documents.csv \
    --column document_chembl_id
```


## Target data aggregation

```bash
python get_target_data.py \
    --input data/input-smoke/targets.csv \
    --column target_chembl_id
```

## Test item data enrichment

```bash
python get_testitem_data.py \
    --input data/input-smoke/testitem.csv \
    --column compound_chembl_id
```

## Input initialisation merging

```bash
python get_input_initialisation.py \
    --same-doc dictionary/classifications/assay_classification.csv \
    --all-doc dictionary/classifications/target_classification.csv \
    --out-dir data/output-smoke
```

Each command supports the common flags `--sep`, `--encoding` and
`--log-level`. The default output path mirrors the input location and can be
manually overridden with `--output` where available.

### Configuration overrides

Command line flags override values from `config.yaml`. Internally the scripts
use `library.cli.apply_config_overrides` to merge provided options into the
runtime configuration. For example, specifying `--sep` or `--encoding`
replaces `io.csv_sep` and `io.csv_encoding` respectively. Likewise
`--chunk-size` maps to `jobs.chunk_size`, `--timeout` to `api.timeout_read` and
`--log-level` to `log.level`.

## Table quality profiler

```bash
python table_quality_main.py \
    --input data/input-smoke/activity.csv \
    --table-name activity
```

The profiler writes `<table-name>_quality_report_table.csv` and
`<table-name>_data_correlation_report_table.csv` in the current directory.

## Running the tests

The repository contains a small regression suite using `pytest`:

```bash
pytest
```

Test data resides in `tests/data` and offers realistic examples of the CSV
structures expected by the command line tools.

## Code style checks

Run the standard formatting, linting and type checking tools before submitting
changes. Install the developer extras first:

```bash
pip install -r requirements-dev.txt
black get_*.py library mapper_main.py table_quality_main.py
ruff check get_*.py library mapper_main.py table_quality_main.py
mypy get_*.py library mapper_main.py table_quality_main.py
```
