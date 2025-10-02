# ChEMBL Data Acquisition Utilities

> **Project version:** 0.2.0 (2025-10-02)
>
> Comprehensive English/Russian documentation lives in the [`docs/`](docs/)
> directory. Each guide mirrors its translated counterpart and is kept in sync
> on every release. The consolidated change history is tracked in
> [`CHANGELOG.md`](CHANGELOG.md).

## Overview

The toolkit orchestrates deterministic downloads, normalisation and export of
ChEMBL-linked datasets (activities, assays, documents, targets, test items) with
strong validation guarantees. Highlights include:

- Unified CLI surface with deterministic chunked IO, schema validation and
  streaming-friendly CSV/Parquet writers.
- Configuration-first design (`config/config.yaml`, environment overrides,
  CLI flags) backed by pydantic models and JSON Schema.
- Extensive helper libraries for interacting with ChEMBL, PubChem, UniProt,
  CrossRef, PubMed, Semantic Scholar and IUPHAR.
- End-to-end QA via pytest, type checks, deterministic regression harnesses and
  packaging verification.

## Documentation map

| Topic | English | Russian |
|-------|---------|---------|
| Configuration | [docs/CONFIG_EN.md](docs/CONFIG_EN.md) | [docs/CONFIG_RU.md](docs/CONFIG_RU.md) |
| CLI utilities | [docs/CLI_TOOLS.md](docs/CLI_TOOLS.md) | [docs/CLI_TOOLS_RU.md](docs/CLI_TOOLS_RU.md) |
| Data schemas | [docs/DATA_SCHEMA_EN.md](docs/DATA_SCHEMA_EN.md) | [docs/DATA_SCHEMA_RU.md](docs/DATA_SCHEMA_RU.md) |
| ETL overview | [docs/ETL_PROCESS_EN.md](docs/ETL_PROCESS_EN.md) | [docs/ETL_PROCESS_RU.md](docs/ETL_PROCESS_RU.md) |
| Pipeline flow | [docs/ETL_DATA_FLOW_EN.md](docs/ETL_DATA_FLOW_EN.md) | [docs/ETL_DATA_FLOW_RU.md](docs/ETL_DATA_FLOW_RU.md) |
| Outputs | [docs/OUTPUT_EN.md](docs/OUTPUT_EN.md) | [docs/OUTPUT_RU.md](docs/OUTPUT_RU.md) |
| QA & release | [docs/QA_PROCESS_EN.md](docs/QA_PROCESS_EN.md) | [docs/QA_PROCESS_RU.md](docs/QA_PROCESS_RU.md) |
| Usage recipes | [docs/USAGE_EN.md](docs/USAGE_EN.md) | [docs/USAGE_RU.md](docs/USAGE_RU.md) |
| Summary | [docs/SUMMARY.EN.md](docs/SUMMARY.EN.md) | [docs/SUMMARY.RU.md](docs/SUMMARY.RU.md) |

## Supported CLI pipelines

Console entry points are installed via `pip install .` and mirror the
`python -m …` modules.

| Console script | Module | Purpose |
|----------------|--------|---------|
| `get-data` | `scripts.get_data` | End-to-end orchestration across all pipelines |
| `get-activity-data` | `scripts.get_activity_data` | Activity export with enrichment and schema validation |
| `get-assay-data` | `scripts.get_assay_data` | Assay metadata collection and post-processing |
| `get-document-data` | `scripts.get_document_data` | Document aggregation for PubMed, CrossRef, Semantic Scholar |
| `get-target-data` | `scripts.get_target_data` | Target aggregation with staged raw snapshots (`--raw-out`) |
| `get-testitem-data` | `scripts.get_testitem_data` | Test item enrichment pipeline |
| `get-document-type` | `library.utils.cli_tools.get_document_type` | Document classification helper |
| `get-input-initialisation` | `library.utils.cli_tools.get_input_initialisation` | Excel workbook merger |
| `csv-utils` | `library.utils.cli_tools.csv_utils_main` | CSV housekeeping utilities |
| `mapper` | `library.utils.cli_tools.mapper_main` | Identifier mapping helpers |
| `table-quality` | `library.utils.cli_tools.table_quality_main` | Profiling and QA dashboards |
| `chunk-io` | `library.utils.cli_tools.chunk_io_main` | Chunked IO smoke tests |
| `get-activities` | `library.utils.cli_tools.get_activities` | Minimal activity fetch example |
| `check-determinism` | `library.utils.cli_tools.check_determinism` | Determinism regression harness |

The `--raw-out`, `--raw-format`, `--id-cols`, `--no-reindex-raw`, and
`--normalize-at-export/--no-normalize-at-export` switches are fully supported by
`get-target-data`. Other pipelines accept the flags but ignore them until their
raw snapshot export layers are implemented.

## Installation

1. **Prepare environment**
   ```bash
   python -m pip install --upgrade pip setuptools wheel
   python -m venv .venv
   source .venv/bin/activate  # Windows: .venv\\Scripts\\activate
   ```

2. **Install dependencies**
   ```bash
   pip install -r requirements-lock.txt
   pip install -e .[dev]
   pre-commit install
   ```

3. **Validate setup**
   ```bash
   pre-commit run --all-files
   pytest
   python -m library.utils.cli_tools.check_determinism --limit 10 --log-level INFO
   ```

The lock file keeps local development aligned with CI. When dependency ranges in
`pyproject.toml` change, regenerate it from a clean environment using
`pip install .[dev]` followed by `pip freeze > requirements-lock.txt`.

## Configuration

- Defaults live in [`config/config.yaml`](config/config.yaml) and conform to
  [`config.schema.json`](config.schema.json).
- Local overrides: create `config/config.local.yaml` or point `--config` to a
  custom file; environment variables use the `CHEMBL_DA__...` prefix.
- Common CLI flags: `--input`, `--output` (primary destination), `--final-out`
  (staged target exports), `--config`, `--print-config`, `--sep`, `--encoding`,
  `--chunk-size` / `--batch-size`, `--log-level`.
- Path macros (`$CHEMBL_DA_BASE_PATH`) resolve to the CLI `--base-path` or the
  `CHEMBL_DA_BASE_PATH` environment variable.

Refer to the configuration manual for the complete list of sections, rate
limiter settings and enrichment toggles.

## Data outputs

`docs/OUTPUT_EN.md` describes the canonical structure of activity, assay,
document, target and test item exports, including column definitions and file
layout. Raw dumps produced via `--raw-out` retain API column order while final
exports enforce normalised schemas validated by [`schemas/`](schemas/).

## Testing & QA

- Unit and integration tests: `pytest`, focus areas include deterministic IO,
  schema validation, mapper pipelines and packaging.
- Static analysis: `pre-commit`, `ruff`, `black`, `mypy`.
- Regression harnesses: `python -m library.utils.cli_tools.check_determinism`,
  `python -m library.utils.cli_tools.mapper_batch_main` for bulk mapping runs.
- QA procedures are listed in [docs/QA_PROCESS_EN.md](docs/QA_PROCESS_EN.md).

## Release management

1. Update version numbers in `pyproject.toml`, READMEs and docs.
2. Record changes in [`CHANGELOG.md`](CHANGELOG.md) and summarise in
   `docs/RELEASE_NOTES.md`.
3. Run the full validation matrix (linting, tests, determinism checks, wheel
   build).
4. Tag the repository (`git tag v0.2.0`) and publish artefacts if required.

For historical details and localisation parity refer to the changelog.
