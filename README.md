# ChEMBL Data Acquisition Utilities

| Language | Overview |
|----------|----------|
| English  | [`docs/README_EN.md`](./docs/README_EN.md) |
| Русский  | [`docs/README_RU.md`](./docs/README_RU.md) |

Use the [documentation index](./docs/index.md) to navigate the full manual set.
All content is maintained in synchronised English (`*_EN.md`) and Russian
(`*_RU.md`) variants.

## Feature highlights

- Pipelines for documents, targets, assays, activities, and test items with
  deterministic exports and YAML sidecars.
- Orchestrator (`get-data`) that chains all pipelines with shared configuration
  and logging.
- Unified CLI layer with shared flags (`--input`, `--final-out`, `--config`,
  `--log-level`) plus staging switches (`--raw-out`, `--raw-format`,
  `--id-cols`, `--no-reindex-raw`, `--normalize-at-export`) for the target
  pipeline.
- Configuration via YAML (`config/config.yaml`), `config.local.yaml`,
  environment variables (`CHEMBL_DA__*`) and CLI overrides. Optional `.env`
  files can be sourced locally.
- Quality gates: schema validation, table-quality reports, deterministic CSV
  writing, strict typing, linting, and unit tests.

## Console entry points

Installing the package (`pip install .` or the published wheel) registers the
following scripts. They mirror the modules in `scripts/` or
`library/cli/commands/`.

| Script | Module | Purpose |
|--------|--------|---------|
| `get-data` | `scripts.get_data:main` | Run all pipelines sequentially. |
| `get-document-data` | `library.cli.commands.get_document_data:main` | Document acquisition and enrichment. |
| `get-target-data` | `library.cli.commands.get_target_data:main` | Target aggregation (ChEMBL, UniProt, IUPHAR). |
| `get-assay-data` | `library.cli.commands.get_assay_data:main` | Assay metadata export. |
| `get-activity-data` | `library.cli.commands.get_activity_data:main` | Activity export with normalisation. |
| `get-testitem-data` | `library.cli.commands.get_testitem_data:main` | Molecule enrichment with PubChem. |
| `get-document-type` | `library.cli.commands.get_document_type:main` | Publication classification helper. |
| `csv-utils` | `library.cli.commands.csv_utils:main` | Deterministic CSV utilities. |
| `mapper` | `library.cli.commands.mapper:main` | UniProt/ChEMBL mapping tool. |
| `table-quality` | `library.cli.commands.table_quality:main` | Column-level quality reports. |
| `chunk-io` | `library.cli.commands.chunk_io:main` | Chunked CSV IO harness. |
| `get-input-initialisation` | `library.cli.commands.get_input_initialisation:main` | Merge Excel initialisation workbooks. |
| `get-activities` | `library.cli.commands.get_activities:main` | Synthetic activity generator for smoke tests. |
| `check-determinism` | `library.cli.commands.check_determinism:main` | Compare CSV hashes across runs. |

Helper utilities under `library.utils.cli_tools` are documented in
[`docs/CLI_TOOLS_EN.md`](./docs/CLI_TOOLS_EN.md) /
[`docs/CLI_TOOLS_RU.md`](./docs/CLI_TOOLS_RU.md).

## Requirements

| Component | Minimum | Latest tested |
|-----------|---------|---------------|
| Python | 3.11 | 3.12 |
| numpy | 1.26 | 2.3.3 |
| pandas | 2.0 | 2.3.3 |
| requests | 2.31 | 2.32.5 |
| PyYAML | 6.0 | 6.0.3 |

Follow the installation, configuration, and QA guidance in the dedicated docs:

- [`docs/README_EN.md`](./docs/README_EN.md) / [`docs/README_RU.md`](./docs/README_RU.md) – project overview.
- [`docs/USAGE_EN.md`](./docs/USAGE_EN.md) / [`docs/USAGE_RU.md`](./docs/USAGE_RU.md) – CLI reference and examples.
- [`docs/CONFIG_EN.md`](./docs/CONFIG_EN.md) / [`docs/CONFIG_RU.md`](./docs/CONFIG_RU.md) – configuration matrix.
- [`docs/OUTPUT_EN.md`](./docs/OUTPUT_EN.md) / [`docs/OUTPUT_RU.md`](./docs/OUTPUT_RU.md) – exported artefacts.
- [`docs/QA_PROCESS_EN.md`](./docs/QA_PROCESS_EN.md) / [`docs/QA_PROCESS_RU.md`](./docs/QA_PROCESS_RU.md) – QA playbook.

## Development quick start

```bash
python -m pip install --upgrade pip setuptools wheel
python -m venv .venv
source .venv/bin/activate  # Windows: .venv\Scripts\activate
pip install -r requirements-lock.txt
pre-commit install
pytest
```

`requirements-lock.txt` pins the dependency set used in CI. Regenerate it after
changing `pyproject.toml` by installing the project in a fresh virtual
environment and running `pip freeze > requirements-lock.txt`.
