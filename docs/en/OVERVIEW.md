# ChEMBL Data Acquisition Toolkit

This guide gives a project-wide overview of the utilities that download,
normalise, and export ChEMBL-derived datasets. Use it as the canonical landing
page for engineers and analysts before diving into the detailed manuals.

## Highlights

- **Entity pipelines** for documents, targets, assays, activities, and test
  items. Each pipeline streams identifiers from CSV, calls external APIs, runs
  deterministic validation, and writes audited exports.
- **Orchestrator** (`get-data`) that chains all entity pipelines with shared
  configuration and consistent logging.
- **Unified CLI layer** with shared options (`--input`, `--final-out`,
  `--log-level`, `--config`, etc.), optional raw snapshots
  (`--raw-out`, `--raw-format`, `--id-cols`, `--no-reindex-raw`) for the target
  pipeline, and console entry points declared in `pyproject.toml`.
- **Configuration triad** – values come from `config/config.yaml`,
  environment variables (prefixed with `CHEMBL_DA__`) and command line flags.
  Optional `.env` files can be sourced with `python-dotenv` during local runs.
- **Deterministic outputs** – CSV writers fix row/column ordering, attach YAML
  sidecars with provenance, compute SHA-256 hashes, and emit table-quality
  reports for every run.
- **Developer tooling** – strict typing (`mypy`), linting (`ruff`), formatting
  (`black`), testing (`pytest`), coverage, and determinism checks in CI.

## Repository layout

| Path | Purpose |
|------|---------|
| `scripts/` | CLI wrappers for each entity pipeline and the `get-data` orchestrator. |
| `library/` | Core modules: HTTP clients, pipeline orchestration, normalisation, validation, post-processing, IO, and metadata helpers. |
| `library/cli/commands/` | Console-script entry points used by installed wheels. |
| `library/utils/cli_tools/` | Lightweight utilities (table profiling, cached target harness, CSV helpers, mapping tools). |
| `config/` | Default YAML configuration, schema definition, and bundled dictionaries. |
| `config/dictionary/` | Reference datasets used by the pipelines (UniProt caches, target taxonomies, QA fixtures). |
| `data/` | Smoke-test inputs and sample exports. |
| `docs/` | Project documentation (English and Russian variants). |
| `tests/` | Unit and integration tests covering pipelines and CLI helpers. |
| `Makefile` | Convenience targets for formatting, tests, packaging, and linting. |

## Supported entry points

Install the project (`pip install .` or `pip install dist/*.whl`) to register the
following console scripts. They correspond to the modules inside
`scripts/`, `library/cli/commands/`, or `library/utils/cli_tools/` and accept the
same arguments as their `python -m …` equivalents.

| Console script | Module | Description |
|----------------|--------|-------------|
| `get-data` | `scripts.get_data:main` | Run all pipelines sequentially with shared configuration. |
| `get-document-data` | `library.cli.commands.get_document_data:main` | Acquire and enrich ChEMBL documents. |
| `get-target-data` | `library.cli.commands.get_target_data:main` | Aggregate ChEMBL, UniProt, and IUPHAR targets. |
| `get-assay-data` | `library.cli.commands.get_assay_data:main` | Export assay metadata. |
| `get-activity-data` | `library.cli.commands.get_activity_data:main` | Export normalised activity records. |
| `get-testitem-data` | `library.cli.commands.get_testitem_data:main` | Enrich molecule records with PubChem details. |
| `get-document-type` | `library.utils.cli_tools.get_document_type:main` | Classify publication types for QA tasks. |
| `csv-utils` | `library.utils.cli_tools.csv_utils_main:main` | Deterministic CSV re-serialisation helpers. |
| `mapper` | `library.utils.cli_tools.mapper_main:main` | Interactive UniProt/ChEMBL mapper. |
| `table-quality` | `library.utils.cli_tools.table_quality_main:main` | Generate column-level quality profiles. |
| `chunk-io` | `library.utils.cli_tools.chunk_io_main:main` | Stream CSV chunks while keeping ordering stable. |
| `get-input-initialisation` | `library.utils.cli_tools.get_input_initialisation:main` | Merge Excel initialisation workbooks. |
| `get-activities` | `library.utils.cli_tools.get_activities:main` | Emit synthetic activity rows and deterministic CSV + `.meta.yaml` artefacts for smoke tests. |
| `check-determinism` | `library.utils.cli_tools.check_determinism:main` | Compare CSV hashes across runs. |

The dedicated reference [`USAGE.md`](./USAGE.md) (and
[`../ru/USAGE.md`](../ru/USAGE.md)) covers arguments, sub-commands, and advanced
usage scenarios in depth.

## Requirements and installation

| Component | Supported range | Latest tested |
|-----------|-----------------|---------------|
| Python | 3.11.x-3.12.x | 3.12.3 |
| numpy | >=2.3.3,<3.0 | 2.3.3 |
| pandas | >=2.3.3,<3.0 | 2.3.3 |
| requests | >=2.32.5,<3.0 | 2.32.5 |
| PyYAML | >=6.0.3,<7.0 | 6.0.3 |
| openpyxl | >=3.1.5,<4.0 | 3.1.5 |
| pyarrow | >=17.0.0,<18.0 | 17.0.0 |
| jsonschema | >=4.25.1,<5.0 | 4.25.1 |
| pandera | >=0.26.1,<0.27 | 0.26.1 |
| pydantic | >=2.11.9,<3.0 | 2.11.9 |

Pinned dependencies live in `requirements-lock.txt`. Regenerate the file only
after editing `pyproject.toml` ranges, using a clean virtual environment and
`pip freeze`.

```bash
python -m pip install --upgrade pip setuptools wheel
python -m venv .venv
source .venv/bin/activate  # Windows: .venv\Scripts\activate
git clone https://github.com/SatoryKono/ChEMBL_data_acquisition.git
cd ChEMBL_data_acquisition
pip install -r requirements-lock.txt
pre-commit install
```

Wheel users can install with `pip install chembl-data-acquisition` once a release
is published; console scripts and the packaged configuration are placed in the
platform-specific user directories listed in [`CONFIG.md`](./CONFIG.md).

## Quick start

1. **Prepare identifier lists.** Use the templates in `data/input` (for example
   `document.csv`, `target.csv`, `assay.csv`, `activity.csv`, `testitem.csv`) or
   export fresh ID lists from your warehouse. Additional, larger smoke samples
   live in `data/input/full`. Each pipeline expects one identifier per row; see
   [`DATA_SCHEMA.md`](./DATA_SCHEMA.md) for column names.
2. **Review configuration.** Copy `config/config.yaml` if you need to override
   API limits, output directories, or staging flags. Environment variable
   overrides follow the `CHEMBL_DA__SECTION__KEY` pattern. Details live in
   [`CONFIG.md`](./CONFIG.md).
3. **Run a pipeline.**

   ```bash
   get-document-data all \
       --input data/input/document.csv \
       --final-out output.documents_$(date +%Y%m%d).csv \
       --config config/config.yaml \
       --log-level INFO
   ```

   Add `--limit 10` for smoke tests. The target pipeline exposes staging
   switches such as `--raw-out` and `--raw-format parquet` for raw snapshots.
4. **Inspect artefacts.** Every CSV is accompanied by `<name>.meta.yaml`,
   quality reports, and (for documents) JSON summaries. See
   [`OUTPUT.md`](./OUTPUT.md) for formats and retention guidance.

## Testing and quality gates

Run the following commands before committing changes or publishing artefacts:

```bash
pre-commit run --all-files
pip check
pytest
pytest --cov=library --cov=scripts --cov-report=term-missing
python -m library.utils.cli_tools.check_determinism --log-level DEBUG \
    --input out/latest.csv --previous out/previous.csv
```

The QA playbook (`QA_PROCESS.md` / `../ru/QA_PROCESS.md`) documents
release gates, smoke checks, and acceptance criteria. Determinism checks rely on
YAML sidecars, so keep them under version control when comparing runs.

## Related documentation

- [`USAGE.md`](./USAGE.md) / [`../ru/USAGE.md`](../ru/USAGE.md) – CLI
  options, sub-commands, and execution recipes.
- [`CONFIG.md`](./CONFIG.md) / [`../ru/CONFIG.md`](../ru/CONFIG.md) –
  configuration sources, environment variables, and staged output locations.
- [`OUTPUT.md`](./OUTPUT.md) / [`../ru/OUTPUT.md`](../ru/OUTPUT.md) –
  artefact layout, metadata sidecars, and raw snapshot handling.
- [`DATA_SCHEMA.md`](./DATA_SCHEMA.md) /
  [`../ru/DATA_SCHEMA.md`](../ru/DATA_SCHEMA.md) – column definitions and
  validation schemas.
- [`ETL_PROCESS.md`](./ETL_PROCESS.md) /
  [`../ru/ETL_PROCESS.md`](../ru/ETL_PROCESS.md) – end-to-end data flow.
- [`CLI_TOOLS.md`](./CLI_TOOLS.md) /
  [`../ru/CLI_TOOLS.md`](../ru/CLI_TOOLS.md) – helper utilities for QA and
  diagnostics.

For architecture diagrams consult [`architecture/ARCHITECTURE.md`](./architecture/ARCHITECTURE.md)
and its Russian counterpart (`../ru/architecture/ARCHITECTURE.md`).
