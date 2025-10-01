# ChEMBL Data Acquisition Utilities

The primary documentation and reference material live in the [docs/](docs/) directory.

## Features

* Unified CLI flags such as `--input`, `--output`, `--log-level`, `--sep`, `--encoding`, `--column`, plus `--config` and
  `--print-config` to manage configuration files. Batch size is controlled via `--chunk-size` or `--batch-size`
  depending on the pipeline.
* Streaming CSV handling with deterministic output for large datasets.
* Schema validators in [`schemas/`](schemas/) and dictionaries in [`dictionary/`](dictionary/) that enforce types,
  ranges and reference data.
* Configuration driven by `config/config.yaml`, environment variables and CLI overrides.
* Logging based on the standard `logging` module with configurable verbosity.
* Full static typing (PEP 484), linting with `ruff`, formatting with `black`, type checking via `mypy` and tests with `pytest`.

## Requirements

| Component | Minimum supported | Latest tested |
|-----------|-------------------|---------------|
| Python    | ≥3.11             | 3.12          |
| numpy     | 2.3.3             | 2.3.3         |
| pandas    | 2.3.2             | 2.3.2         |
| requests  | 2.32.5            | 2.32.5        |
| PyYAML    | 6.0.2             | 6.0.2         |

See `requirements-dev.txt` or `pyproject.toml` for the full list. Runtime dependencies follow compatible release ranges so patch
updates within each minor version remain supported. Continuous integration validates both the minimum and latest rows above.

### Runtime environment

* Linux or macOS shell with Bash or PowerShell support (Windows users should rely on WSL2).
* Up-to-date versions of `pip`, `setuptools` and `wheel`; follow the steps in [Installation](#installation).
* Network access to ChEMBL/PubChem/UniProt APIs (port 443).

## Installation

### Runtime environment

* Upgrade packaging tools before installing the project.

  ```bash
  python -m pip install --upgrade pip setuptools wheel
  ```

* Create an isolated virtual environment to keep dependencies local.

  ```bash
  python -m venv .venv
  source .venv/bin/activate  # Windows: .venv\\Scripts\\activate
  ```

### Steps

Clone the repository, change into it and install the package with development extras. Afterwards enable pre-commit hooks so
formatting, linting, type checking and unit tests run automatically.

```bash
git clone https://github.com/<org>/ChEMBL_data_acquisition.git
cd ChEMBL_data_acquisition
pip install .[dev]
pre-commit install
```

Sensitive configuration such as API tokens belongs in a local `.env` file – see [Configuration via `.env`](#configuration-via-env)
for usage guidelines.

## Quick Start

1. **Install dependencies** – follow the steps in [Installation](#installation).
2. **Initialise pre-commit hooks**

   ```bash
   pre-commit install
   ```

   Git hooks ensure formatting, linting, static type checks and tests run before each commit.

3. **Run a sample script**

   ```bash
   python -m library.utils.cli_tools.get_activities --limit 10 --log-level INFO
   ```

   This lightweight helper only emits structured log messages describing the dummy activity rows it would generate; it neither
   reads input files nor writes outputs. Use it to verify logging configuration and CLI wiring before launching full pipelines.
   Common CLI flags include `--limit` to cap processed records, `--log-level` for verbosity, `--sep` for CSV delimiters and
   `--encoding` for file encoding. For end-to-end exports that create files, run one of the data pipelines, for example:

   ```bash
   python -m library.utils.cli_tools.mapper_main --input tests/data/chembl_targets_min.csv \
       --column target_chembl_id --output out/targets_mapped.csv --log-level DEBUG
   python -m library.utils.cli_tools.table_quality_main --input tests/data/chembl_targets_min.csv \
       --output out/quality --table-name chembl_targets --log-level INFO
   ```

   In the second example the `--output` argument must point to a directory where the report files will be created.

4. **Run the tests** – see [Tests](#tests).

## Tests

The `pre-commit` suite runs formatting, linting and static type checks. Execute `pytest` for unit tests and add coverage flags
when required. Determinism and smoke checks are available through dedicated CLI helpers.

```bash
pre-commit run --all-files
pip check
pytest
pytest --cov=library --cov=scripts --cov-report=term-missing --cov-report=xml
python -m library.utils.cli_tools.check_determinism --log-level DEBUG
python -m library.utils.cli_tools.mapper_batch_main --input chembl_ids.csv \
    --output out/mapped.csv --log-level INFO
```

Before running the smoke command, create a `chembl_ids.csv` file with a `chembl_id` header and the required identifiers.

## Data generation

Five production pipelines live in [`scripts/`](scripts/) and write CSV outputs to `data/output/`:

* `get_activity_data.py` — retrieves activity data from ChEMBL and enriches it with derived value ranges.
* `get_assay_data.py` — exports assay descriptions.
* `get_document_data.py` — merges publication metadata from ChEMBL and aggregators (PubMed, Semantic Scholar, OpenAlex, Crossref).
* `get_target_data.py` — collects target information from ChEMBL, UniProt and IUPHAR.
* `get_testitem_data.py` — enriches compounds with structural attributes and PubChem data.

The cached harness `library.utils.cli_tools.pipeline_targets_main`, located under
[`library/utils/cli_tools/`](library/utils/cli_tools/), reuses the CLI contract of the
production target pipeline while operating solely on local files and prepared identifier
batches without network calls.

Example full pipeline execution:

```bash
python -m scripts.get_activity_data --input tests/data/activity_ids_small.csv \
    --output data/output/activities.csv --limit 10 --log-level INFO
```

The command reads data from the ChEMBL API, writes the CSV table and the accompanying `*.meta.yaml`. Development utilities are in
`library/utils/cli_tools/`; for instance, the `get_activities` module focuses on demo logging and performs no file operations. See
[`docs/CLI_TOOLS.md`](docs/CLI_TOOLS.md) for descriptions and command patterns. The output directory is ignored by Git and exposed
as a CI artifact.

> **Note.** The legacy `activity_extraction_main.py` entry point has been superseded by the modular
> `python -m scripts.get_activity_data`, improving maintainability and virtual environment compatibility.

## Usage

The examples below demonstrate how to run the main CLI tools with common options such as `--input`, `--output` and `--limit`.

### `scripts/get_document_data.py`

Retrieve document metadata for a list of PubMed IDs using the bundled sample file:

```bash
python -m scripts.get_document_data pubmed \
    --input tests/data/pmids.csv \
    --output out/documents.csv \
    --limit 5 \
    --log-level INFO
```

The `tests/data/pmids.csv` file contains a small set of PMIDs for experimentation.

You can also run the PubMed pipeline directly via the library module:

```bash
python -m library.pubmed_library \
    --input-csv tests/data/pmids.csv \
    --output out/documents.csv \
    --log-level INFO
```

### `scripts/get_target_data.py`

Fetch basic target information from ChEMBL:

```bash
python -m scripts.get_target_data chembl \
    --input path/to/targets.csv \
    --output out/targets.csv \
    --limit 5 \
    --log-level INFO
```

Replace `path/to/targets.csv` with a CSV containing a `target_chembl_id` column. Input and output use the same column to align with
validation schemas.

### `library.utils.cli_tools.pipeline_targets_main`

Exercise the chunking and batch configuration used by the production target pipeline without contacting remote services:

```bash
python -m library.utils.cli_tools.pipeline_targets_main \
    --input tests/data/chembl_targets_min.csv \
    --output out/targets_cached.csv \
    --chunk-size 25 \
    --batch-size 25 \
    --limit 100
```

The command reads target identifiers from the CSV, chunks them according to `--chunk-size` and `--limit`, forwards the batch size to
`library.pipeline_targets.run_pipeline` and writes the cached ChemBL table with pipeline metadata via `write_csv`. Use it to verify
CLI overrides, logging and deterministic output before launching the network-backed `get_target_data` pipeline.

### `library/utils/cli_tools/get_activities.py`

Generate dummy activity entries without contacting external services:

```bash
python -m library.utils.cli_tools.get_activities --limit 500 --dry-run
```

The command logs that it would generate 500 activity rows and exits without creating any files.

## Updating dependencies

Refresh pinned dependencies periodically and confirm compatibility:

```bash
pip install -U .[dev]
pre-commit run --all-files
```

The first command upgrades runtime and development dependencies to the newest minor releases allowed by the version ranges. The
second command formats code, lints, runs static type checks and executes the test suite to confirm nothing broke during the upgrade.

## Configuration via `.env`

Some utility parameters can be provided through environment variables. To avoid exporting them manually each time, place
`NAME=value` pairs in a `.env` file and load them with [`python-dotenv`](https://pypi.org/project/python-dotenv/).

Example file:

```dotenv
CHEMBL_DA_LOG_LEVEL=INFO
CHEMBL_DA_BASE=https://www.ebi.ac.uk/chembl/api/data
```

Use either the short alias `CHEMBL_DA_BASE` or the fully qualified
`CHEMBL_DA__SOURCES__CHEMBL__API__CHEMBL_BASE`; both expand to the same setting.
See the [alias table](library/config.py#L1471-L1498) in `library/config.py` for
the complete mapping list.

See `.env.example` for typical contact e-mail variables.

Run a script with automatic configuration loading:

```bash
python -m dotenv run -- python -m scripts.get_assay_data --input assay_ids.csv \
    --output out/assays.csv
```

The `assay_ids.csv` file must contain an `assay_chembl_id` column with the required identifiers, for example:

```csv
assay_chembl_id
CHEMBL1234567
CHEMBL2345678
```

Utilities read environment variables automatically, so values from `.env` are available to all CLI tools without extra arguments.

## `api.user_agent`

The `api.user_agent` parameter identifies the application in API requests and must contain **real** contact details. Keeping the
placeholder `contact@example.org` causes configuration validation to fail. Default value:

```yaml
api:
  user_agent: "chembl-da/0.1 (mailto:contact@example.org)"
```

Override the parameter in `config/config.yaml` or via the environment variable `CHEMBL_DA__SOURCES__CHEMBL__API__USER_AGENT`. The
`contact@example.org` value is a placeholder and the validator rejects it, so update the address before running the tools. There is
no dedicated CLI flag (see `library/cli/parser.py`), so configuration is limited to files or environment variables.


## Configuration validation

`library.config.load_config` checks the values inside `config/config.yaml`. An invalid URL raises `ValueError` during loading:

```yaml
api:
  chembl_base: https://
```

```
ValueError: api.chembl_base must be a valid URL
```

Correct configuration must specify the full service URL:

```yaml
api:
  chembl_base: https://www.ebi.ac.uk/chembl/api/data
```

## Configuration errors

Invalid values in `config/config.yaml` raise `ValidationError`. Example:

```yaml
api:
  rps: -1
```

```python
from library.config import load_config
load_config("config/config.yaml")
```

Output:

```
pydantic_core._pydantic_core.ValidationError: 1 validation error for Config
api.rps
  Input should be greater than or equal to 1 [type=greater_than_equal, input_value=-1, input_type=int]
    For further information visit https://errors.pydantic.dev/2.11/v/greater_than_equal
```

Fix the value to a positive number:

```yaml
api:
  rps: 5  # or any value >= 1
```

Allowed ranges are documented in [`config.schema.json`](./config.schema.json), which is exported from the Pydantic models and serves
as a reference artifact with the minimum `1` specified for `api.rps`.

## Logging

CLI helpers configure structured JSON logging via `library.logging_setup.configure_logger`. Use environment variables or CLI flags to
adjust verbosity. The JSON structure is fixed and cannot be customised.

Set the log level via the `--log-level` flag or the `CHEMBL_DA_LOG_LEVEL` environment variable:

```bash
CHEMBL_DA_LOG_LEVEL=DEBUG python -m scripts.get_assay_data --input assay_ids.csv \
    --output out/assays.csv
```

Sample log entry:

```json
{"ts":"2024-05-01T12:00:00Z","level":"INFO","event":"pipeline_start","run_id":"abc123","stage":"pipeline"}
```

Key fields:

* `ts` – UTC timestamp in ISO 8601 format.
* `level` – severity such as `DEBUG`, `INFO`, `WARN` or `ERROR`.
* `event` – short machine-readable event name.
* `run_id` – unique identifier for the current run.
* `stage` – optional pipeline stage.
* `msg` – optional human-readable message.
* Additional keys – event-specific context such as `elapsed`, `url` or `rows`.

Dry-run executions emit logs with `event` set to `dry_run`, enabling easy filtering, for example:

```bash
jq 'select(.event=="dry_run")' log.jsonl
```

The run identifier is generated by the CLI helpers using `uuid.uuid4().hex` and passed to the logger, which includes it with every
record. Override the value before calling `configure_logger` to inject a custom identifier.

Secrets are automatically redacted: values for keys ending in `token`, `key`, `secret` or `password` are replaced with `"***"`.
Log-level filtering drops records below the configured `--log-level` or `CHEMBL_DA_LOG_LEVEL` setting.

Typical log entries look like:

```json
{"ts":"2024-05-01T12:00:00Z","level":"INFO","event":"pipeline_start","run_id":"abc123","stage":"pipeline"}
{"ts":"2024-05-01T12:00:01Z","level":"INFO","event":"request_ok","run_id":"abc123","stage":"fetch","url":"https://api.example.org","status":200}
{"ts":"2024-05-01T12:00:02Z","level":"INFO","event":"validate_done","run_id":"abc123","stage":"validate","rows":42}
{"ts":"2024-05-01T12:00:03Z","level":"INFO","event":"pipeline_done","run_id":"abc123","stage":"pipeline","elapsed":3.2}
```

## Reproducibility

Deterministic CSV writers in `library.io` keep outputs and metadata stable across runs. The function
`library.csv_utils.write_csv_deterministic` normalises column ordering and metadata so repeated executions yield identical files.
All CLI scripts share a common set of flags:

```bash
python -m library.utils.cli_tools.table_quality_main --input data.csv --table-name data
```

`--output` defaults to `output.<input_name>_YYYYMMDD.csv` in the directory defined by `local.io.output_dir`. For additional
examples see [`docs/USAGE_EN.md`](docs/USAGE_EN.md) (Russian version: [`docs/USAGE_RU.md`](docs/USAGE_RU.md)).

## Project structure

```
ChEMBL_data_acquisition/
├── config/config.yaml
├── dictionary/
├── library/
│   ├── __init__.py
│   ├── chembl_client.py
│   ├── csv_utils.py
│   ├── config.py
│   └── ...
├── scripts/
│   ├── get_activity_data.py
│   ├── get_assay_data.py
│   ├── ...
├── tests/
│   └── data/
└── docs/
    ├── CONFIG_EN.md / CONFIG_RU.md
    ├── OUTPUT_EN.md / OUTPUT_RU.md
    └── USAGE_EN.md / USAGE_RU.md
```

## Configuration

Parameters are read from `config/config.yaml`, environment variables (`CHEMBL_DA__...`) and CLI flags. Details are documented in
[`docs/CONFIG_EN.md`](docs/CONFIG_EN.md) (Russian version: [`docs/CONFIG_RU.md`](docs/CONFIG_RU.md)).

### Environment variables

Environment variables override values from the YAML file. Names follow the `CHEMBL_DA__...` pattern with double underscores
separating each nested section. For example, to enable debug logging:

```bash
export CHEMBL_DA__LOG__LEVEL=DEBUG
```

Most options also provide short aliases for backwards compatibility. The table lists every supported alias and the canonical key it
maps to:

| Alias | Equivalent key |
|-------|----------------|
| `CHEMBL_DA_BASE` | `CHEMBL_DA__SOURCES__CHEMBL__API__CHEMBL_BASE` |
| `CHEMBL_DA_BURST` | `CHEMBL_DA__SOURCES__CHEMBL__API__BURST` |
| `CHEMBL_DA_CACHE_DIR` | `CHEMBL_DA__LOCAL__IO__CACHE_DIR` |
| `CHEMBL_DA_CACHE_MAXSIZE` | `CHEMBL_DA__SOURCES__CHEMBL__CACHE__CACHE_MAXSIZE` |
| `CHEMBL_DA_CACHE_TTL` | `CHEMBL_DA__SOURCES__CHEMBL__CACHE__CACHE_TTL` |
| `CHEMBL_DA_CROSSREF_BASE` | `CHEMBL_DA__SOURCES__CROSSREF__BASE` |
| `CHEMBL_DA_CROSSREF_MAILTO` | `CHEMBL_DA__SOURCES__CROSSREF__MAILTO` |
| `CHEMBL_DA_CROSSREF_TIMEOUT_CONNECT` | `CHEMBL_DA__SOURCES__CROSSREF__TIMEOUT_CONNECT` |
| `CHEMBL_DA_CROSSREF_TIMEOUT_READ` | `CHEMBL_DA__SOURCES__CROSSREF__TIMEOUT_READ` |
| `CHEMBL_DA_DICT_DIR` | `CHEMBL_DA__LOCAL__RESOURCES__DICTIONARY_DIR` |
| `CHEMBL_DA_GLOBAL_BURST` | `CHEMBL_DA__SYSTEM__RATE__GLOBAL_BURST` |
| `CHEMBL_DA_GLOBAL_RPS` | `CHEMBL_DA__SYSTEM__RATE__GLOBAL_RPS` |
| `CHEMBL_DA_IUPHAR_BASE` | `CHEMBL_DA__SOURCES__IUPHAR__BASE` |
| `CHEMBL_DA_IUPHAR_FAMILY_CSV` | `CHEMBL_DA__LOCAL__RESOURCES__IUPHAR_FAMILY_CSV` |
| `CHEMBL_DA_IUPHAR_TARGET_CSV` | `CHEMBL_DA__LOCAL__RESOURCES__IUPHAR_TARGET_CSV` |
| `CHEMBL_DA_LIMITER_CACHE_MAXSIZE` | `CHEMBL_DA__SYSTEM__RATE__LIMITER_CACHE_MAXSIZE` |
| `CHEMBL_DA_LIMITER_CACHE_TTL` | `CHEMBL_DA__SYSTEM__RATE__LIMITER_CACHE_TTL` |
| `CHEMBL_DA_LOG_LEVEL` | `CHEMBL_DA__SYSTEM__LOG__LEVEL` |
| `CHEMBL_DA_MOLECULE_CATALOG_CACHE` | `CHEMBL_DA__SOURCES__CHEMBL__MOLECULE_CATALOG__CACHE_PATH` |
| `CHEMBL_DA_OPENALEX_BASE` | `CHEMBL_DA__SOURCES__OPENALEX__BASE` |
| `CHEMBL_DA_OPENALEX_MAILTO` | `CHEMBL_DA__SOURCES__OPENALEX__MAILTO` |
| `CHEMBL_DA_OPENALEX_TIMEOUT_CONNECT` | `CHEMBL_DA__SOURCES__OPENALEX__TIMEOUT_CONNECT` |
| `CHEMBL_DA_OPENALEX_TIMEOUT_READ` | `CHEMBL_DA__SOURCES__OPENALEX__TIMEOUT_READ` |
| `CHEMBL_DA_OUTDIR` | `CHEMBL_DA__LOCAL__IO__OUTPUT_DIR` |
| `CHEMBL_DA_PUBCHEM_BASE` | `CHEMBL_DA__SOURCES__PUBCHEM__BASE` |
| `CHEMBL_DA_RETRY_BACKOFF_FACTOR` | `CHEMBL_DA__SYSTEM__RETRY__BACKOFF_FACTOR` |
| `CHEMBL_DA_RETRY_MAX_ATTEMPTS` | `CHEMBL_DA__SYSTEM__RETRY__MAX_ATTEMPTS` |
| `CHEMBL_DA_RPS` | `CHEMBL_DA__SOURCES__CHEMBL__API__RPS` |
| `CHEMBL_DA_TARGETS_TYPE_CSV` | `CHEMBL_DA__LOCAL__RESOURCES__TARGETS_TYPE_CSV` |
| `CHEMBL_DA_TIMEOUT_CONNECT` | `CHEMBL_DA__SOURCES__CHEMBL__API__TIMEOUT_CONNECT` |
| `CHEMBL_DA_TIMEOUT_READ` | `CHEMBL_DA__SOURCES__CHEMBL__API__TIMEOUT_READ` |
| `CHEMBL_DA_UNIPROT_BASE` | `CHEMBL_DA__SOURCES__UNIPROT__API__BASE` |
| `CHEMBL_DA_UNIPROT_DATA_DIR` | `CHEMBL_DA__LOCAL__RESOURCES__UNIPROT_DATA_DIR` |
| `CHEMBL_DA__IO__CACHE_DIR` | `CHEMBL_DA__LOCAL__IO__CACHE_DIR` |
| `CHEMBL_DA__IO__EXIST_OK` | `CHEMBL_DA__LOCAL__IO__EXIST_OK` |

See [`docs/CONFIG_EN.md`](docs/CONFIG_EN.md) for a complete overview of all configuration options (Russian version —
[`docs/CONFIG_RU.md`](docs/CONFIG_RU.md)).

### Schema validation

Configuration values are validated by :func:`library.config.load_config`, which calls :meth:`Config.model_validate <pydantic.BaseModel.model_validate>`.
Validation follows the model definitions and produces detailed error messages for nested fields. The accompanying `config.schema.json`
is generated from the same Pydantic model for IDE hints and documentation; it is not executed at runtime and should not be passed to
`jsonschema`.

Command-line flags have the highest priority. All utilities accept `--config` to point at a configuration file and `--print-config`
to show the effective values after all overrides have been applied. The final precedence is:

```
YAML < environment variables < CLI options
```

Only the top-level command line scripts read the configuration file. Modules under `library/` expect a `Config` (or one of its
subsections) to be passed explicitly, making dependencies clear and avoiding hidden global state.

The directories referenced by `local.io.output_dir` and `local.io.cache_dir` are checked but not created when loading the
configuration. Scripts that need these paths can call `library.config.ensure_dirs` after `load_config` to create them if they are
missing and `local.io.exist_ok` permits it.

Path values such as `local.io.output_dir`, `local.io.cache_dir` and the `local.init` workbook paths are exposed as `pathlib.Path`
objects. String values in `config/config.yaml` or overrides from the environment and command line are automatically converted.

```bash
python -m library.utils.cli_tools.table_quality_main \
    --input tests/data/activity.csv \
    --table-name activity
```

`--output` defaults to `output.<input_name>_YYYYMMDD.csv` in the directory specified by `local.io.output_dir`.
For additional examples see [`docs/USAGE_EN.md`](docs/USAGE_EN.md) (Russian version: [`docs/USAGE_RU.md`](docs/USAGE_RU.md)).

## Output and metadata

Pipelines persist deterministic CSV tables via `library.io.write_csv` and store accompanying `*.meta.yaml` sidecars in
`data/output`. Each sidecar records the Git commit, launch parameters, SHA-256 checksum and row/column statistics. See
[`docs/OUTPUT_EN.md`](docs/OUTPUT_EN.md) / [`docs/OUTPUT_RU.md`](docs/OUTPUT_RU.md) for layout details.

## Dtype Inspector

The `dtype_inspector` utility executes each `get_*_data` helper on a small identifier set and logs the resulting pandas dtypes. Run
it periodically to spot schema drift when upstream services change responses.

```bash
python -m library.utils.cli_tools.dtype_inspector_main --log-level INFO
```

Consider wiring the script into CI to detect dtype changes early. The tool is lightweight and makes only a handful of requests per
dataset.

## Development and testing

Individual tools such as `ruff`, `black` and `mypy` are wired into `pre-commit` but can be executed manually when iterating on
specific modules.

```bash
ruff check scripts library library/utils/cli_tools
black scripts library library/utils/cli_tools
mypy scripts library
pytest
```

Test datasets live in `tests/data`; `library.utils.cli_tools.check_determinism` verifies repeatable CSV output.

## Licence

MIT License. See the `LICENSE` file (if present).

Add new reference materials directly to the `docs` directory when updating documentation.
