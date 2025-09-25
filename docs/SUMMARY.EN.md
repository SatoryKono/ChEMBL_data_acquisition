# Project Overview
- ChEMBL Data Acquisition is a suite of Python 3.12 CLI scripts and libraries for downloading, normalizing, and exporting biodata from the ChEMBL, PubChem, UniProt, PubMed, and related APIs with deterministic CSV/Parquet outputs.


- The project ships unified launch flags, streaming CSV handling, validation schemas, reference lookups, and an execution metadata/logging stack.


- Documentation covers installation, pre-commit usage, smoke commands, and tests to guarantee a reproducible pipeline for new data team members.


- Configuration is managed via a YAML file with overrides through environment variables and CLI flags; the layout is documented in CONFIG.md.


- Output tables are accompanied by sidecar files with hashes, configuration, and statistics, simplifying data quality audit and monitoring.



- The codebase emphasizes strict typing, pandera-based schemas, and a comprehensive QA toolchain (black, ruff, mypy, pytest, deterministic CSV).



- Current docs omit a roadmap or known issues — the README focuses on execution, testing, and data generation without referencing future work.



## Architecture and Modules
- **scripts/** — CLI wrappers for individual pipelines (activities, assays, targets, etc.) handling errors, normalization, pandera validation, and sidecar generation.


- **library/** — core logic:
  - `io.py` — reads identifiers/CSVs and writes deterministic CSVs while creating metadata and directory structures automatically.


  - `chembl_client.py` — HTTP client with caching, RPS throttling, and exponential backoff for the ChEMBL API.


  - `document_pipeline.py`, `target_postprocessing.py`, `input_initialisation_library.py` — aggregation and post-processing modules mirroring the production Power Query flow with metric calculations, dictionary mapping, and filtering statistics.



  - `cli/` — centralized parser assembly, configuration loading, logger setup, and shared arguments (input/output/encoding/chunk-size, etc.).


  - `logging_setup.py`, `log.py` — structured JSON logging with secret redaction and a global logger tracking status and RPS.



  - `metadata.py`, `sidecar.py`, `pipeline_metadata.py` — generate sidecar YAML, compute SHA-256, and append pipeline version and timestamp fields.




- **schemas/** — pandera schemas and normalizers (activities, targets, documents, testitems) enforcing value and type checks.


- **docs/** — reference guides for execution, configuration, and output artifacts.



- **tests/** — extensive coverage across CLI, config, validation, post-processing, and determinism; includes activity smoke tests and target pipeline/config checks.




## Data Flow (Diagram)
```
[Identifier CSV]
        │ read_ids()
        ▼
[ID iterator] --limit--> [ID batches]
        │ ChemblClient/request_json() + rate limiting
        ▼
[Raw API payloads] --normalize_*/postprocess--> [DataFrame]
        │ Pandera validation + SidecarErrors
        ▼
[Cleaning and aggregation]
        │ add_pipeline_metadata()
        ▼
write_csv() ──► <table>.csv + <table>.csv.meta.yaml + (opt.) failure_cases.csv
```
- Identifier ingestion with empty-value handling and column checks is performed via `read_ids`.


- API retrieval uses `ChemblClient` with cache, RPS throttling, and exponential backoff.


- Normalization, validation, and sidecar error capture live in CLI scripts leveraging pandera and `SidecarErrors`.



- CSV writing and metadata serialization rely on `write_csv`, `write_meta_yaml`, and injecting version/timestamp columns.




## Configuration
- The primary `config.yaml` contains sections such as `api`, `chembl`, `openalex`, `crossref`, `uniprot`, `pubchem`, `document`, `target`, `resources`, `io`, `jobs`, `batch`, etc., with defaults for URLs, timeouts, RPS limits, chunk sizes, and dictionary paths.


- Pydantic models enforce typing, URL validation, and a mandatory `user_agent` with an email; precedence order: YAML < environment < CLI.


- Documentation explains overrides through `CHEMBL_DA__SECTION__KEY` variables, short aliases, and `--section.key=value` flags.


- Tests confirm that empty YAML applies default RPS settings, environment variables and aliases override values correctly, and missing files raise `ConfigError`.



## Dictionaries and Reference Data
- Configuration points to the root `dictionary` directory, IUPHAR and UniProt references, plus `targets_type.csv` and related CSVs for target classification.


- `input_initialisation_library` consumes `dictionary/_Curation/citation_fraction.csv` and `dictionary/_Target/targets_type.csv` to compute citation metrics, target types, and other attributes; missing files trigger explicit errors with expected paths.


- Enrichment stages read local UniProt JSON and classification CSVs to reconcile external references with the ChEMBL export.




## Dependencies and Installation
- Minimum versions: Python 3.12, pandas 2.1, requests 2.31, PyYAML 6.0; the complete list of runtime and dev tools (black, ruff, mypy, pytest, responses, hypothesis, etc.) is declared in pyproject.toml.



- Install via `pip install .[dev]` or `pip install -e .[dev]`; use a virtual environment and upgrade pip/setuptools/wheel.



- Enable `pre-commit install` afterward to run formatting, linting, and static checks automatically before commits.



## How to Run (Steps)
1. Clone the repository, create a Python 3.12 virtual environment, and activate it.


2. Install project and dev dependencies with `pip install .[dev]`.


3. Run `pre-commit install` and optionally `pre-commit run --all-files` for an initial quality sweep.


4. Prepare an input CSV with the required identifier column (defaults: `input.csv`/`activity_id`).


5. Execute the desired CLI script, e.g. `python -m scripts.get_activities --input tests/data/activity_ids_small.csv --output out/activities.csv --limit 10 --log-level INFO`.


6. Alternative pipelines: `get_assay_data`, `get_target_data`, `get_document_data`, `get_testitem_data`, `get_input_initialisation`, `table_quality_main`.


7. Inspect the generated CSV and `.meta.yaml` files under `data/output` (or the configured directory).



## Tests and Verification
- Core quality checks: `pre-commit run --all-files`, `pytest`, `scripts/check_determinism.py --log-level DEBUG`; these commands cover formatting, linting, typing, unit tests, and deterministic output.


- CLI smoke tests validate successful CSV generation with API stubs and file hashing.


- Config tests cover aliasing, environment overrides, and missing-file scenarios, helping diagnose setup issues before pipeline execution.


- Pipeline-oriented tests verify batch size propagation, optional-stage handling, and DataFrame column expectations.


- Deterministic CSV writing is enforced via a standalone script and unit tests comparing SHA-256 hashes across reruns.




## Output Artifacts (Tables/Files, Fields, Examples)
- Primary results are CSV files under `io.output_dir` (`data/output` by default) with optional subfolders reflecting source or processing stage.


- Each CSV is paired with `<name>.csv.meta.yaml` containing the git SHA, command line, configuration, row/column statistics, and SHA-256 hash of the main file.


- Validation failures yield dedicated `*_failure_cases.csv` via `SidecarErrors`, easing debugging of problematic records.



- `table_quality_main.py` and `get_input_initialisation.py` also emit quality reports (`<table>_quality_report_table.csv`, `<table>_data_correlation_report_table.csv`).


- Data outputs automatically receive `pipeline_version` and `timestamp_utc` columns for release traceability.




## Logging and Diagnostics
- Logging relies on a custom JSON logger with secret redaction and global context; level and run_id are configurable via CLI/config.




- CLI scripts log key events: input reads, limits, requests, API errors, validation outcomes, and file writes, streamlining incident investigation.


- Metadata is emitted into sidecar YAML alongside hashes and configuration, providing an additional diagnostic trail for run history.




## Constraints and Assumptions
- Requires Python 3.12, outbound network access (port 443), and a valid API `user_agent` containing an email; otherwise Pydantic validation fails.



- Default RPS, timeout, and backoff values may need adjustment under stricter API limits; throttling is enforced via `ChemblClient` and its rate limiter.



- Reference data must follow the `dictionary/…` layout; missing key CSVs raise exceptions, so synchronize dictionaries beforehand.



- README does not list roadmap or known issues — clarify additional constraints or plans separately.



## Common Errors and Fixes
- **Missing input column/file** — `read_ids` and CLI emit `read_fail`; check the column name (`--column`) and CSV path.


- **Invalid limits/parameters** — negative `activity.limit` logs `invalid_limit`; update the value in `config.yaml` or CLI.


- **API failure** — `ChemblClient` logs `request_fail` and retries; persistent issues require network, timeout, and user-agent checks.


- **Missing dictionaries** — absence of `targets_type.csv` or `citation_fraction.csv` raises FileNotFoundError with path hints; sync `dictionary` or pass `--dictionary`.


- **Missing config** — loading a nonexistent YAML raises `ConfigError`; create the file or rely on the default path.


- **Data validation issues** — schema mismatches create `*_failure_cases.csv`; review the file, adjust dictionaries or cleaning, and rerun.



## Glossary (ChEMBL/ETL)
- **ActivitiesSchema** — pandera schema for activity tables (columns, range checks).


- **Sidecar** — auxiliary CSV with validation errors or YAML with run metadata, generated automatically during output writes.



- **Pipeline metadata** — `pipeline_version` and `timestamp_utc` columns appended to all outputs for version/time tracking.


- **Deterministic CSV** — approach combining row/column ordering and SHA-256 checks to ensure identical outputs across reruns.



- **Input initialisation** — process combining raw Excel/CSV inputs, filtering, and distributing across entities (activity, assay, target, document, testitem) with metric calculations.



## Quality Requirements
- Run formatting, linting, static typing, and pytest through pre-commit or directly before release.


- Data must pass pandera validation; recorded errors in sidecars are blocking until resolved.


- CSV determinism must be preserved (validated via unit tests and the standalone script).



- Each output must include metadata with hashes and configuration; missing sidecar files violate reproducibility requirements.




## Launch Checklist
1. Python 3.12 installed; outbound port 443 confirmed.


2. Virtual environment created/activated; pip/setuptools/wheel upgraded.


3. `pip install .[dev]` completed without dependency errors.


4. `pre-commit install` executed; optionally `pre-commit run --all-files`.


5. Input CSV prepared with the expected identifier column (or supplied via `--column`).


6. Confirmed presence of `dictionary/…` assets or supplied `--dictionary` path.



7. Target command (`python -m scripts.get_activities ...`, etc.) executed with required flags and log level.



8. Verified CSV and `.meta.yaml` creation under `data/output` (or configured directory).


9. Reviewed JSON logs for errors and presence of run_id.


10. If needed, ran `pytest` and `scripts/check_determinism.py --log-level DEBUG`.


11. Investigated any `*_failure_cases.csv` to remediate data issues.


12. Documented CLI/env configuration overrides alongside run artifacts.




## One-liner Setup & Run
```
python -m venv .venv && source .venv/bin/activate && pip install .[dev] && python -m scripts.get_activities --input tests/data/activity_ids_small.csv --output data/output/activities.csv --limit 10 --log-level INFO
```



