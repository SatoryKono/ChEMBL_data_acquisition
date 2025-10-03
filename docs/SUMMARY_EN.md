# Project Summary

This document provides an at-a-glance description of the ChEMBL data acquisition
utilities, their structure, shared services and supporting workflows.

## Repository layout

* `scripts/` – command-line entry points for the activity, assay, document,
  target and test-item pipelines.
* `library/utils/cli_tools/` – cached target harness helpers for offline smoke
  tests.
* `library/` – reusable modules covering API clients, rate limiting,
  normalisation, enrichment, validation, deterministic I/O, logging and
  metadata helpers.
* `library/schemas/` – `pandera` validation schemas and normalisers that keep column
  ordering, data types and canonical values consistent across exports.
* `dictionary/` & `data/` – local lookup tables, cached API responses and input
  workbooks consumed by the pipelines.
* `docs/` – the reference documentation you are reading.
* `tests/` – unit and integration coverage for configuration overrides,
  enrichment logic, deterministic exports and CLI wiring.

## Data flow overview

```
[Identifier CSV] --read_ids--> [Iterator] --chunk/limit--> [Batches]
        │                                       │
        ▼                                       ▼
  ChemblClient / external services      Normalise & enrich
        │                                       │
        ▼                                       ▼
    Pandera validation  ──►  SidecarErrors (optional failure CSV)
        │
        ▼
add_pipeline_metadata → write_csv_deterministic →
  <table>.csv + <table>.csv.meta.yaml + quality reports
```

* `io.read_ids` streams identifiers while filtering empty values and enforcing
  the required column names.
* API access is centralised in `ChemblClient` and companion clients that obey
  rate limits, retries and timeouts declared in `config/config.yaml`.
* Normalisation and enrichment happen inside the scripts, relying on helper
  modules such as `document_pipeline`, `target_postprocessing`,
  `testitem_enrichment` and `activity_bounds`. Validation failures are captured
  by `SidecarErrors` and saved next to the export for troubleshooting.
* Deterministic CSV writing plus metadata sidecars are handled by
  `write_csv_deterministic`, `write_meta_yaml`, `add_pipeline_metadata` and the
  table-quality analyser.

## Configuration

* Defaults live in [`config/config.yaml`](../config/config.yaml) and are validated by
  `library.config.load_config` via Pydantic's `Config.model_validate`. The generated
  [`config.schema.json`](../config.schema.json) mirrors the structure for IDEs and tooling.
* Key sections:
  * `sources.*` – base URLs, retry policy, rate limiting and pipeline defaults
    for ChEMBL, UniProt, IUPHAR, PubMed, Semantic Scholar, OpenAlex, CrossRef
    and PubChem.
  * `local.*` – filesystem layout, CSV formatting defaults and reference
    workbooks.
  * `activity_enrichment` / `activity_bounds` – enrichment and bound-derivation
    toggles for the activity pipeline.
  * `testitem_molecule_enrichment` – optional salt/catalogue augmentation for
    test-item exports.
  * `system.*` – logging, global rate limiting, retry configuration and document
    classification weights.
* Overrides follow the precedence `config/config.yaml` < environment variables < CLI
  arguments. Short aliases such as `CHEMBL_DA_RPS` (`sources.chembl.api.rps`)
  and `CHEMBL_DA_OUTDIR` (`local.io.output_dir`) are exposed for convenience.
  The full matrix is documented in `docs/CONFIG_EN.md`.

## External services

* **ChEMBL REST API** supplies activities, assays, targets, documents and
  molecule metadata. Requests are chunked and retried according to the pipeline
  configuration.
* **PubMed, Semantic Scholar, OpenAlex, CrossRef** enrich document exports with
  bibliographic data and DOI coverage.
* **UniProt** provides protein annotations and ID mapping for the target
  pipeline.
* **IUPHAR** contributes receptor classifications from local CSV snapshots.
* **PubChem** augments test items with canonical identifiers and chemical
  descriptors.

## Installation & tooling

1. Create and activate a Python ≥3.11 virtual environment.
2. Install dependencies from the lock file: `pip install -r requirements-lock.txt`.
3. Enable the quality gate: `pre-commit install`.
4. Recommended ad-hoc checks (see Quality Assurance Process — [English](./QA_PROCESS_EN.md) / [Русский](./QA_PROCESS_RU.md) — for the living checklist):
   * `pre-commit run --all-files`
   * `pytest` / `pytest --cov=library --cov=scripts`
   * `ruff check`, `black --check .`, `mypy`

Dependency versions and optional tooling are declared in `pyproject.toml` and
mirrored in `requirements-lock.txt` for reproducible environments.

## Usage highlights

* Every `scripts/get_*_data.py` command accepts standard flags such as
  `--config`, `--print-config`, `--input`, `--final-out`, `--log-level`, `--sep`,
  `--encoding`, `--column` and either `--batch-size` or `--chunk-size`.
* Pipelines expose additional switches (for example, `--timeout`, `--limit`,
  `--dry-run`, sub-commands for documents and targets). CLI options are merged
  back into the configuration via `apply_config_overrides` before execution so
  downstream helpers receive consistent settings.
* Deterministic logging uses JSON lines with `run_id`, `event`, `stage` and
  per-pipeline counters, allowing easy monitoring with `jq` or log collectors.

Detailed command walkthroughs live in `docs/USAGE_EN.md` and
`docs/USAGE_RU.md`.

## Outputs

* Primary CSV exports live under `local.io.output_dir` (default `~/.local/share/chembl-da/output`).
* Every run writes a `<name>.csv.meta.yaml` sidecar with configuration, command
  line, row counts and SHA-256 digests. Validation errors go into
  `<name>_failure_cases.csv`, and table-quality reports are produced for each
  dataset.
* Pipeline-specific extras include document quality JSON files and intermediate
  target exports for the `all` workflow.

Refer to `docs/OUTPUT_EN.md` / `docs/OUTPUT_RU.md` for field-level details and
examples.

## Testing and determinism

* Deterministic CSV writers and metadata hashing are verified by dedicated unit
  tests and utility CLIs (for example,
  `library.utils.cli_tools.check_determinism`).
* Smoke fixtures under `tests/data/` cover activities, targets, documents and
  test items for offline experimentation.
* `python -m library.utils.cli_tools.table_quality_main` and related tools aid
  validation of third-party datasets before ingestion.

This summary should orient newcomers and serve as the high-level entry point to
more detailed reference material in the `docs/` directory.
