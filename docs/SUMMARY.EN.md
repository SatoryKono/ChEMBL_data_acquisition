# Project Summary

This document provides an at-a-glance description of the ChEMBL data acquisition
utilities, their structure, shared services and supporting workflows.

## Repository layout

* `scripts/` – command-line entry points for the activity, assay, document,
  target and test-item pipelines, plus the cached target harness used for
  offline smoke tests.【F:scripts/get_activity_data.py†L1-L1160】【F:scripts/pipeline_targets_main.py†L1-L141】
* `library/` – reusable modules covering API clients, rate limiting,
  normalisation, enrichment, validation, deterministic I/O, logging and
  metadata helpers.【F:library/__init__.py†L1-L50】【F:library/io.py†L1-L236】
* `schemas/` – `pandera` validation schemas and normalisers that keep column
  ordering, data types and canonical values consistent across exports.【F:schemas/__init__.py†L1-L16】
* `dictionary/` & `data/` – local lookup tables, cached API responses and input
  workbooks consumed by the pipelines.【F:config.yaml†L96-L154】
* `docs/` – the reference documentation you are reading.
* `tests/` – unit and integration coverage for configuration overrides,
  enrichment logic, deterministic exports and CLI wiring.【F:tests/test_activity_pipeline.py†L1-L220】

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
  the required column names.【F:library/io.py†L87-L160】
* API access is centralised in `ChemblClient` and companion clients that obey
  rate limits, retries and timeouts declared in `config.yaml`.
* Normalisation and enrichment happen inside the scripts, relying on helper
  modules such as `document_pipeline`, `target_postprocessing`,
  `testitem_enrichment` and `activity_bounds`. Validation failures are captured
  by `SidecarErrors` and saved next to the export for troubleshooting.【F:library/sidecar.py†L1-L154】
* Deterministic CSV writing plus metadata sidecars are handled by
  `write_csv_deterministic`, `write_meta_yaml`, `add_pipeline_metadata` and the
  table-quality analyser.【F:library/csv_utils.py†L451-L603】【F:library/metadata.py†L29-L133】

## Configuration

* Defaults live in [`config.yaml`](../config.yaml) and are validated against
  [`config.schema.json`](../config.schema.json).
* Key sections:
  * `sources.*` – base URLs, retry policy, rate limiting and pipeline defaults
    for ChEMBL, UniProt, IUPHAR, PubMed, Semantic Scholar, OpenAlex, CrossRef
    and PubChem.【F:config.yaml†L11-L258】
  * `local.*` – filesystem layout, CSV formatting defaults and reference
    workbooks.【F:config.yaml†L108-L154】
  * `activity_enrichment` / `activity_bounds` – enrichment and bound-derivation
    toggles for the activity pipeline.【F:config.yaml†L155-L238】
  * `testitem_molecule_enrichment` – optional salt/catalogue augmentation for
    test-item exports.【F:config.yaml†L239-L269】
  * `system.*` – logging, global rate limiting, retry configuration and document
    classification weights.【F:config.yaml†L270-L315】
* Overrides follow the precedence `config.yaml` < environment variables < CLI
  arguments. Short aliases such as `CHEMBL_DA_RPS` (`sources.chembl.api.rps`)
  and `CHEMBL_DA_OUTDIR` (`local.io.output_dir`) are exposed for convenience.
  The full matrix is documented in `docs/CONFIG_EN.md`.

## External services

* **ChEMBL REST API** supplies activities, assays, targets, documents and
  molecule metadata. Requests are chunked and retried according to the pipeline
  configuration.【F:library/chembl_client.py†L1-L286】
* **PubMed, Semantic Scholar, OpenAlex, CrossRef** enrich document exports with
  bibliographic data and DOI coverage.【F:scripts/get_document_data.py†L242-L533】
* **UniProt** provides protein annotations and ID mapping for the target
  pipeline.【F:library/uniprot_library.py†L1-L357】
* **IUPHAR** contributes receptor classifications from local CSV snapshots.【F:library/target_postprocessing.py†L1-L599】
* **PubChem** augments test items with canonical identifiers and chemical
  descriptors.【F:library/testitem_enrichment.py†L17-L216】

## Installation & tooling

1. Create and activate a Python ≥3.12 virtual environment.
2. Install the project with development extras: `pip install .[dev]`.
3. Enable the quality gate: `pre-commit install`.
4. Recommended ad-hoc checks:
   * `pre-commit run --all-files`
   * `pytest` / `pytest --cov=library --cov=scripts`
   * `ruff check`, `black --check .`, `mypy`

Dependency versions and optional tooling are declared in `pyproject.toml` and
`requirements-dev.txt`.

## Usage highlights

* Every `scripts/get_*_data.py` command accepts standard flags such as
  `--config`, `--print-config`, `--input`, `--output`, `--log-level`, `--sep`,
  `--encoding`, `--column` and either `--batch-size` or `--chunk-size`.
* Pipelines expose additional switches (for example, `--timeout`, `--limit`,
  `--dry-run`, sub-commands for documents and targets). CLI options are merged
  back into the configuration via `apply_config_overrides` before execution so
  downstream helpers receive consistent settings.【F:library/cli.py†L1-L322】
* Deterministic logging uses JSON lines with `run_id`, `event`, `stage` and
  per-pipeline counters, allowing easy monitoring with `jq` or log collectors.

Detailed command walkthroughs live in `docs/USAGE_EN.md` and
`docs/USAGE_RU.md`.

## Outputs

* Primary CSV exports live under `local.io.output_dir` (default `data/output`).
* Every run writes a `<name>.csv.meta.yaml` sidecar with configuration, command
  line, row counts and SHA-256 digests. Validation errors go into
  `<name>_failure_cases.csv`, and table-quality reports are produced for each
  dataset.【F:library/metadata.py†L29-L133】【F:library/table_quality.py†L1-L192】
* Pipeline-specific extras include document quality JSON files and intermediate
  target exports for the `all` workflow.

Refer to `docs/OUTPUT_EN.md` / `docs/OUTPUT_RU.md` for field-level details and
examples.

## Testing and determinism

* Deterministic CSV writers and metadata hashing are verified by dedicated unit
  tests and utility CLIs (for example,
  `library.utils.cli_tools.check_determinism`).【F:library/utils/cli_tools/check_determinism.py†L1-L145】
* Smoke fixtures under `tests/data/` cover activities, targets, documents and
  test items for offline experimentation.
* `python -m library.utils.cli_tools.table_quality_main` and related tools aid
  validation of third-party datasets before ingestion.【F:library/utils/cli_tools/table_quality_main.py†L1-L171】

This summary should orient newcomers and serve as the high-level entry point to
more detailed reference material in the `docs/` directory.
